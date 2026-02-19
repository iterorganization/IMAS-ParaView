import logging
import sys
from collections import OrderedDict
from pathlib import Path

import click
import imas
import imas.backends.imas_core.imas_interface
from imas.backends.imas_core.imas_interface import ll_interface
from rich import box, console, traceback
from rich.table import Table

import imas_paraview
from imas_paraview.convert import Converter
from imas_paraview.util import find_closest_indices, load_vtpc
from imas_paraview.vtk2ggd import VTK2GGDConverter

logger = logging.getLogger(__name__)


def _excepthook(type_, value, tb):
    logger.debug("Suppressed traceback:", exc_info=(type_, value, tb))
    # Only display the last traceback frame:
    if tb is not None:
        while tb.tb_next:
            tb = tb.tb_next
    rich_tb = traceback.Traceback.from_exception(type_, value, tb, extra_lines=0)
    console.Console(stderr=True).print(rich_tb)


@click.group("imas-paraview", invoke_without_command=True, no_args_is_help=True)
@click.option("-v", "--version", is_flag=True, help="Show version information")
def cli(version):
    """IMAS-ParaView command line interface.

    Please use one of the available commands listed below. You can get help for each
    command by executing:

        imas-paraview <command> --help
    """
    # Limit the traceback to 1 item: avoid scaring CLI users with long traceback prints
    # and let them focus on the actual error message
    sys.excepthook = _excepthook

    if version:
        print_version()


def print_version():
    """Print version information of IMAS-ParaView."""
    cons = console.Console()
    grid = Table(
        title="IMAS-ParaView version info", show_header=False, title_style="bold"
    )
    grid.box = box.HORIZONTALS
    if cons.size.width > 120:
        grid.width = 120
    grid.add_row("IMAS-ParaView version:", imas_paraview.__version__)
    grid.add_section()
    grid.add_row("IMAS-Python version:", imas.__version__)
    grid.add_section()
    grid.add_row("Default data dictionary version:", imas.IDSFactory().dd_version)
    dd_versions = ", ".join(imas.dd_zip.dd_xml_versions())
    grid.add_row("Available data dictionary versions:", dd_versions)
    grid.add_section()
    grid.add_row("Access Layer core version:", ll_interface.get_al_version() or "N/A")
    console.Console().print(grid)


@cli.command("ggd2vtk")
@click.argument("uri", type=str)
@click.argument("output_dir", type=Path)
@click.option(
    "--index",
    "-i",
    type=str,
    help="Specify an index, a list of indices or a range of indices to convert.",
)
@click.option(
    "--time",
    "-t",
    type=str,
    help="Specify a time step, a list of time steps or a time step range in seconds to "
    "convert.",
)
@click.option(
    "--all-times", "-a", is_flag=True, help="Convert all available time steps."
)
@click.option(
    "--format",
    "-f",
    type=click.Choice(["xml", "vtkhdf"], case_sensitive=True),
    default="xml",
    help="VTK output format: 'xml' (default) for standard VTK XML files, "
    "or 'vtkhdf' for VTK files using HDF5 format.",
)
@click.option("--lazy", is_flag=True, help="Enable lazy loading.")
@click.option("--no-lazy", is_flag=True, help="Disable lazy loading.")
def convert_ggd_to_vtk(
    uri,
    output_dir,
    index,
    time,
    all_times,
    format,
    lazy,
    no_lazy,
):
    """
    Convert a GGD structure in an IDS and write the converted VTK file to disk.
    Optionally, specific time values or time indices can be specified to be converted.

    \b
    Arguments:
    \b
    uri             URI of the Data Entry. This should contain the IDS and its fragment.
    output_dir      Output directory to store the VTK files.

    ------------------------ Examples ------------------------

    Without time slicing:

     \b
     To convert the middle time step (default), no time options need to be given:
         $ ggd2vtk imas:hdf5?path=testdb#edge_profiles test_dir
     \b
     To convert all time slices in the IDS:
         $ ggd2vtk imas:hdf5?path=testdb#edge_profiles test_dir -a
     \b
     To convert occurrence number 1 (default=0):
         $ ggd2vtk imas:hdf5?path=testdb#edge_profiles:1 test_dir

    Index-based slicing:

     \b
     To convert index 5:
         $ ggd2vtk imas:hdf5?path=testdb#edge_profiles test_dir -i 5
     \b
     To convert indices 5, 8, and 9:
         $ ggd2vtk imas:hdf5?path=testdb#edge_profiles test_dir -i 5,8,9
     \b
     To convert a range of indices, such as 2,3,4,5,6:
         $ ggd2vtk imas:hdf5?path=testdb#edge_profiles test_dir -i 2:6

    Time-based slicing:

     \b
     To convert time step at 5.5s:
         $ ggd2vtk imas:hdf5?path=testdb#edge_profiles test_dir -t 5.5
     \b
     To convert time steps 5.5s, 8s, and 9.1s:
         $ ggd2vtk imas:hdf5?path=testdb#edge_profiles test_dir -t 5.5,8,9.1
     \b
     To convert all time steps that fall between 2.2s and 6.6s:
         $ ggd2vtk imas:hdf5?path=testdb#edge_profiles test_dir -t 2.2:6.6

        \b
        Note: If the specified time step is not found in the IDS, the closest earlier
        time step will be used by default.
    """

    sys.excepthook = _excepthook
    uri, ids_name, occurrence = parse_uri(uri)
    click.echo(f"Loading {ids_name} from {uri} with occurrence {occurrence}...")
    with imas.DBEntry(uri, "r") as entry:
        click.echo("Loading IDS...")
        use_lazy = is_lazy(all_times, lazy, no_lazy)
        ids = entry.get(
            ids_name,
            lazy=use_lazy,
            occurrence=occurrence,
            autoconvert=False,
        )
        index_list = parse_time_options(ids.time, index, time, all_times)

        click.echo("Converting GGD to a VTK file...")

        # TODO: Add time-dependent VTKHDF conversion
        if format == "xml":
            converter = Converter(ids, dbentry=entry)
            converter.write_to_xml(output_dir, index_list)

        elif format == "vtkhdf":
            raise NotImplementedError("vtkhdf format is not yet implemented.")


@cli.command("vtk2ggd")
@click.argument("path", type=Path)
@click.argument("uri", type=str)
@click.option(
    "--ids_name",
    type=str,
    help="Name of the IDS, if not provided, it will be inferred from filenames.",
)
@click.option(
    "--dd_version",
    type=str,
    help="Specify the version of the Data Dictionary to export the IDS to.",
)
@click.option(
    "--cylindrical_coordinates",
    is_flag=True,
    help="Store GGD grid in cylindrical coordinates instead of cartesian.",
)
def convert_vtk_to_ggd(path, uri, ids_name, dd_version, cylindrical_coordinates):
    """Convert a file series of .vtpc files, or a single .vtpc containing
    unstructured grids to a GGD grid in an IDS.

    \b
    Arguments:
    \b
    path            Path to a single .vtpc file, or to the directory containing a series
                    of .vtpc files.
    uri             Output URI to store the IDS.
    """
    sys.excepthook = _excepthook

    vtpc_files = []

    if path.is_dir():  # The input is a file series
        click.echo(f"Scanning directory {path} for VTK files...")
        # Find all .vtpc files and sort them by the trailing index (e.g., _0, _1)
        string_match = f"{ids_name}_*.vtpc" if ids_name else "*.vtpc"

        vtpc_files = sorted(
            path.glob(string_match),
            key=lambda x: int(x.stem.split("_")[-1]) if "_" in x.stem else 0,
        )
        if not vtpc_files:
            raise click.UsageError(f"No .vtpc files found in {path}")

        # Infer ids_name from the file names (edge_profiles_0.vtpc -> edge_profiles)
        if not ids_name:
            ids_names = set()
            for vtpc_file in vtpc_files:
                ids_names.add(vtpc_file.stem.rsplit("_", 1)[0])
                if len(ids_names) > 1:
                    raise click.UsageError(
                        "Could not infer IDS name from file names, because there are "
                        "multiple IDS names in the directory. Use '--ids_name' to set "
                        "the required IDS name"
                    )

            (ids_name,) = ids_names
    else:
        if not ids_name:
            raise click.UsageError(
                "Argument '--ids_name' is required for single file conversion."
            )
        vtpc_files = [path]

    if ids_name not in imas.IDSFactory().ids_names():
        raise click.UsageError(f"{ids_name} is not a valid IDS name")

    click.echo(f"Converting {len(vtpc_files)} time steps for IDS '{ids_name}'...")

    vtk_objects = [load_vtpc(f) for f in vtpc_files]
    converter = VTK2GGDConverter(
        vtk_objects,
        ids_name,
        dd_version=dd_version,
        cylindrical_coordinates=cylindrical_coordinates,
    )
    ids = converter.convert()

    with imas.DBEntry(uri, "x", dd_version=dd_version) as entry:
        entry.put(ids)

    click.echo(f"Successfully wrote {len(vtpc_files)} time steps to {uri}")


def is_lazy(all_times, lazy_flag, no_lazy_flag):
    """Determine whether to enable or disable lazy loading, based on
    if the all_times flag is enabled. Optionally, the default behaviour can be
    overridden by using the --lazy or --no-lazy flag.
    Args:
        all_times: Flag to convert all time steps in an IDS.
        lazy_flag: Flag to enable lazy loading.
        no_lazy_flag: Flag to disable lazy loading.
    """
    if lazy_flag and no_lazy_flag:
        click.echo("Both --lazy and --no-lazy flag were provided. Ignoring...")
        lazy_flag = None
        no_lazy_flag = None

    # Lazy loading is enabled by default
    use_lazy_loading = True
    if all_times:
        use_lazy_loading = False

    # Override default if flag is provided
    if lazy_flag:
        use_lazy_loading = True
    elif no_lazy_flag:
        use_lazy_loading = False

    if use_lazy_loading:
        click.echo("Lazy loading is enabled.")
    else:
        click.echo("Lazy loading is disabled.")
    return use_lazy_loading


def parse_uri(uri):
    """Parses the URI according to the IMAS URI Scheme Documentation v0.3.
    Parses the URI and extracts the fragment part of the URI.

    Example:
        >>> uri = "imas:hdf5?path=testdb#edge_profiles:1"
        >>> parse_uri(uri)
        ("imas:hdf5?path=testdb", "edge_profiles", 1)

    Args:
        uri: URI to parse, should contain a fragment denoting the IDS name.

    Returns:
        A tuple with three elements.

        - uri_no_frag: URI without fragment part
        - ids_name: Name of the IDS.
        - occurrence: Occurrence number of the IDS.
    """
    if "#" in uri:
        split_uri = uri.split("#")
        uri_no_frag = split_uri[0]
        fragment = split_uri[1]
        if "/" in fragment:
            raise click.UsageError(
                "It is currently not possible to select an IDS subset for conversion."
                "It is only possible to convert the entire IDS."
            )
        elif ":" in fragment:
            split_fragment = fragment.split(":")
            ids_name = split_fragment[0]
            occurrence = int(split_fragment[1])
        else:
            ids_name = fragment
            occurrence = 0
    else:
        raise click.UsageError(
            "The IDS must be provided as a fragment to the URI. For example: "
            'uri = "imas:hdf5?path=testdb#edge_profiles"'
        )
    return uri_no_frag, ids_name, occurrence


def parse_time_options(ids_time, index, time, all_times):
    """Parses the provided time options. Only a single time option (index, time or
    all_times) may be used. Either a single value ("5"), a list of values ("2,3,4") or
    a range of values ("2:4") may be provided for the index or time parameters.

    Args:
        index: String containing either a single index, a list of indices or a range of
            indices. Alternatively, it can be None.
        time: String containing either a single time step, a list of time steps or a
            range of time steps. Alternatively, it can be None.
        all_times: Boolean to convert all time steps

    Returns:
        A list containing the time indices to convert.
    """
    # Check if more than a single time option is provided
    if sum([index is not None, time is not None, all_times]) > 1:
        raise click.UsageError(
            "You can only provide one time-related argument: either --time, --index, "
            "or --all-times."
        )
    if index:
        index_list = parse_index(index)
    elif time:
        index_list = parse_time(ids_time, time)
    elif all_times:
        index_list = list(range(len(ids_time)))
        click.echo(
            "Converting all time steps in the IDS. Depending on the number of time "
            f"steps, this could take a while. Converting a total of {len(index_list)} "
            "time steps."
        )
    else:
        index = len(ids_time) // 2
        middle_time = ids_time[index]
        click.echo(
            "No time options were set, so only converting the middle time step: "
            f"t = {middle_time}  at index {index}"
        )
        index_list = [index]

    if index_list == []:
        raise click.UsageError(
            "Could not find any time steps for the provided time steps."
        )
    return index_list


def parse_index(index):
    """Parses the index parameter and returns a list of IDS time indices to convert.
    Either a single integer ("5"), a list of integers ("2,3,4") or a range of integers
    ("2:4") may be provided.

    Args:
        index: String containing either a single index, a list of indices or a range of
            indices.

    Returns:
        A list containing the time indices to convert.
    """
    # Single index
    if index.isdigit():
        index_list = [int(index)]
    # List of indices
    elif "," in index:
        for input_index in index.split(","):
            if not input_index.strip().isdigit():
                raise click.UsageError("All indices in given list must be integers.")
        index_list = [int(x.strip()) for x in index.split(",")]
        indices_dict = OrderedDict.fromkeys(index_list)
        if len(index_list) != len(indices_dict):
            click.echo(
                "Duplicate time step indices were detected. Any duplicate indices will "
                "be ignored."
            )

            index_list = list(indices_dict)
    # Range of indices
    elif ":" in index:
        if index.count(":") > 1:
            raise click.UsageError("Only a single range may be provided.")
        start_str, end_str = index.split(":")
        if not start_str.strip().isdigit() or not end_str.strip().isdigit():
            raise click.UsageError(
                "The lower and upper bound of indices must be integers."
            )
        start = int(start_str)
        end = int(end_str)
        if end < start:
            raise click.UsageError(
                "The final time index in range must be greater than the first."
            )
        index_list = list(range(start, end + 1))
    else:
        raise click.UsageError(
            "Could not determine which indices should be converted.\n"
            "Provide either a single integer ('-i 5'), "
            "a list of integers ('-i 2,3,4') or "
            "a range of integers ('-i 2:4')"
        )
    click.echo(f"Converting the following indices: {index_list}")
    return index_list


def parse_time(ids_times, time):
    """Parses the time parameter and returns a list of IDS time indices to convert.
    Either a single float ("5.0"), a list of integers ("2.0,3,4.4") or a range of
    floats ("2.3:4.3") may be provided. If the specified time step is not found in the
    IDS, the closest earlier time step will be converted instead.

    Args:
        time: String containing either a single time step, a list of time steps or a
            range of time steps.

    Returns:
        A list containing the time indices to convert.
    """
    # List of time steps
    if "," in time:
        for input_time in time.split(","):
            try:
                float(input_time)
            except ValueError as err:
                raise click.UsageError(
                    "All time steps in given list must be valid floats."
                ) from err
        time_list = [float(x.strip()) for x in time.split(",")]
        index_list = find_closest_indices(time_list, ids_times)
        indices_dict = OrderedDict.fromkeys(index_list)
        if len(index_list) != len(indices_dict):
            click.echo(
                "Duplicate time steps were detected. Note that provided time steps "
                "will be rounded down to a time in the IDS, and the duplicate time "
                "steps will be ignored. For example, if the IDS contains times "
                "[1.0, 2.0] and the provided time steps are t=1.1,1.8 both are rounded "
                "down to time step 1.0, and as a result only the time step at t=1.0 is "
                "converted."
            )
            index_list = list(indices_dict)
    # Range of time steps
    elif ":" in time:
        if time.count(":") > 1:
            raise click.UsageError("Only a single range may be provided.")
        start_str, end_str = time.split(":")
        try:
            start = float(start_str)
            end = float(end_str)
        except ValueError as err:
            raise click.UsageError(
                "The minimum and maximum range values must be valid floats."
            ) from err
        if end < start:
            raise click.UsageError(
                "The final time index in range must be greater than the first."
            )
        index_list = [
            index for index, value in enumerate(ids_times) if start <= value <= end
        ]
    # Single time step
    else:
        try:
            time_list = [float(time)]
            index_list = find_closest_indices(time_list, ids_times)
        except ValueError as err:
            raise click.UsageError(
                "Could not determine which time steps should be converted.\n"
                "Provide either a single float ('-t 5.0'), "
                "a list of floats ('-t 2.1,3.5,4') or "
                "a range of floats ('-t 2.2:4.4')"
            ) from err
    click.echo(f"Converting the following time steps: {ids_times[index_list]}")
    return index_list


if __name__ == "__main__":
    cli()
