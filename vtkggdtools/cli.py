import logging
import sys
from collections import OrderedDict
from pathlib import Path

import click
import imaspy
import imaspy.backends.imas_core.imas_interface
import numpy as np
from imaspy.backends.imas_core.imas_interface import ll_interface
from rich import box, console, traceback
from rich.table import Table

import vtkggdtools
from vtkggdtools.convert import convert_to_xml

logger = logging.getLogger(__name__)


def _excepthook(type_, value, tb):
    logger.debug("Suppressed traceback:", exc_info=(type_, value, tb))
    # Only display the last traceback frame:
    if tb is not None:
        while tb.tb_next:
            tb = tb.tb_next
    rich_tb = traceback.Traceback.from_exception(type_, value, tb, extra_lines=0)
    console.Console(stderr=True).print(rich_tb)


@click.group("vtkggdtools", invoke_without_command=True, no_args_is_help=True)
@click.option("-v", "--version", is_flag=True, help="Show version information")
def cli(version):
    """vtkggdtools command line interface.

    Please use one of the available commands listed below. You can get help for each
    command by executing:

        vtkggdtools <command> --help
    """
    # Limit the traceback to 1 item: avoid scaring CLI users with long traceback prints
    # and let them focus on the actual error message
    sys.excepthook = _excepthook

    if version:
        print_version()


def print_version():
    """Print version information of vtkggdtools."""
    cons = console.Console()
    grid = Table(
        title="vtkggdtools version info", show_header=False, title_style="bold"
    )
    grid.box = box.HORIZONTALS
    if cons.size.width > 120:
        grid.width = 120
    grid.add_row("vtkggdtools version:", vtkggdtools.__version__)
    grid.add_section()
    grid.add_row("IMASPy version:", imaspy.__version__)
    grid.add_section()
    grid.add_row("Default data dictionary version:", imaspy.IDSFactory().dd_version)
    dd_versions = ", ".join(imaspy.dd_zip.dd_xml_versions())
    grid.add_row("Available data dictionary versions:", dd_versions)
    grid.add_section()
    grid.add_row("Access Layer core version:", ll_interface.get_al_version() or "N/A")
    console.Console().print(grid)


@cli.command("ggd2vtk")
@click.argument("uri", type=str)
@click.argument("output", type=Path)
@click.option("--index", "-i", type=str, help="Specify a single index to convert.")
@click.option("--time", "-t", type=str, help="Specify a specific time step in seconds.")
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
def convert_ggd_to_vtk(
    uri,
    output,
    index,
    time,
    all_times,
    format,
):
    """Convert a GGD structure in an IDS to a VTK file.

    \b
    uri         URI of the Data Entry (e.g. "imas:mdsplus?path=testdb").
    output      Name of the output VTK directory.
    """
    uri, ids_name, occurrence = parse_uri(uri)
    click.echo(f"Loading {ids_name} from {uri} with occurrence {occurrence}...")
    entry = imaspy.DBEntry(uri, "r")
    ids = entry.get(ids_name, occurrence=occurrence, autoconvert=False)
    index_list = parse_time_options(ids.time, index, time, all_times)

    click.echo("Converting GGD to a VTK file...")

    # TODO: Add time-dependent VTKHDF conversion
    if format == "xml":
        convert_to_xml(ids, output, index_list)

    elif format == "vtkhdf":
        click.echo("vtkhdf format is not yet implemented.")


def parse_uri(uri):
    """_summary_

    Args:
        uri: _description_

    Returns:
        _description_
    """
    if "#" in uri:
        split_uri = uri.split("#")
        uri_path = split_uri[0]
        fragment = split_uri[1]
        if "/" in fragment:
            raise click.UsageError(
                "It is currently not possible to select an IDS subset for conversion."
                "It is only possible to convert the entire IDS."
            )
        elif ":" in fragment:
            split_fragment = fragment.split(":")
            ids_name = split_fragment[0]
            occurrence = split_fragment[1]
        else:
            ids_name = fragment
            occurrence = 0
    else:
        raise click.UsageError(
            "The IDS must be provided as a fragment to the URI. For example: "
            'uri = "imas:hdf5?path=testdb#edge_profiles"'
        )
    return uri_path, ids_name, occurrence


def parse_time_options(ids_time, index, time, all_times):
    """_summary_

    Args:
        index: _description_
        time: _description_
        all_times: _description_

    Returns:
        _description_
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

    return index_list


def parse_index(index):
    """_summary_

    Args:
        index: _description_

    Returns:
        _description_
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
                "Duplicate time steps were detected. Note that provided time steps "
                "will be rounded down to the nearest found time in the IDS. All "
                "duplicates time steps will be ignored."
            )
            index_list = list(indices_dict)
    # Range of indices
    elif ":" in index:
        if index.count(":") > 1:
            raise click.UsageError("Only a single range may be provided.")
        start_str, end_str = index.split(":")
        if not start_str.strip().isdigit() or not end_str.strip().isdigit():
            raise click.UsageError(
                "The lower and upper bound of indices must be " "integers."
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
    # List of time steps
    if "," in time:
        for input_time in time.split(","):
            try:
                float(input_time)
            except ValueError:
                raise click.UsageError(
                    "All time steps in given list must be valid floats."
                )
        time_list = [float(x.strip()) for x in time.split(",")]
        index_list = find_closest_indices(time_list, ids_times)
        indices_dict = OrderedDict.fromkeys(index_list)
        if len(index_list) != len(indices_dict):
            click.echo(
                "Duplicate time steps were detected. Note that provided time steps "
                "will be rounded down to the nearest found time in the IDS. All "
                "duplicates time steps will be ignored."
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
        except ValueError:
            raise click.UsageError(
                "The minimum and maximum range values must be valid floats."
            )
        if end < start:
            raise click.UsageError(
                "The final time index in range must be greater than the first."
            )
        index_list = [
            index for index, value in enumerate(ids_times) if start <= value <= end
        ]
        if index_list == []:
            raise click.UsageError(
                "Could not find any time steps between in provided range."
            )
    # Single time step
    else:
        try:
            time_list = [float(time)]
            index_list = find_closest_indices(time_list, ids_times)
        except ValueError:
            raise click.UsageError(
                "Could not determine which time steps should be converted.\n"
                "Provide either a single float ('-t 5.0'), "
                "a list of floats ('-t 2.1,3.5,4') or "
                "a range of floats ('-t 2.2:4.4')"
            )
    click.echo(f"Converting the following time steps: {ids_times[index_list]}")
    return index_list


def find_closest_indices(values_to_extract, source_array):
    closest_indices = []
    for value in values_to_extract:
        candidates = source_array[source_array <= value]
        if candidates.size > 0:
            closest_value = candidates.max()
            closest_index = np.where(source_array == closest_value)[0][0]
            closest_indices.append(closest_index)
    return closest_indices


if __name__ == "__main__":
    cli()
