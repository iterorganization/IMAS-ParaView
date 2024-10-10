import logging
import sys
from pathlib import Path

import click
import imaspy
import imaspy.backends.imas_core.imas_interface
from imaspy.backends.imas_core.imas_interface import ll_interface
from rich import box, console, traceback
from rich.table import Table
from vtk import vtkXMLPartitionedDataSetCollectionWriter

import vtkggdtools
from vtkggdtools.convert import ggd_to_vtk

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
@click.argument("ids", type=str)
@click.argument("output", type=str)
@click.argument(
    "occurrence",
    type=int,
    default=0,
)
@click.option(
    "--time-mode",
    type=click.Choice(["index", "value"], case_sensitive=True),
    default="index",
    help="Select time mode: 'index' (default) for slicing by index or 'value' for "
    "slicing by time value.",
)
@click.option("--time", type=float, help="Convert single time slice.")
@click.option("--time-range", type=str, help="Convert time range in format start:end.")
@click.option("--all-times", is_flag=True, help="Convert all available times.")
@click.option(
    "--vtk-mode",
    type=click.Choice(["xml", "vtkhdf"], case_sensitive=True),
    default="xml",
    help="VTK output format: 'xml' (default) for standard VTK XML files, "
    "or 'vtkhdf' for VTK files using HDF5 format.",
)
def convert_ggd_to_vtk(
    uri, ids, output, occurrence, time_mode, time, time_range, all_times, vtk_mode
):
    """Convert a GGD structure in an IDS to a VTK file.

    \b
    uri         URI of the Data Entry (e.g. "imas:mdsplus?path=testdb").
    ids         Name of the IDS to print (e.g. "edge_profiles").
    output      Name of the output VTK file/directory.
    occurrence  Which occurrence to print (defaults to 0).
    """
    time_options = validate_time_options(time, time_range, all_times, time_mode)

    click.echo(f"Loading {ids} from {uri}...")
    entry = imaspy.DBEntry(uri, "r")
    ids = entry.get(ids, occurrence=occurrence, autoconvert=False)

    click.echo("Converting GGD to a VTK file...")

    # TODO: Add time-dependent VTKHDF conversion
    if vtk_mode == "xml":
        convert_to_xml(ids, output, time_options)

    elif vtk_mode == "vtkhdf":
        click.echo("vtkhdf format is not yet implemented.")


def convert_to_xml(ids, output, time_options):
    (time, time_range, all_times, time_mode) = time_options
    vtk_object = None
    if time is not None:
        if time_mode == "index":
            click.echo(f"Converting time step closest to {time}")
            vtk_object = ggd_to_vtk(ids, time=time)
        elif time_mode == "value":
            click.echo(f"Converting time step at index {time}")
            vtk_object = ggd_to_vtk(ids, time_idx=time)
        write_vtk(vtk_object, output)


def write_vtk(vtk_object, output):
    if vtk_object is None:
        click.echo("Could not convert GGD to VTK file.")
        raise RuntimeError
    click.echo("Writing VTK file to disk...")
    writer = vtkXMLPartitionedDataSetCollectionWriter()
    writer.SetInputData(vtk_object)
    output_file = Path(output).with_suffix(".vtpc")
    writer.SetFileName(output_file)
    writer.Write()
    click.echo(f"Successfully wrote VTK object to {output_file}.")


def validate_time_options(time, time_range, all_times, time_mode):
    """_summary_

    Args:
        time: _description_
        time_range: _description_
        all_times: _description_
        time_mode: _description_
    """

    if sum([time is not None, time_range is not None, all_times]) > 1:
        raise click.UsageError(
            "You must provide only one of --time, --time-range, or --all-times."
        )

    if all_times:
        click.echo(
            "Converting all time steps in the IDS. Depending on the number of time "
            "steps, this could take a while."
        )
    elif time_range:
        try:
            if time_mode == "index":
                start, end = map(int, time_range.split(":"))
                click.echo(f"Converting time indices {list(range(start, end+1))}")
            elif time_mode == "value":
                start, end = map(float, time_range.split(":"))
                click.echo(
                    f"Converting time steps in the range {list(range(start, end+1))}"
                )
            if end < start:
                raise click.UsageError(
                    "The final time index must be greater than or equal to the first."
                )
        except Exception:
            raise click.UsageError(
                "Time range must be in the format 'start:end' with valid numbers."
            )

    elif time is not None:
        if time_mode == "index":
            if not time.is_integer():
                raise click.UsageError(
                    "The time-mode is set to 'index', so an integer must be provided."
                )
            time = int(time)
        click.echo(f"Converting time {time}")
    else:
        click.echo("No time options were set, so only converting the first time step.")
        time_mode = "index"
        time = 0

    return (time, time_range, all_times, time_mode)


if __name__ == "__main__":
    cli()
