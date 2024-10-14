import logging
import sys
from pathlib import Path

import click
import imaspy
import imaspy.backends.imas_core.imas_interface
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
@click.argument("ids", type=str)
@click.argument("output", type=Path)
@click.argument(
    "occurrence",
    type=int,
    default=0,
)
@click.option("--index", "-i", type=int, help="Specify a single index to convert.")
@click.option(
    "--index-range", "-ir", type=str, help="Specify a range of indices as 'start:end'."
)
@click.option(
    "--time", "-t", type=float, help="Specify a specific time step in seconds."
)
@click.option(
    "--time-range",
    "-tr",
    type=str,
    help="Specify a time range as 'start:end' in seconds.",
)
@click.option(
    "--all-times", "-a", is_flag=True, help="Convert all available time steps."
)
def convert_ggd_to_vtk(
    uri,
    ids,
    output,
    occurrence,
    index,
    index_range,
    time,
    time_range,
    all_times,
    vtk_mode,
):
    """Convert a GGD structure in an IDS to a VTK file.

    \b
    uri         URI of the Data Entry (e.g. "imas:mdsplus?path=testdb").
    ids         Name of the IDS to print (e.g. "edge_profiles").
    output      Name of the output VTK file/directory.
    occurrence  Which occurrence to print (defaults to 0).
    """
    index, index_range, time, time_range, all_times = validate_time_options(
        index, index_range, time, time_range, all_times
    )
    click.echo(f"Loading {ids} from {uri}...")
    entry = imaspy.DBEntry(uri, "r")
    ids = entry.get(ids, occurrence=occurrence, autoconvert=False)

    click.echo("Converting GGD to a VTK file...")

    # TODO: Add time-dependent VTKHDF conversion
    if vtk_mode == "xml":
        convert_to_xml(ids, output, index, index_range, time, time_range, all_times)

    elif vtk_mode == "vtkhdf":
        click.echo("vtkhdf format is not yet implemented.")


def validate_time_options(index, index_range, time, time_range, all_times):
    """_summary_

    Args:
        time: _description_
        time_range: _description_
        all_times: _description_
        time_mode: _description_
    """
    # Check if at most one of the time options are provided
    if (
        sum(param is not None for param in [index, index_range, time, time_range])
        + all_times
        > 1
    ):
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
            start, end = map(float, time_range.split(":"))
            click.echo(f"Converting time steps between {start} and {end}")
            if end < start:
                raise click.UsageError("The final time must be greater than the first.")
            time_range = [start, end]
        except Exception:
            raise click.UsageError(
                "Time range must be in the format 'start:end' with valid numbers."
            )
    elif index_range:
        try:
            start, end = map(int, index_range.split(":"))
            click.echo(f"Converting time indices {list(range(start, end+1))}")
            if end < start:
                raise click.UsageError(
                    "The final time index must be greater than the first."
                )
            index_range = [start, end]
        except Exception:
            raise click.UsageError(
                "Time range must be in the format 'start:end' with valid integers."
            )
    elif time is not None:
        click.echo(f"Converting time {time}")
    elif index is not None:
        click.echo(f"Converting time at index {index}")
    else:
        click.echo("No time options were set, so only converting the first time step.")
        index = 0

    return index, index_range, time, time_range, all_times


if __name__ == "__main__":
    cli()
