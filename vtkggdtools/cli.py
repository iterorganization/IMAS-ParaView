import logging
import sys

import click
import imaspy
import imaspy.backends.imas_core.imas_interface
from imaspy.backends.imas_core.imas_interface import ll_interface
from rich import box, console, traceback
from rich.table import Table

import vtkggdtools

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
def cli():
    """vtkggdtools command line interface.

    Please use one of the available commands listed below. You can get help for each
    command by executing:

        vtkggdtools <command> --help
    """
    # Limit the traceback to 1 item: avoid scaring CLI users with long traceback prints
    # and let them focus on the actual error message
    sys.excepthook = _excepthook


@cli.command("version")
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


@cli.command("convert_ggd_to_vtk")
@click.argument("uri", type=str)
@click.argument("ids", type=str)
@click.argument("output_name", type=str)
@click.argument("occurrence", type=int, default=0)
def convert_ggd_to_vtk(uri, ids, output_name, occurrence):
    """Convert a GGD structure to a VTK file.

    \b
    uri         URI of the Data Entry (e.g. "imas:mdsplus?path=testdb").
    ids         Name of the IDS to print (e.g. "core_profiles").
    output      Filename of the output vtk file.
    occurrence  Which occurrence to print (defaults to 0).
    """

    click.echo(
        f"Converting IDS {ids} from uri {uri} with occurrence {occurrence} "
        "to a VTK file."
    )

    # TODO add conversion function for ggd to vtk

    click.echo(f"Successfully converted to {output_name}.")


if __name__ == "__main__":
    cli()
