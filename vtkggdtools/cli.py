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
@click.argument("time", type=int, default=None)
@click.argument("occurrence", type=int, default=0)
def convert_ggd_to_vtk(uri, ids, output, time, occurrence):
    """Convert a GGD structure in an IDS to a VTK file.

    \b
    uri         URI of the Data Entry (e.g. "imas:mdsplus?path=testdb").
    ids         Name of the IDS to print (e.g. "edge_profiles").
    output      Name of the output VTK file/directory.
    time        Time slice to convert (defaults to first time step).
    occurrence  Which occurrence to print (defaults to 0).
    """

    entry = imaspy.DBEntry(uri, "r")
    click.echo(f"Loading {ids} from {uri}...")
    ids = entry.get(ids, occurrence=occurrence, autoconvert=False)

    click.echo("Converting GGD to a VTK file...")
    # TODO: convert timesteps other than first
    vtk_object = ggd_to_vtk(ids, time=time)
    if vtk_object is None:
        click.echo("Could not convert GGD to VTK file.")
        return

    click.echo("Writing VTK file to disk...")
    writer = vtkXMLPartitionedDataSetCollectionWriter()
    writer.SetInputData(vtk_object)
    output_file = Path(output).with_suffix(".vtpc")
    writer.SetFileName(output_file)
    writer.Write()

    click.echo(f"Successfully converted GGD to {output_file}.")


if __name__ == "__main__":
    cli()
