"""Support module to run IMAS-Paraview as a module:

.. code-block:: bash
    :caption: Options to run IMAS-Paraview CLI interface

    # Run as a module (implemented in imas_paraview/__main__.py)
    python -m imas_paraview

    # Run as "program" (see project.scripts in pyproject.toml)
    imas-paraview
"""

from imas_paraview.cli import cli

cli()
