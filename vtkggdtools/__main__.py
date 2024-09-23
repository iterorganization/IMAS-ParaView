"""Support module to run vtkggdtools as a module:

.. code-block:: bash
    :caption: Options to run vtkggdtools CLI interface

    # Run as a module (implemented in vtkggdtools/__main__.py)
    python -m vtkggdtools

    # Run as "program" (see project.scripts in pyproject.toml)
    vtkggdtools
"""

from vtkggdtools.cli import cli

cli()
