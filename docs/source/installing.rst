.. _`installing`:

Installing GGD-VTK
==================

SDCC installation
-----------------

* Setup a project folder and clone git repository

  .. code-block:: bash

    mkdir projects
    cd projects
    git clone ssh://git@git.iter.org/imex/ggd-vtk.git
    cd ggd-vtk


* To run a plugin in Paraview

  .. code-block:: bash

    # Load compatible IMASPy, IMAS and ParaView modules, like:
    # AL5 and ParaView 5.12 (recommended on RHEL9):
    module load IMASPy/1.1.0-foss-2023b IMAS-AL-Python/5.3.0-foss-2023b-DD-3.42.0 ParaView/5.12.0-foss-2023b
    # export environment variables, assumes the current working directory is the root of the repository
    export PV_PLUGIN_PATH=$PWD/vtkggdtools/plugins:$PV_PLUGIN_PATH PYTHONPATH=$PWD:$PYTHONPATH
    # Use LD_PRELOAD to work around a VTK bug: https://gitlab.kitware.com/vtk/vtk/-/issues/19373
    export LD_PRELOAD=$HDF5_DIR/lib/libhdf5.so.310
    # Run paraview
    paraview
    # You should now see plugins appear in Paraview under 'Sources/VTKGGDTools'

* To use the command-line interface, setup a python virtual environment and install python dependencies

  .. code-block:: bash

    # Load IMAS and IMASPy before install
    module load IMAS IMASPy
    cd projects/ggd-vtk
    python3 -m venv ./venv
    . venv/bin/activate
    pip install --upgrade pip
    pip install --upgrade wheel setuptools
    # For development install in editable mode
    pip install -e .[all]
    # Run CLI with help information
    vtkggdtools --help

* Run each session if tool is already installed

  .. code-block:: bash

    # Load modules every time you use vtkggdtools
    module load IMAS IMASPy
    # And activate the Python virtual environment
    . venv/bin/activate

* Test the installation

  .. code-block:: bash

    vtkggdtools --version
    python -m pytest

* To build the GGD-VTK documentation, execute:

  .. code-block:: bash

    make -C docs html
