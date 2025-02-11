.. _`installing`:

Installing GGD-VTK
==================

If the GGD-VTK module is available on your HPC environment, you can ignore the further 
steps below and simply load the module as follows:

.. code-block:: bash

  module load GGD-VTK
  paraview


SDCC installation
-----------------

* Setup a project folder and clone git repository

.. code-block:: bash

  mkdir projects
  cd projects
  git clone ssh://git@git.iter.org/imex/ggd-vtk.git
  cd ggd-vtk


* To run a plugin in Paraview, run the following at the root of the project directory.

.. code-block:: bash

  # Load compatible IMASPy, IMAS and ParaView modules, like:
  # AL5 and ParaView 5.12 (recommended on RHEL9):
  module load IMASPy/1.1.1-foss-2023b \
  IMAS-AL-Python/5.3.0-foss-2023b-DD-3.42.0 ParaView/5.12.0-foss-2023b
  # export environment variables, this assumes the current
  # working directory is the root of the repository
  export PV_PLUGIN_PATH=$PWD/vtkggdtools/plugins:$PV_PLUGIN_PATH
  export PYTHONPATH=$PWD:$PYTHONPATH
  # Use LD_PRELOAD to work around a VTK bug:
  # https://gitlab.kitware.com/vtk/vtk/-/issues/19373
  export LD_PRELOAD=$HDF5_DIR/lib/libhdf5.so.310
  # Run paraview
  paraview
  # Open the "Sources" tab in the top left, if you see "IMAS Tools" 
  # in the drop down, it is installed correctly.

* To use the command-line interface, setup a python virtual environment and install python dependencies

.. code-block:: bash

  # Load compatible IMASPy, IMAS and ParaView modules, like:
  module load IMASPy/1.1.1-foss-2023b \
  IMAS-AL-Python/5.3.0-foss-2023b-DD-3.42.0 ParaView/5.12.0-foss-2023b
  # create virtual environment and install dependencies
  python3 -m venv ./venv
  . venv/bin/activate
  pip install --upgrade pip
  pip install --upgrade wheel setuptools
  # For development install in editable mode
  pip install -e .[all]
  # Run CLI with help information
  vtkggdtools --help
  # If you see the help page of the vtkggdtools, it is installed correctly.

* Every time that a new session is started, ensure the correct modules are loaded, 
  the python virtual environment is activated, and the environment variables are set.

.. code-block:: bash

  # Load the required modules
  module load IMASPy/1.1.1-foss-2023b \
  IMAS-AL-Python/5.3.0-foss-2023b-DD-3.42.0 ParaView/5.12.0-foss-2023b
  # Export the environment variables
  export PV_PLUGIN_PATH=$PWD/vtkggdtools/plugins:$PV_PLUGIN_PATH
  export PYTHONPATH=$PWD:$PYTHONPATH
  export LD_PRELOAD=$HDF5_DIR/lib/libhdf5.so.310
  # And activate the Python virtual environment
  . venv/bin/activate
  # Validate if it is working as intended
  vtkggdtools --version

* To run the unit and integration tests, make sure the install is working using the 
  code block above. Also ensure the optional test dependencies are pip installed (or 
  simply use all, to install all optional dependencies).

.. code-block:: bash

  # The integration tests require X virtual framebuffer to be installed
  module load Xvfb/21.1.9-GCCcore-13.2.0
  python -m pytest
  # Alternatively, if you want to skip running the integration tests
  python -m pytest -m "not integration"

* To build the GGD-VTK documentation, ensure the optional docs dependencies are pip 
  installed (or simply use all, to install all optional dependencies).

.. code-block:: bash

  make -C docs html
  # You can now open ./docs/_build/html/index.html

..
  TODO: add local installing documentation, maybe wait until ggd-vtk goes open source?
  As it needs to be installed with IMASPy.
