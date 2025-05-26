.. _`installing`:

Installing IMAS-ParaView
========================

If the IMAS-ParaView module is available on your HPC environment, you can ignore the further 
steps below and simply load the module as follows:

# TODO: Change module name when new IMAS-ParaView easybuild module is available
.. code-block:: bash

  module load GGD-VTK
  # Launch Paraview
  paraview
  # Open the "Sources" tab in the top left, and ensure you see "IMAS Tools" in 
  # the drop down menu.
  # The command-line interface can be started using:
  imas-paraview --help


SDCC installation
-----------------

* Setup a project folder and clone git repository

.. code-block:: bash

  mkdir projects
  cd projects
  git clone git@github.com:iterorganization/IMAS-ParaView.git
  cd IMAS-ParaView


* To run a plugin in Paraview, run the following at the root of the project directory.

.. code-block:: bash

  # Load compatible IMAS-Python, IMAS-Core and ParaView modules, like:
  module load IMAS-AL-Core/5.4.3-foss-2023b IMAS-Python/2.0.0-foss-2023b \
  ParaView/5.12.0-foss-2023b
  # export environment variables, this assumes the current
  # working directory is the root of the repository
  export PV_PLUGIN_PATH=$PWD/imas_paraview/plugins:$PV_PLUGIN_PATH
  export PYTHONPATH=$PWD:$PYTHONPATH
  # Run paraview (add vglrun to enable hardware acceleration)
  vglrun paraview
  # Open the "Sources" tab in the top left, if you see "IMAS Tools" 
  # in the drop down, it is installed correctly.

* To use the command-line interface, setup a python virtual environment and install python dependencies

.. code-block:: bash

  # Load compatible IMAS-Python, IMAS-Core and ParaView modules, like:
  module load IMAS-AL-Core/5.4.3-foss-2023b IMAS-Python/2.0.0-foss-2023b \
  ParaView/5.12.0-foss-2023b
  # create virtual environment and install dependencies
  python3 -m venv ./venv
  . venv/bin/activate
  pip install --upgrade pip
  pip install --upgrade wheel setuptools
  # For development install in editable mode
  pip install -e .[all]
  # Run CLI with help information
  imas-paraview --help
  # If you see the help page of IMAS-ParaView, it is installed correctly.

* Every time that a new session is started, ensure the correct modules are loaded, 
  the python virtual environment is activated, and the environment variables are set.

.. code-block:: bash

  # Load the required modules
  module load IMAS-AL-Core/5.4.3-foss-2023b IMAS-Python/2.0.0-foss-2023b \
  ParaView/5.12.0-foss-2023b
  # Export the environment variables
  export PV_PLUGIN_PATH=$PWD/imas_paraview/plugins:$PV_PLUGIN_PATH
  export PYTHONPATH=$PWD:$PYTHONPATH
  # And activate the Python virtual environment
  . venv/bin/activate
  # Validate if it is working as intended
  imas-paraview --version

* To run the unit and integration tests, make sure the install is working using the 
  code block above. Also ensure the optional test dependencies are pip installed (or 
  simply use all, to install all optional dependencies).

.. code-block:: bash

  # The integration tests require X virtual framebuffer to be installed
  module load Xvfb/21.1.9-GCCcore-13.2.0
  python -m pytest
  # Alternatively, if you want to skip running the integration tests
  python -m pytest -m "not integration"

* To build the IMAS-ParaView documentation, ensure the optional docs dependencies are pip 
  installed (or simply use all, to install all optional dependencies).

.. code-block:: bash

  make -C docs html
  # You can now open ./docs/_build/html/index.html

..
  TODO: add local installing documentation, maybe wait until ggd-vtk goes open source?
  As it needs to be installed with IMAS-Python.
