.. _`installing`:

Installing IMAS-ParaView
========================

If the IMAS-ParaView module is available on your HPC environment, you can ignore the further 
steps below and simply load the module as follows:

.. code-block:: bash

  module load IMAS-ParaView
  # Launch Paraview
  paraview
  # Open the "Sources" tab in the top left, and ensure you see "IMAS Tools" in 
  # the drop down menu.
  # The command-line interface can be started using:
  imas-paraview --help


Installation on Ubuntu 24.04
----------------------------

.. caution::
  The paraview package in the Ubuntu 24.04 apt repository has compatibility issues with
  Python 3.12 (which is the Python version distributed with Ubuntu 24.04).

  Instead, we will download and use the program from the ParaView website.

1.  Download ParaView 5.13.3 from the Paraview website:
    https://www.paraview.org/download/?version=v5.13&filter=Linux
2.  Unzip the paraview package:

    .. code-block:: bash
    
      tar -xzvf ParaView-5.13.3-MPI-Linux-Python3.10-x86_64.tar.gz

3.  Install ``uv``:

    .. code-block:: bash
      
      python -m venv uvvenv
      . uvvenv/bin/activate
      pip install uv

4.  Use ``uv`` to create a virtual environment with Python 3.10 (which is the version
    used by ParaView) and install the required packages.

    Note: you will need the ``imas_core`` package to load IMAS data that is stored in
    the HDF5 or MDSplus backends. Unfortunately, this component is not publicly
    available, but can be installed provided you have access to its git repository.
    See
    https://sharepoint.iter.org/departments/POP/CM/IMDesign/Code%20Documentation/ACCESS-LAYER-doc/python/dev/building_installing.html#prerequisites
    (behind login wall) for the prerequisites.

    .. code-block:: bash

      uv venv --python 3.10
      uv pip install --python 3.10 imas-paraview
      # Uninstall VTK, this would conflict with ParaView's built-in VTK module
      uv pip uninstall --python 3.10 vtk
      # Optional: install imas-core package from the ITER git:
      git clone ssh://git@git.iter.org/imas/al-core.git -b main
      uv pip install --python 3.10 ./al-core

5.  Set required environment variables and run ParaView:

    .. code-block:: bash

      export PYTHONPATH="$PWD/.venv/lib/python3.10/site-packages/"
      export PV_PLUGIN_PATH="$PWD/.venv/lib/python3.10/site-packages/imas_paraview/plugins"
      # This LD_PRELOAD is required when loading data with the imas_core HDF5 backend.
      # Note that it may break any built-in ParaView HDF5 support...
      export LD_PRELOAD="/lib/x86_64-linux-gnu/libhdf5_serial.so.103"
      ./ParaView-5.13.3-MPI-Linux-Python3.10-x86_64/bin/paraview


Installation on Ubuntu 22.04
----------------------------

1.  Download ParaView 5.13.3 from the Paraview website:
    https://www.paraview.org/download/?version=v5.13&filter=Linux

2.  Unzip the paraview package:

    .. code-block:: bash
    
      tar -xzvf ParaView-5.13.3-MPI-Linux-Python3.10-x86_64.tar.gz

3.  Create a virtual environment and the install the required packages.

    Note: you will need the ``imas_core`` package to load IMAS data that is stored in
    the HDF5 or MDSplus backends. Unfortunately, this component is not publicly
    available, but can be installed provided you have access to its git repository.
    See
    https://sharepoint.iter.org/departments/POP/CM/IMDesign/Code%20Documentation/ACCESS-LAYER-doc/python/dev/building_installing.html#prerequisites
    (behind login wall) for the prerequisites.

    .. code-block:: bash

      python -m venv venv
      . venv/bin/activate
      pip install imas-paraview
      # Uninstall VTK, this would conflict with ParaView's built-in VTK module
      pip uninstall vtk
      # Optional: install imas-core package from the ITER git:
      git clone ssh://git@git.iter.org/imas/al-core.git -b main
      pip install ./al-core

4.  Set required environment variables and run ParaView:

    .. code-block:: bash

      export PYTHONPATH="$PWD/venv/lib/python3.10/site-packages/"
      export PV_PLUGIN_PATH="$PWD/venv/lib/python3.10/site-packages/imas_paraview/plugins"
      # This LD_PRELOAD is required when loading data with the imas_core HDF5 backend.
      # Note that it may break any built-in ParaView HDF5 support...
      export LD_PRELOAD="/lib/x86_64-linux-gnu/libhdf5_serial.so.103"
      ./ParaView-5.13.3-MPI-Linux-Python3.10-x86_64/bin/paraview

Development installation on SDCC
--------------------------------

This section provides instructions for installing and using a development version of the
IMAS-Paraview plugins on the ITER SDCC cluster.

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
