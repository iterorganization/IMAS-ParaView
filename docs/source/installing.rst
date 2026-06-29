.. _`installing`:

Installing IMAS-ParaView
========================

If the IMAS-ParaView module is available on your HPC environment, you can ignore the further 
steps below and simply load the module as follows:

.. code-block:: bash

  module load IMAS-ParaView
  # Run paraview (add vglrun to enable hardware acceleration)
  vglrun paraview
  # Open the "Sources" tab in the top left, and ensure you see "IMAS Tools" in 
  # the drop down menu.
  # The command-line interface can be started using:
  imas-paraview --help


Installation from PyPI
----------------------

IMAS-ParaView is published on `PyPI <https://pypi.org/project/imas-paraview/>`_ and can be
installed into a virtual environment that uses ParaView's own Python interpreter
(``pvpython``).

1.  Download ParaView 5.13.3 from the ParaView website:
    https://www.paraview.org/download/?version=v5.13&filter=Linux

2.  Unzip the ParaView package and enter the extracted directory:

    .. code-block:: bash

      tar -xzvf ParaView-5.13.3-MPI-Linux-Python3.10-x86_64.tar.gz
      cd ParaView-5.13.3-MPI-Linux-Python3.10-x86_64

3.  Create a virtual environment using ParaView's bundled ``pvpython`` interpreter and
    install IMAS-ParaView from PyPI. This example uses
    `uv <https://docs.astral.sh/uv/>`_ (``pip install uv``), but ``python -m venv`` and
    ``pip`` work too.

    .. code-block:: bash

      uv venv -p ./bin/pvpython
      uv pip install imas-paraview

4.  Activate the virtual environment and set the required environment variables.

    .. code-block:: bash

      . .venv/bin/activate
      export PYTHONPATH="$(python -c 'import sysconfig; print(sysconfig.get_path("purelib"))')"
      export PV_PLUGIN_PATH="$(python -c 'import imas_paraview; print(imas_paraview.__path__[0])')/plugins"
      # This LD_PRELOAD is required when loading data with the imas_core HDF5 backend.
      # Note that it may break any built-in ParaView HDF5 support...
      export LD_PRELOAD="$(ls "$PYTHONPATH"/imas_core.libs/libhdf5-*.so.*)"

5.  Run ParaView:

    .. code-block:: bash

      ./bin/paraview
      # Open the "Sources" tab in the top left, and ensure you see "IMAS Tools" in
      # the drop down menu.

    The command-line interface is available in the activated environment as well:

    .. code-block:: bash

      imas-paraview --help

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

  # Load compatible IMAS-Python and ParaView modules, like:
  module load IMAS-Python/2.3.0-foss-2023b ParaView/5.12.0-foss-2023b
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

  # Load compatible IMAS-Python and ParaView modules, like:
  module load IMAS-Python/2.3.0-foss-2023b ParaView/5.12.0-foss-2023b
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
  module load IMAS-Python/2.3.0-foss-2023b ParaView/5.12.0-foss-2023b
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

  # The tests require downloading external IDSs from Zenodo
  wget -P data/ -i zenodo_datasets.txt

  # Run the tests
  python -m pytest

* Alternatively, if you want to skip running the integration tests and the 
  tests that contain external data:

.. code-block:: bash

  python -m pytest -m "not integration and not external_data"

* To build the IMAS-ParaView documentation, ensure the optional docs dependencies are pip 
  installed (or simply use all, to install all optional dependencies).

.. code-block:: bash

  make -C docs html
  # You can now open ./docs/_build/html/index.html
