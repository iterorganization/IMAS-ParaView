.. _`ci configuration`:

CI configuration
================

IMAS-Paraview uses `GitHub Actions <https://github.com/features/actions>`_ for CI. This page provides an overview
of the CI Plan and deployment processes. Some of the jobs in the CI Plan can also be run manually,
examples provided below.

CI Plan
-------

The GitHub Actions workflow file is defined in
`.github/workflows/ci.yml <https://github.com/iterorganization/IMAS-ParaView/blob/develop/.github/workflows/ci.yml>`_.
It includes the following jobs:

Linting
    Runs ``ruff`` for both lint checks and code formatting.

    This is automatically triggered on push. Example manual run:

    .. code-block:: console

        $ ruff check --output-format=github imas_paraview
        $ ruff format --diff imas_paraview

Testing
    Executes unit tests with ``pytest`` using the `kitware/paraview-for-ci` container.

    Tests are run with coverage reports and JUnit XML output. Example:

    .. code-block:: console

        $ /opt/paraview/install/bin/pvpython -m pytest \
            --cov=imas_paraview \
            --cov-report=term-missing \
            --cov-report=html \
            --junit-xml=junit-3.11.xml

Benchmarking
    Runs ASV benchmarks on the CI server. Benchmarks compare current changes against
    `develop`, `master`, and HEAD.

    Includes setup of a reproducible machine configuration and publishes HTML results.

    Example:

    .. code-block:: console

        $ asv machine --yes
        $ asv run -v --show-stderr --skip-existing-successful HEAD^!
        $ asv compare develop HEAD
        $ asv publish

Documentation
    Builds the Sphinx documentation and uploads the HTML output as an artifact.

    Example:

    .. code-block:: console

        $ make -C docs html
