.. _`ci configuration`:

CI configuration
================

GGD-VTK uses `ITER Bamboo <https://ci.iter.org/>`_ for CI. This page provides an overview
of the CI Plan and deployment projects.

CI Plan
-------

The `GGD-VTK CI plan <https://ci.iter.org/browse/VIS-GGDVTK28>`_ consists of 3 types of jobs:

Linting 
    Run ``black``, ``flake8``, and ``isort`` on the vtkggdtools code base.
    See :ref:`code style and linting`.

    The CI script executed in this job is ``ci/linting.sh``.

Testing
    This runs all unit tests with pytest.

    The CI script executed in this job is ``ci/run_pytest.sh``, which expects the
    modules it needs to load as arguments. 

Benchmark
    This job runs the ASV benchmarks on the CI server. It
    is configured such that it can only run on a single CI agent
    (`io-ls-bamboowk11.iter.org`). There are two reasons for this:

    1.  Simplify the data I/O of the script - we can avoid file locks because there will
        only be a single Benchmark Job running globally. This is the main reason.
    2.  Benchmarks should be more reproduceable when always run on the same machine.
        Although the agents are virtualized, so the performance will always depend to
        some extent to the load on the CI cluster.

    The CI script executed in this job is: ``ci/run_benchmark.sh``.

Build docs
    This job builds the Sphinx documentation.

    The CI script executed in this job is: ``ci/build_docs.sh``, which expects the
    modules it needs to load as arguments.


Deployment projects
-------------------

TODO
.. There is currently one Bamboo deployment project for GGD-VTK:
..
.. `Deploy IDS-Validator-Doc <https://ci.iter.org/deploy/viewDeploymentProjectEnvironments.action?id=1908899843>`_
..     Deploy the documentation created in the `Build docs` job to `Sharepoint
..     <https://sharepoint.iter.org/departments/POP/CM/IMDesign/Code%20Documentation/IDS-Validator/index.html#>`_.
..
..     This deployment project runs for after each successful CI build of the IDS-validator main
..     branch.
