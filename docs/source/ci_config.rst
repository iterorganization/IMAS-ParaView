.. _`ci configuration`:

CI configuration
================

GGD-VTK uses `ITER Bamboo <https://ci.iter.org/>`_ for CI. This page provides an overview
of the CI Plan and deployment projects. Some of the jobs in the CI Plan can also be run manually,
examples provided below.

CI Plan
-------

The `GGD-VTK CI plan <https://ci.iter.org/browse/VIS-GGDVTK>`_ consists of 3 types of jobs:

Linting 
    Run ``black``, ``flake8``, and ``isort`` on the vtkggdtools code base.
    See :ref:`code style and linting`.

    The CI script executed in this job is ``ci/linting.sh``.

    **Example:**
    ``$ ci/linting.sh Python/3.11.5-GCCcore-13.2.0``

Testing
    This runs all unit tests with pytest.

    The CI script executed in this job is ``ci/run_pytest.sh``, which expects the
    modules it needs to load as arguments.

    It needs a virtual framebuffer implementation. This is currently
    not working reliably on the CI nodes, and so this job as been
    disabled for now. To run the complete suite, including the virtual
    framebuffer tests, make sure you have all the environment modules
    loaded and run ``pytest`` in the root folder of the repository.

    **Example:**
    
    ``$ ci/run_pytest.sh IMASPy IMAS-AL-Python ParaView Xvfb # identical run as in the CI``

    ``$ module load IMASPy IMAS-AL-Python ParaView Xvfb``
    ``$ pytest # enforcing the vfb tests``

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

    **Example:**
    ``$ ci/run_benchmark.sh IMASPy``

Build docs
    This job builds the Sphinx documentation.

    The CI script executed in this job is: ``ci/build_docs.sh``, which expects the
    modules it needs to load as arguments.

    **Example:**
    ``$ ci/build_docs.sh IMASPy``


Deployment projects
-------------------

There is currently one Bamboo deployment project for GGD-VTK:

`Deploy GGD-VTK Docs <https://ci.iter.org/deploy/viewDeploymentProjectEnvironments.action?id=1942093825>`_
    Deploy the documentation created in the `Build docs` job to `Sharepoint
    <https://sharepoint.iter.org/departments/POP/CM/IMDesign/Code%20Documentation/GGD-VTK/index.html#>`_.

    This deployment project runs for after each successful CI build of the GGD-VTK main
    branch.
