.. _`code style and linting`:

Code style and linting
======================

Code style
----------

IMAS-Paraview uses `Ruff <https://docs.astral.sh/ruff/>`_ for both code formatting and linting. All Python files must be formatted and linted using the ``ruff`` command-line tool. This is enforced in :ref:`CI <ci configuration>`.

Why Ruff?
'''''''''

We use Ruff to ensure that code style is uniform across all Python files, regardless of the developer who wrote the code.

This improves the efficiency of developers working on the project:

- Uniform code style makes it easier to read, review, and understand others' code.
- Autoformatting code reduces time spent on style decisions, allowing developers to focus on logic and functionality.
- Static analysis detects common issues before runtime, preventing certain classes of bugs.

Using Ruff
''''''''''

The easiest way to work with Ruff is via editor integration. See https://docs.astral.sh/ruff/editors/#editor-integrations for details.

You can also install Ruff and run it manually before committing:

.. code-block:: console

    $ ruff format imas_paraview       # Format code
    $ ruff check imas_paraview        # Lint code

Docstring style
---------------

While not enforced, we recommend using `Napoleon-style docstrings <https://sphinxcontrib-napoleon.readthedocs.io/en/latest/>`_.
