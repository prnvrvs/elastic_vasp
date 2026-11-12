Maintainer notes
================

Development install
-------------------

.. code-block:: bash

   pip install -e ".[test,doc]"

Tests
-----

.. code-block:: bash

   pytest

The suite runs without VASP: it generates the strained POSCARs and writes synthetic
OUTCARs whose energies follow a known stiffness matrix. The main groups:

* ``test_elastic.py`` – every crystal system end to end; energy reading and run checks.
* ``test_reference.py`` – compares generated structures and constants with
  ``tests/data/reference.json`` (recorded from earlier code). Regenerate only for an
  intended change: ``python tests/test_reference.py``.
* ``test_systems.py`` – consistency of the deformation tables.
* ``test_inputs.py`` – POSCAR forms, ``strain.dat``, volumes, OUTCAR checks.
* ``test_orientation.py`` – orientation checks and ``POSCAR.reoriented``.
* ``test_auto.py`` – the automatic mode, including rotated crystals and noisy energies.
* ``test_cli.py`` – ``check_symmetry``, the ``strains`` argument and the ``elastic-vasp``
  command line.
* ``test_spacegroups.py`` – all 230 space groups: constant counts, rotations, error
  amplification, and the symmetry check for every named table.

Documentation
-------------

.. code-block:: bash

   make -C docs html          # output in docs/_build/html

The deformation tables in :doc:`guide/systems` are generated from
``elastic_vasp.systems`` at build time (``docs/_generated``). The ``Docs`` GitHub workflow
builds the documentation on every push and deploys it to GitHub Pages from ``main``
(the repository's Pages source must be set to *GitHub Actions*).

Releasing
---------

1. Update ``__version__`` in ``elastic_vasp/__init__.py`` (``setup.py`` reads it) and
   :doc:`release_notes`.
2. Commit, tag ``vX.Y.Z`` and push both.
3. Create a GitHub release from the tag. The ``Python Package`` workflow runs the tests
   and publishes to PyPI only for releases, and only if the tests pass.
