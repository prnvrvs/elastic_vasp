Legacy version
==============

The original code is still available, so results produced with it can be reproduced
exactly – for example those of the papers listed in :doc:`credit`.

Where to get it
---------------

* **Git branch** ``legacy`` of the repository (version 0.0.4):

  .. code-block:: bash

     git clone -b legacy https://github.com/prnvrvs/elastic_vasp.git

  or install it directly:

  .. code-block:: bash

     pip install "git+https://github.com/prnvrvs/elastic_vasp.git@legacy"

* **PyPI**: the earlier releases remain installable, e.g.

  .. code-block:: bash

     pip install "elastic-vasp==0.0.4"

Use a separate virtual environment if you need both versions side by side.

The function names are the same as today (``cubic_pre`` / ``cubic_post``, …), and every
crystal system uses the same folder names and strain patterns, with one exception: the
monoclinic folder ``C44_C66_C46_II`` (see below).

**Re-analysing legacy calculations.** Because the deformations are the same, finished
calculations made with the legacy version can be analysed with the current ``*_post``
functions to get corrected constants – for trigonal II and triclinic crystals nothing has to
be recalculated; for monoclinic crystals only the ``C44_C66_C46_II`` folder must be
regenerated and recalculated. The current version needs a finished ``OUTCAR`` in every
strain folder (see :doc:`guide/postprocessing`).

Known issues
------------

These were found and fixed in later versions. **Use the current version for new work**, and
re-check earlier results for the affected crystal systems:

.. list-table::
   :header-rows: 1
   :widths: 22 50 28

   * - Crystal system
     - Problem in the legacy version
     - Fixed in
   * - Trigonal II (``trigonal_2``)
     - C15 has the wrong sign and is left out of the stiffness matrix (C15, C25, C46)
     - 0.0.6
   * - Monoclinic
     - C15 and C46 have the wrong sign; the ``C44_C66_C46_II`` deformation strains
       e5/e6 instead of e4/e6, so C66 and C46 are wrong
     - 0.0.6 (folder ``C44_C66_C46_II`` must be recalculated)
   * - Triclinic
     - the fifth row of the printed/returned matrix has C25 and C35 where C35 and C45
       belong, so the matrix is not symmetric
     - 0.0.5
   * - Tetragonal II (75–88)
     - no dedicated workflow; the orthorhombic one cannot measure C16
     - current (``tetragonal_2_*``)

Cubic, hexagonal, orthorhombic and trigonal I results are not affected by these issues.

Other differences
-----------------

The legacy version has none of the later checks and conveniences:

* energies are read with ``grep``/``awk`` from the last line of ``OSZICAR`` (the E0 value,
  the same quantity as today's ``energy(sigma->0)``), so it needs a Unix shell, and
  unfinished or unconverged runs are not detected;
* the reference volume is read from ``strain_0.00/OUTCAR``, so that folder must exist
  under exactly that name;
* existing deformation folders are deleted without asking;
* POSCAR must use Direct coordinates and a single positive scale factor;
* no symmetry or orientation checks, no tetragonal or automatic workflows.

See :doc:`release_notes` for everything that changed since.
