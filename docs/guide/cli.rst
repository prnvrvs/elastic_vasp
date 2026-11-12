Command line
============

Installing the package provides the ``elastic-vasp`` command (also available as
``python -m elastic_vasp``). Every step of the workflow can be run from the shell:

.. code-block:: bash

   elastic-vasp check                                   # what is this structure?
   elastic-vasp pre cubic                               # strained POSCARs
   # ... run VASP in every strain_* folder ...
   elastic-vasp post cubic --output C.txt               # constants; save the 6x6 matrix

Commands
--------

``elastic-vasp check [POSCAR] [--symprec S]``
   Space group, number of constants, and which workflow works with the structure as
   oriented. Read-only. See :doc:`check`.

``elastic-vasp pre SYSTEM [options]``
   Pre-processing. ``SYSTEM`` is one of ``auto``, ``cubic``, ``hexagonal``,
   ``tetragonal_1``, ``tetragonal_2``, ``trigonal_1``, ``trigonal_2``, ``orthogonal``,
   ``monoclinic``, ``triclinic``.

   ``--strains S [S ...]``
      the strain values, e.g. ``--strains -0.02 -0.01 0 0.01 0.02``
   ``--max-strain S`` and ``--points N``
      ``N`` evenly spaced strains from ``-S`` to ``+S`` (default ``N`` = 7)
   ``--overwrite``
      replace existing deformation folders (deletes their contents)
   ``--symprec S``
      spglib tolerance for ``auto``

   Without ``--strains`` / ``--max-strain`` the strains are read from ``strain.dat``.

``elastic-vasp post SYSTEM [-o FILE] [--frame input|standard]``
   Post-processing with the same ``SYSTEM`` as ``pre``. ``-o FILE`` / ``--output FILE``
   also saves the 6×6 stiffness matrix (GPa) as plain text – six rows, no header – readable
   with ``numpy.loadtxt`` and directly by elastic-vis (:doc:`visualise`). For ``auto``,
   ``--frame standard`` saves the matrix in the crystal's
   standard setting instead of your POSCAR's orientation (both are always printed).

``elastic-vasp systems``
   Lists every crystal-system workflow with its constants and deformations.

Global options
--------------

``-C DIR`` / ``--directory DIR``
   run in ``DIR`` instead of the current folder, e.g. ``elastic-vasp -C calc/Zr post hexagonal``
``--version``
   print the installed version
``-h`` / ``--help``
   show help; also works per command, e.g. ``elastic-vasp pre --help``

Messages and exit codes
-----------------------

Warnings (for example from the symmetry check) are printed as ``warning: …`` and the
command continues. Problems that stop a step – a missing POSCAR, existing folders without
``--overwrite``, bad strain values, unfinished VASP runs – are printed as ``error: …`` and
the command exits with status 1, so it can be used in job scripts:

.. code-block:: bash

   elastic-vasp post hexagonal --output C.txt || echo "some runs are not finished yet"

Example: a whole calculation
----------------------------

.. code-block:: bash

   cd my_structure                # POSCAR, INCAR, KPOINTS, POTCAR
   elastic-vasp check
   elastic-vasp pre auto --max-strain 0.02 --points 5
   for d in */strain_*/; do (cd "$d" && ln -sf ../../INCAR ../../KPOINTS ../../POTCAR . && sbatch ../../job.sh); done
   # ... when all jobs are done:
   elastic-vasp post auto --output C.txt
