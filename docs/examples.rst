Examples
========

The ``example/`` folder of the repository contains ready-to-run inputs (without
``POTCAR``, which you must copy from your VASP installation):

.. list-table::
   :header-rows: 1
   :widths: 20 40 40

   * - Folder
     - Structure
     - Functions
   * - ``Cubic``
     - fcc Al, conventional cell (Fm-3m)
     - ``cubic_pre`` / ``cubic_post``
   * - ``Hexagonal``
     - hcp Zr (P6\ :sub:`3`/mmc)
     - ``hexagonal_pre`` / ``hexagonal_post``
   * - ``Orthogonal``
     - fcc Al – a cubic crystal run through the orthorhombic workflow (9 constants;
       C11 = C22 = C33 etc. should come out equal)
     - ``orthogonal_pre`` / ``orthogonal_post``
   * - ``Trigonal``
     - α-Al\ :sub:`2`\ O\ :sub:`3` (corundum, R-3c), 10-atom primitive cell, c ∥ z, a ∥ x
     - ``trigonal_1_pre`` or ``trigonal_2_pre``
   * - ``Triclinic``
     - fcc Al – a cubic crystal run through the general triclinic workflow (21 constants)
     - ``triclinic_pre`` / ``triclinic_post``

The ``Orthogonal`` and ``Triclinic`` folders use a cubic crystal on purpose: a workflow for
a lower symmetry works for any higher-symmetry crystal, and comparing the results with the
cubic run is a good consistency check.

Each folder has:

* ``POSCAR``, ``INCAR``, ``KPOINTS``, ``strain.dat``;
* ``test.py`` – uncomment the pre-processing call, run it, then (after VASP) comment it and
  uncomment the post-processing call;
* ``copy.sh`` – links ``INCAR``, ``POTCAR``, ``KPOINTS`` and ``run.pbs`` into every strain
  folder and submits the jobs (edit the submission command for your cluster);
* ``run.pbs`` – an example job script (not in ``Hexagonal``; copy one from another folder).

Instead of editing ``test.py`` you can use the command line, e.g. ``elastic-vasp check``,
``elastic-vasp pre cubic`` and ``elastic-vasp post cubic`` in the example folder (see
:doc:`guide/cli`).

Running an example
------------------

.. code-block:: bash

   cd example/Cubic
   cp /path/to/your/POTCAR .
   python3 test.py            # with eg.cubic_pre() uncommented
   bash -i copy.sh            # submit the VASP runs
   # ... after all runs finished, switch test.py to constants.cubic_post()
   python3 test.py

The ``INCAR`` and ``KPOINTS`` files are starting points; converge ``ENCUT`` and the
k-mesh for your own material.
