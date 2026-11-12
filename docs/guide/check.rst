Check your structure
====================

Start every calculation here. One command tells you the space group, how many elastic
constants the crystal has, and which workflow to use with the POSCAR as it is:

.. code-block:: python

   from elastic_vasp import elastic_gen as eg

   report = eg.check_symmetry()          # reads ./POSCAR; or check_symmetry("path/to/POSCAR")

or on the command line:

.. code-block:: bash

   elastic-vasp check

Example output for the trigonal example (corundum):

.. code-block:: text

   POSCAR: Al4O6 (10 atoms)
   Space group:        R-3c (167) – trigonal I (Laue class -3m), spglib symprec = 0.001 A
   Elastic constants:  6 independent (C11, C12, C13, C14, C33, C44 in this orientation)
   Recommended:        trigonal_1_pre() / trigonal_1_post()  – works in this orientation
   Also works:         trigonal_2, triclinic (more calculations)
   Any orientation:    auto_pre() / auto_post() – 6 deformations

and for a hexagonal crystal whose c axis lies along x:

.. code-block:: text

   Space group:        P6_3/mmc (194) – hexagonal (Laue class 6/mmm), spglib symprec = 0.001 A
   Elastic constants:  5 independent (C11, C12, C22, C23, C55 in this orientation)
   Recommended:        hexagonal_pre() / hexagonal_post()  – but NOT in this orientation: ...
   Works as it is:     orthogonal_pre(), monoclinic_pre(), triclinic_pre()
   To use hexagonal_pre(), rotate the structure first (its warning writes a rotated
   POSCAR.reoriented), or use the automatic mode.

What it reports
---------------

* **Space group** (with spglib, from the lattice and the atoms), the crystal system and the
  Laue class, which decides the form of the stiffness matrix.
* **The independent elastic constants** for the structure **as oriented in the file**.
* **The recommended named workflow** for that space group (see :doc:`systems`; for the
  trigonal "31m" groups this is ``trigonal_2``), and whether it can measure the stiffness
  matrix with the crystal in its current orientation.
* **Every named workflow that works as it is**, from the fewest calculations up.
* The automatic mode, which always works.

``check_symmetry`` only reads the file; it never writes anything. It returns a
``SymmetryReport`` with the same information as attributes (``space_group``, ``number``,
``crystal_system``, ``laue_class``, ``constants``, ``recommended``, ``recommended_works``,
``working``); pass ``verbose=False`` to skip the printout.

Tolerance
---------

``symprec`` (default 1e-3 Å) is spglib's distance tolerance. Relaxed structures are
symmetric only to about 10\ :sup:`-7` Å; such numerical noise is ignored. Use a larger
value for structures that are slightly distorted on purpose, but not so large that spglib
finds symmetry the crystal does not have.

.. admonition:: Next
   :class: tip

   :doc:`systems` – what the recommendation means and when to pick something else.
