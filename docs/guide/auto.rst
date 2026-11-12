Automatic mode
==============

``auto_pre`` and ``auto_post`` take the crystal as it is in your POSCAR – in any
orientation, for any space group – and work out what to calculate:

.. code-block:: python

   from elastic_vasp import elastic_gen as eg, constants

   eg.auto_pre()              # pre-processing
   C = constants.auto_post()  # after the VASP runs

What ``auto_pre`` does
----------------------

1. **Finds the crystal symmetry** with spglib from the lattice and the atoms, in your
   Cartesian frame. Nothing is rotated.
2. **Derives the independent constants.** The stiffness matrices that are unchanged by
   every symmetry operation form a linear subspace; its reduced row-echelon basis over the
   entries C11, C12, …, C66 gives the independent constants and how every other entry
   depends on them. In standard orientations these are exactly the textbook sets
   (C11, C12, C44 for cubic, …); in other orientations other entries appear (e.g. C16 for
   a monoclinic crystal with its unique axis along z).
3. **Chooses a minimal set of deformations** – one per independent constant – from simple
   patterns (one strain component, or two with equal or opposite sign). They are chosen to
   keep the error amplification small: an error of 1 GPa in a fitted curvature changes
   the stiffness matrix by at most about 6 GPa for all 230 space groups in any
   orientation, similar to the hand-made tables (1.4–5.1).
4. **Prints** the space group, the constants, the matrix pattern (for simple
   orientations), the deformations and the error amplification, creates the folders and
   saves everything to ``elastic_vasp_auto.json``.

Example output for a hexagonal crystal in the standard orientation:

.. code-block:: text

   Crystal symmetry (spglib, symprec=0.001 A): P6/mmm (191)
   5 independent elastic constants in this orientation: C11, C12, C13, C33, C44
   Stiffness matrix pattern:
       C11   C12   C13   0     0     0
       C12   C11   C13   0     0     0
       C13   C13   C33   0     0     0
       0     0     0     C44   0     0
       0     0     0     0     C44   0
       0     0     0     0     0     0.5*C11 - 0.5*C12
   Deformations (5):
       e3       (0, 0, 1, 0, 0, 0) * delta                              a = 0.5*C33
       e1+e2    (1, 1, 0, 0, 0, 0) * delta                              a = C11 + C12
       e1-e2    (1, -1, 0, 0, 0, 0) * delta, e3 = delta^2/(1-delta^2)   a = C11 - C12
       e4       (0, 0, 0, 1, 0, 0) * delta, e1 = delta^2/(4-delta^2)    a = 0.5*C44
       e1+e3    (1, 0, 1, 0, 0, 0) * delta                              a = 0.5*C11 + C13 + 0.5*C33
   Error amplification: at most 3.3 GPa in C per GPa error in a fitted coefficient.

What ``auto_post`` does
-----------------------

Reads ``elastic_vasp_auto.json`` – so exactly the deformations that were generated are
used, even if POSCAR changed since – and then does the same checks, fits and reports as the
named ``*_post`` functions. It prints the stiffness matrix twice:

* **in the orientation of your POSCAR**, and
* **in the crystal's standard setting** – cube axes along x, y, z; hexagonal and trigonal
  c along z with a along x; monoclinic unique axis b along y; … – together with the
  textbook constants of its crystal system, exactly as the named ``*_post`` function would
  report them.

Both describe the same crystal: the rotation between the two frames comes from spglib and
is exact, and the stability check and the polycrystalline averages are the same in either
frame. For example, fcc Al in a cell with x ∥ [1-10], y ∥ [11-2], z ∥ [111]:

.. code-block:: text

   ==== stiffness matrix ====                 (your frame)
   [[265.02 114.98 109.97  -7.09   0.     0.  ]
    ...
   ==== in the crystal's standard setting ====
   space group 225 - cubic (Laue class m-3m); rotation from your cell to the standard frame:
   ...
   cubic constants (as from cubic_post()):
   C11= 249.9880
   C12= 119.9899
   C44= 80.0300

The return value is the matrix in your frame (``frame="input"``, default) or in the
standard setting:

.. code-block:: python

   C_mine = constants.auto_post()
   C_std = constants.auto_post(frame="standard")

So with a single set of 3 (cubic) … 21 (triclinic) deformations on your own oriented cell
you get both the matrix for that cell and the standard constants to compare with tables
and papers. For low-symmetry crystals the standard constants depend on the standard choice
of axes (e.g. which lattice vector is a); for cubic and hexagonal crystals they are unique.

When to use it
--------------

* The structure is not in the standard orientation, and you don't want to rotate it.
* You are unsure which crystal system or setting applies (e.g. trigonal "31m" groups,
  tetragonal II).
* Low-symmetry or unusual settings (c-unique monoclinic, rhombohedral cells along [111]).

Notes
-----

* Folder names describe the strain pattern: ``e1``, ``e4``, ``e1+e2``, ``e5-e6``, …
* ``symprec`` (default 1e-3 Å) is the spglib tolerance; pass ``auto_pre(symprec=...)`` for
  slightly distorted structures. Too large a tolerance can find symmetry the crystal does
  not have, which would make some constants wrongly zero.
* Magnetic order is not considered by the symmetry search.
* For tables and papers, use the standard-setting constants that ``auto_post`` prints
  (or ``frame="standard"``); the matrix in your frame has the same physics but other
  labels (e.g. the hexagonal C33 appears as C11 when c lies along x).

.. admonition:: Next
   :class: tip

   Back to the main path: :doc:`inputs` and :doc:`workflow` work the same way with ``auto``.
