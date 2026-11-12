Orientation
===========

The strains are applied along the Cartesian x, y, z axes of the POSCAR lattice, and the
named functions assume the standard orientation of each crystal system:

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Functions
     - Required orientation
   * - ``cubic_*``
     - cube axes along x, y, z (conventional or the usual primitive fcc/bcc cells)
   * - ``orthogonal_*``
     - a ∥ x, b ∥ y, c ∥ z
   * - ``tetragonal_1_*``
     - c (4-fold axis) ∥ z and a ∥ x (or a at 45°, along [110] – x must lie on a 2-fold axis)
   * - ``tetragonal_2_*``
     - c (4-fold axis) ∥ z, any in-plane orientation (C16 is fitted)
   * - ``hexagonal_*``, ``trigonal_2_*``
     - c (6-fold / 3-fold axis) ∥ z, any in-plane orientation
   * - ``trigonal_1_*``
     - c ∥ z and x along a 2-fold axis (normal to a mirror plane)
   * - ``monoclinic_*``
     - unique axis (b, the 2-fold axis) ∥ y; a and c in the xz plane
   * - ``triclinic_*``
     - any

Monoclinic unique axis
----------------------

A monoclinic crystal has one special direction, the 2-fold rotation axis (or the normal to
its mirror plane): the *unique axis*. The lattice vector along it is perpendicular to the
other two. Most databases (and spglib) use the **b-unique** setting (α = γ = 90°, β ≠ 90°),
which is what ``monoclinic_*`` expects, with b along y:

.. code-block:: text

   a        0   0              <- a along x
   0        b   0              <- b (unique axis) along y
   c·cosβ   0   c·sinβ         <- c in the xz plane

With the 2-fold axis along y the non-zero off-diagonal constants are C15, C25, C35 and
C46; with a c-unique cell they would be C16, C26, C36 and C45 instead.

Automatic checks
----------------

The named ``*_pre`` functions check the orientation with **spglib**, using the lattice and
the atom positions (so a crystal whose lattice looks more symmetric than the crystal
itself, e.g. a monoclinic crystal in a cell with all angles 90°, is still judged
correctly):

* ``hexagonal_pre`` and ``trigonal_*_pre`` – the crystal's 3-fold / 6-fold axis must be
  along z. This also catches rhombohedral cells with the 3-fold axis along [111].
* ``monoclinic_pre`` – the crystal's 2-fold axis must be along y.
* every named function – the symmetry check described in :doc:`systems`.

If the axis is misaligned, the warning says where it actually is and a rotated copy is
written to ``POSCAR.reoriented``. The rotation is pure: fractional coordinates and atom
order are unchanged, so ``POTCAR`` and ``MAGMOM`` stay valid. Pre-processing then
**continues with the original orientation**, which gives wrong constants. To use the
rotated cell:

.. code-block:: bash

   cp POSCAR.reoriented POSCAR

and rerun the pre-processing with ``overwrite=True``.

Magnetic order is not taken into account: spglib sees only the atoms.

Example: cubic crystal along [1-10], [11-2], [111]
--------------------------------------------------

A cell often used for surfaces and dislocations has x ∥ [1-10], y ∥ [11-2] and
z ∥ [111]. For a cubic crystal the cube axes are then not along x, y, z, so
``cubic_pre`` does not apply (it warns, and its constants would be wrong). But the 3-fold
axis [111] is along z and the 2-fold axis [1-10] along x, so in this frame the crystal looks
**trigonal I**: ``trigonal_1_pre`` (6 deformations) measures it correctly, and
``auto_pre`` needs only 3. ``check_symmetry()`` lists exactly these options:

.. code-block:: text

   Space group:        Fm-3m (225) – cubic (Laue class m-3m)
   Elastic constants:  3 independent (C11, C12, C13 in this orientation)
   Recommended:        cubic_pre() / cubic_post()  – but NOT in this orientation ...
   Works as it is:     trigonal_1_pre(), trigonal_2_pre(), triclinic_pre()

``trigonal_1_post`` gives the matrix in this frame. ``auto_post`` gives it in this frame
**and** rotated back to the standard cubic setting, with the usual C11, C12, C44 – from the
same 3 deformations (see :doc:`auto`).

Not rotating at all
-------------------

:doc:`auto` determines the constants and deformations for the orientation you have, so
nothing needs to be rotated. It reports the stiffness matrix in your Cartesian frame (for
a hexagonal crystal with c along x the usual C33 appears as C11) and also in the crystal's
standard setting, with the textbook constants.

.. admonition:: Next
   :class: tip

   :doc:`auto` – skip the rotation entirely and get both the matrix in your frame and the textbook constants.
