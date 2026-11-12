Choose a workflow
=================

.. admonition:: The short answer
   :class: important

   1. Run ``elastic-vasp check`` (:doc:`check`).
   2. If it says the recommended workflow **works in this orientation**, use it – e.g.
      ``hexagonal_pre`` / ``hexagonal_post``.
   3. Otherwise use the **automatic mode**, ``auto_pre`` / ``auto_post`` (:doc:`auto`) –
      it works for every crystal in any orientation, and also reports the textbook constants.

The rest of this page explains the workflows. Each one belongs to a range of space groups;
the number of deformations (VASP calculations per strain value) equals the number of
independent elastic constants.

.. list-table::
   :header-rows: 1
   :widths: 18 18 12 30 22

   * - Crystal system
     - Space groups
     - Constants
     - Functions
     - Laue class
   * - Triclinic
     - 1–2
     - 21
     - ``triclinic_pre`` / ``triclinic_post``
     - -1
   * - Monoclinic
     - 3–15
     - 13
     - ``monoclinic_pre`` / ``monoclinic_post``
     - 2/m
   * - Orthorhombic
     - 16–74
     - 9
     - ``orthogonal_pre`` / ``orthogonal_post``
     - mmm
   * - Tetragonal II
     - 75–88
     - 7
     - ``tetragonal_2_pre`` / ``tetragonal_2_post``
     - 4/m
   * - Tetragonal I
     - 89–142
     - 6
     - ``tetragonal_1_pre`` / ``tetragonal_1_post``
     - 4/mmm
   * - Trigonal II
     - 143–148
     - 7
     - ``trigonal_2_pre`` / ``trigonal_2_post``
     - -3
   * - Trigonal I
     - 149–167
     - 6
     - ``trigonal_1_pre`` / ``trigonal_1_post`` (see below)
     - -3m
   * - Hexagonal
     - 168–194
     - 5
     - ``hexagonal_pre`` / ``hexagonal_post``
     - 6/mmm
   * - Cubic
     - 195–230
     - 3
     - ``cubic_pre`` / ``cubic_post``
     - m-3m

If you are unsure, use :doc:`auto` – it works for every space group and orientation.

Things to know
--------------

**Tetragonal II (75–88) needs its own functions.** Point groups 4, -4 and 4/m have no
2-fold axes or mirrors along x and y, so the coupling C16 (with C26 = −C16) between
stretching and shearing in the basal plane is allowed. The orthorhombic deformations never
stretch and shear at the same time, so ``orthogonal_pre`` cannot measure C16 and would
return it as zero – the other constants come out right, but the matrix as a whole (and its
stability) does not. ``tetragonal_2_pre`` adds the deformation ``C11_C66_C16``.
C16 cannot be removed by choosing the orientation in advance: the rotation that makes it
vanish depends on the elastic constants themselves.

**Tetragonal I (89–142)** can also be done with ``orthogonal_pre`` (9 instead of 6
calculations); the extra constants then come out equal (C22 = C11, …) within the noise.

**Trigonal I "31m" groups.** In the standard setting of P312, P3\ :sub:`1`\ 12,
P3\ :sub:`2`\ 12, P31m, P31c, P-31m and P-31c (149, 151, 153, 157, 159, 162, 163) the
2-fold axes are perpendicular to a, so the constant that appears is C15 instead of C14.
Use ``trigonal_2_pre`` (it fits C14 and C15 together) for these, or rotate the cell by 90°
about c.

**Higher symmetry with a lower system is fine.** Any table that contains your crystal's
constants works, at the price of more calculations; e.g. ``triclinic_pre`` works for
everything.

**The functions check this for you.** Every named ``*_pre`` function uses spglib to check
whether its deformations can represent your crystal's stiffness matrix in its current
orientation. If not, it warns, names the constants that would be missed, and suggests
functions that work:

.. code-block:: text

   UserWarning: orthogonal_pre() cannot measure this crystal's stiffness matrix: POSCAR
   (space group P4/m (83)) has 7 independent elastic constants in this orientation, and C16
   can be non-zero, which the orthogonal deformations assume to be zero. Use
   tetragonal_2_pre() or monoclinic_pre() or auto_pre(), which works in any orientation.

.. admonition:: Next
   :class: tip

   :doc:`inputs` – the POSCAR and the strain values. (The full list of deformations for every workflow follows below, for reference.)

Deformations used
-----------------

Generated from the code (``elastic_vasp.systems``), so this is exactly what the functions
do. For each workflow: its independent constants, the form of the stiffness matrix
:math:`C`, and its deformations.

* **Deformation** – the matrix :math:`F = I + \varepsilon` by which the lattice vectors
  are multiplied in that folder (rows of the POSCAR lattice :math:`L` become
  :math:`L' = L\,F`), for each strain value :math:`\delta`. Its off-diagonal entries are
  half the engineering shear strains. A diagonal term like
  :math:`\delta^2/(1-\delta^2)` keeps the cell volume constant (see :doc:`../theory`).
* **Fitted curvature** – the coefficient :math:`a` of the fit
  :math:`E(\delta)/V_0 = a\,\delta^2`, which equals
  :math:`\tfrac12\,\mathbf e^{\mathsf T} C\,\mathbf e` for the Voigt strain
  :math:`\mathbf e` per unit :math:`\delta`. Post-processing solves these equations for the
  independent constants.

.. include:: ../_generated/system_tables.rst
