Theory
======

Energy–strain relation
----------------------

For small homogeneous strains of a stress-free crystal the energy per reference volume is

.. math::

   \frac{E(\boldsymbol\varepsilon) - E_0}{V_0} = \frac12 \sum_{i,j=1}^{6} C_{ij}\,\varepsilon_i\,\varepsilon_j
   + \mathcal{O}(\varepsilon^3),

with the strain in Voigt notation, :math:`\varepsilon_1 = \varepsilon_{xx}`, …,
:math:`\varepsilon_3 = \varepsilon_{zz}` and engineering shear strains
:math:`\varepsilon_4 = 2\varepsilon_{yz}`, :math:`\varepsilon_5 = 2\varepsilon_{xz}`,
:math:`\varepsilon_6 = 2\varepsilon_{xy}`.

Each deformation uses a fixed strain pattern :math:`\mathbf e` scaled by the strain values
:math:`\delta` (from ``strain.dat`` or ``--strains``), :math:`\boldsymbol\varepsilon = \delta\,\mathbf e`, so

.. math::

   \frac{E(\delta) - E_0}{V_0} = a\, \delta^2 + \mathcal{O}(\delta^3), \qquad
   a = \tfrac12\, \mathbf e^{\mathsf T} C\, \mathbf e .

Fitting :math:`E(\delta)/V_0` with a quadratic gives :math:`a` (converted from eV/Å³ to GPa
with 1 eV/Å³ = 160.2177 GPa). For example, for a cubic crystal the pattern
:math:`\mathbf e = (1, -1, 0, 0, 0, 0)` gives :math:`a = C_{11} - C_{12}`.

Applying the strain
-------------------

With the lattice vectors as rows of :math:`L`, the deformed lattice is :math:`L' = L\,F` with

.. math::

   F = I + \boldsymbol\varepsilon =
   \begin{pmatrix}
   1+\varepsilon_1 & \varepsilon_6/2 & \varepsilon_5/2\\
   \varepsilon_6/2 & 1+\varepsilon_2 & \varepsilon_4/2\\
   \varepsilon_5/2 & \varepsilon_4/2 & 1+\varepsilon_3
   \end{pmatrix},

the symmetric strain tensor built from the Voigt strain (off-diagonal entries are half the
engineering shear strains). The atoms keep their fractional coordinates, so they move with
the cell (with ``ISIF = 2`` VASP then relaxes them inside the fixed cell). The test suite
checks for every deformation of every crystal system that the written POSCAR lattice is
exactly :math:`L\,F`.

Volume-conserving corrections
-----------------------------

Some deformations contain a second-order diagonal term that keeps the cell volume
:math:`\det F` exactly 1. For opposite normal strains,

.. math::

   F = \begin{pmatrix} 1+\delta & 0 & 0 \\ 0 & 1-\delta & 0 \\ 0 & 0 & 1+\varepsilon_3 \end{pmatrix},
   \qquad
   \det F = (1+\delta)(1-\delta)(1+\varepsilon_3) = 1
   \;\Rightarrow\; \varepsilon_3 = \frac{\delta^2}{1-\delta^2},

and for a pure shear :math:`\varepsilon_4 = \delta`,

.. math::

   F = \begin{pmatrix} 1+\varepsilon_1 & 0 & 0 \\ 0 & 1 & \tfrac{\delta}{2} \\ 0 & \tfrac{\delta}{2} & 1 \end{pmatrix},
   \qquad
   \det F = (1+\varepsilon_1)\Bigl(1 - \frac{\delta^2}{4}\Bigr) = 1
   \;\Rightarrow\; \varepsilon_1 = \frac{\delta^2}{4-\delta^2}.

The correction is of second order in :math:`\delta`, so it changes the energy only at order
:math:`\delta^3` and above and leaves the curvature :math:`a` unchanged. Its purpose is to
remove the effect of a small residual pressure :math:`P` in the reference cell: a volume
change :math:`\Delta V \propto \delta^2` would add :math:`P\,\Delta V` directly to the
:math:`\delta^2` term. Deformations meant to change the volume (e.g.
:math:`F = (1+\delta)\,I` for the bulk modulus) have no such term.

Solving for the constants
-------------------------

The independent constants :math:`\mathbf p` of a crystal system (or of the symmetry-adapted
basis in :doc:`guide/auto`) enter :math:`C` linearly, :math:`C = \sum_k p_k B_k`, so every
deformation gives one linear equation

.. math::

   a_j = \sum_k \left(\tfrac12\, \mathbf e_j^{\mathsf T} B_k\, \mathbf e_j\right) p_k ,

and the set of deformations is chosen so that this square system has full rank. The
largest factor by which an error in the :math:`a_j` can grow in :math:`C` (the *error
amplification*) is 1.4–5.1 for the named tables and at most about 6 for the automatic mode
(checked for all 230 space groups in random orientations).

Symmetry
--------

A rotation :math:`R` of the crystal transforms the Voigt stiffness matrix as
:math:`C' = M C M^{\mathsf T}`, with the 6×6 matrix

.. math::

   M_{IJ} = R_{ik} R_{jl} \;(J \le 3), \qquad
   M_{IJ} = R_{ik} R_{jl} + R_{il} R_{jk} \;(J > 3),

for the index pairs :math:`I \leftrightarrow (ij)`, :math:`J \leftrightarrow (kl)`. The
allowed stiffness matrices of a crystal are those with :math:`M C M^{\mathsf T} = C` for
every operation of its point group. This is how the symmetry checks and the automatic mode
work; it gives 21, 13, 9, 7, 6, 7, 6, 5 and 3 independent constants for the Laue classes
-1, 2/m, mmm, 4/m, 4/mmm, -3, -3m, 6/mmm and m-3m.

Mechanical stability
--------------------

A crystal is mechanically stable (Born criteria) when the strain energy is positive for
every small strain, i.e. :math:`C` is positive definite: all its eigenvalues are positive.
For cubic crystals this is :math:`C_{11} - C_{12} > 0`, :math:`C_{11} + 2C_{12} > 0`,
:math:`C_{44} > 0`; for hexagonal crystals :math:`C_{11} > |C_{12}|`,
:math:`2C_{13}^2 < C_{33}(C_{11}+C_{12})`, :math:`C_{44} > 0`, :math:`C_{66} > 0`.

Polycrystalline averages
------------------------

With the compliance :math:`S = C^{-1}`:

.. math::

   9K_V &= (C_{11}+C_{22}+C_{33}) + 2(C_{12}+C_{13}+C_{23})\\
   15G_V &= (C_{11}+C_{22}+C_{33}) - (C_{12}+C_{13}+C_{23}) + 3(C_{44}+C_{55}+C_{66})\\
   1/K_R &= (S_{11}+S_{22}+S_{33}) + 2(S_{12}+S_{13}+S_{23})\\
   15/G_R &= 4(S_{11}+S_{22}+S_{33}) - 4(S_{12}+S_{13}+S_{23}) + 3(S_{44}+S_{55}+S_{66})

The Hill values are the means :math:`K_H = (K_V+K_R)/2`, :math:`G_H = (G_V+G_R)/2`, and

.. math::

   E_H = \frac{9K_H G_H}{3K_H + G_H}, \qquad
   \nu_H = \frac{3K_H - 2G_H}{2(3K_H + G_H)}, \qquad
   \text{Pugh ratio} = K_H / G_H .

Limitations
-----------

* The reference structure must be relaxed; residual stress biases deformations that change
  the volume.
* Anharmonic terms bias the quadratic fit at larger strains (a few tenths of a percent at
  ±3 % for typical materials).
* The symmetry analysis ignores magnetic order.
