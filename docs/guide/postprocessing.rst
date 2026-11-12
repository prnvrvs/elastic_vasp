Analyse the results
===================

Run the ``*_post`` function that matches the ``*_pre`` function you used, in the same
folder:

.. code-block:: python

   from elastic_vasp import constants

   C = constants.hexagonal_post()     # or constants.auto_post()

Checks on the VASP runs
-----------------------

Before anything is fitted, every ``strain_*`` folder of every deformation is checked. All
problems are collected and reported together, and post-processing stops
(``MissingEnergyError``):

* a deformation folder without ``strain_*`` folders, or with fewer than 3;
* a missing ``OUTCAR``;
* no ``energy(sigma->0)`` line (the run crashed before finishing an SCF cycle);
* no ``General timing and accounting`` footer – the run did not terminate normally
  (killed, out of walltime, still running);
* an ionic relaxation (``NSW > 0`` and ``IBRION`` = 1, 2 or 3) without
  ``reached required accuracy`` – it ran out of steps. Increase ``NSW`` and restart from
  ``CONTCAR``.

.. code-block:: text

   MissingEnergyError: Cannot post-process cubic: 2 problem(s) found:
     Run did not finish in C44/strain_0.01/OUTCAR: no 'General timing and accounting' footer ...
     OUTCAR not found in C11_C12_I/strain_-0.03

Only ``OUTCAR`` is read; ``OSZICAR`` is not used. Each OUTCAR is read once.

Energy and fit
--------------

* The energy is ``energy(sigma->0)`` of the final ionic step.
* The reference volume ``V0`` is the cell volume of the zero-strain run (or of the
  unstrained POSCAR one level up if there is no zero-strain folder).
* In every deformation folder :math:`(E - E_\mathrm{min})/V_0` is fitted with a quadratic
  in the strain :math:`\delta`; the coefficient of :math:`\delta^2`, converted from eV/Å³ to
  GPa (× 160.2177), is the curvature :math:`a`.
  ``Result.txt`` (strain, energy) and ``plot.jpg`` (data, fit and R²) are written into the
  folder. Identical energies (a flat curve) trigger a warning.

From curvatures to constants
----------------------------

Each deformation gives one equation
:math:`a = \tfrac12\,\mathbf e^{\mathsf T} C\,\mathbf e` that is linear in the independent
constants (:math:`\mathbf e` is the strain per unit :math:`\delta` in Voigt notation). The equations are solved together; the methodology section of the output lists
every deformation with its equation.

Output
------

* the fitted curvatures and the methodology;
* the independent constants and the full stiffness matrix (GPa);
* **mechanical stability**: the matrix must be positive definite (all eigenvalues
  positive), which is equivalent to the Born criteria of each crystal system; stable is
  printed in green, unstable in red;
* **polycrystalline averages**: Voigt, Reuss and Hill bulk and shear moduli, Hill Young's
  modulus, Poisson ratio and Pugh ratio (K/G).

The return value is the 6×6 stiffness matrix as a NumPy array in GPa (Voigt notation,
x, y, z of your POSCAR). On the command line, ``--output C.txt`` also saves it as a plain
6×6 text file.

``auto_post`` additionally prints the matrix in the crystal's **standard setting** with the
textbook constants of its crystal system, and can return that one instead
(``frame="standard"``); see :doc:`auto`.

.. admonition:: Next
   :class: tip

   Optional: turn the matrix into directional plots with elastic-vis – :doc:`visualise`.
