Generate and run
================

With the workflow chosen and the inputs ready, generate the strained structures, run VASP
in each of them, and come back for the analysis.

1. Generate the strained structures
-----------------------------------

.. code-block:: bash

   elastic-vasp pre hexagonal                    # strains from strain.dat
   elastic-vasp pre hexagonal --max-strain 0.02 --points 5

or in Python, ``eg.hexagonal_pre()``. This

* checks POSCAR, the strain values and the crystal symmetry – and stops **before touching
  anything** if POSCAR or the strains are wrong;
* warns (and names a working alternative) if the workflow does not fit the crystal as
  oriented;
* refuses to overwrite existing deformation folders, which may hold finished runs, unless
  you add ``--overwrite``;
* writes one folder per deformation and one ``strain_<s>`` folder per strain value, each
  with its ``POSCAR``:

.. code-block:: text

   C11_C12_I/strain_-0.02/POSCAR
   C11_C12_I/strain_-0.01/POSCAR
   ...
   C44/strain_0.02/POSCAR

2. Run VASP
-----------

Copy (or link) the same ``INCAR``, ``KPOINTS`` and ``POTCAR`` into every ``strain_*``
folder and run VASP there. The ``copy.sh`` script in each example folder does this; a
minimal loop:

.. code-block:: bash

   for d in */strain_*/; do (cd "$d" && ln -sf ../../INCAR ../../KPOINTS ../../POTCAR . && mpirun -np 4 vasp_std > vasp.out); done

Settings that matter
--------------------

* ``ISIF = 2`` – the cell shape is the imposed strain; only the atoms relax. (Use
  ``NSW = 0`` / ``IBRION = -1`` for clamped-ion constants.)
* The **same** ``ENCUT``, k-mesh and smearing as for relaxing the POSCAR, and a reference
  POSCAR relaxed with exactly these settings.
* A **Γ-centred** k-mesh: it is defined relative to the crystal, so every strained cell is
  sampled consistently.

How accurate do the energies need to be?
----------------------------------------

The method measures small energy changes. For a strain :math:`\delta` the change per cell
is about :math:`V C \delta^2`: for a 4-atom fcc Al cell (66 Å³) and C11 − C12 ≈ 48 GPa,
that is only **2 meV at δ = 1 %** and 8 meV at 2 %. Numerical noise of the same size makes the fitted
curvatures meaningless. Signs of trouble in the ``Result.txt`` / ``plot.jpg`` of a
deformation folder:

* the energy at zero strain is *not* the lowest point (a W-shaped instead of a U-shaped
  curve);
* the curve is clearly lopsided, or R² is far from 1;
* the constants change a lot with the k-mesh or the strain range.

How to get clean curves:

* **Converge the k-mesh for energy differences of ~0.1 meV per cell**, not just for the
  total energy. Metals need far more k-points than insulators: a 12×12×12 mesh with
  Methfessel–Paxton smearing still left ~1 meV of noise for the 4-atom Al cell above, and a
  single Γ point gives nonsense for metals (strain splits degenerate states at the Fermi
  level and *lowers* the energy).
* **Use larger strains** for soft or noisy cases: ±3–4 % increases the signal 4–16 times
  compared with ±1 %, while anharmonic errors stay at a few tenths of a percent.
* **For simple metals whose atoms do not move under strain** (e.g. fcc/bcc with one atom
  per primitive cell), ``ISMEAR = -5`` (tetrahedron method with Blöchl corrections) gives
  much smoother total energies at the same k-mesh. Do not use it when atoms relax inside
  the strained cell (forces are not accurate with it in metals) – use ``ISMEAR = 1`` with a
  denser mesh instead.

  The :doc:`../worked_example` shows this in practice for Al.
* Tight electronic convergence: ``EDIFF = 1E-7`` or smaller.

VASP builds on macOS: if VASP crashes with a segmentation fault at the first electronic
step and raising the stack limit (``ulimit -s``) is not possible, try
``LSCALAPACK = .FALSE.`` – a ScaLAPACK problem in the build is a common cause.

.. admonition:: Next
   :class: tip

   When all runs have finished: :doc:`postprocessing`.
