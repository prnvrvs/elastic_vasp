Release notes
=============

The original code (0.0.4) is kept on the ``legacy`` branch; see :doc:`legacy`.

0.1.2
-----

* New video tutorial (`YouTube <https://youtu.be/kXaeIoNdHzw>`_): the full workflow on
  fcc aluminium with real VASP results; the earlier video is kept as the legacy tutorial
  (:doc:`tutorial`).
* ``tutorial_video/`` contains the script that builds the video and its editable narration.

0.1.1
-----

* Works with spglib 2.8: its new deprecation warning about error handling
  (``OLD_ERROR_HANDLING``) is no longer shown, and spglib errors are handled the same way
  on every version from 2.4 to 3.x. This also fixes the tests, which failed with spglib 2.8.

0.1.0
-----

* **Automatic mode for any orientation**: ``auto_pre`` / ``auto_post`` derive the
  independent constants and a well-conditioned set of deformations from the crystal's
  symmetry (spglib) in the POSCAR's own orientation.
* **Tetragonal workflows**: ``tetragonal_1_*`` (4/mmm, 6 constants) and
  ``tetragonal_2_*`` (4, -4, 4/m, 7 constants including C16). Previously tetragonal II
  crystals were pointed to the orthorhombic workflow, which cannot measure C16.
* **Symmetry check in every named** ``*_pre`` **function**: warns when the deformations
  cannot represent the crystal's stiffness matrix (wrong system, setting or orientation)
  and names functions that work.
* Orientation checks use the crystal symmetry from spglib (lattice and atoms), so
  pseudo-symmetric cells are judged correctly.
* The reference volume :math:`V_0` is read from the exact POSCAR of the zero-strain folder;
  OUTCAR prints the cell volume with only two decimals (an error of up to ~0.03 % in every
  constant). OUTCAR is still used when no POSCAR is present.
* New dependency: spglib ≥ 2.4.
* Tests for all 230 space groups.
* ``check_symmetry()`` / ``elastic-vasp check``: space group, constants and which workflow
  works with the structure as oriented.
* Command line ``elastic-vasp`` with ``check``, ``pre``, ``post`` and ``systems``.
* Strains can be passed directly (``strains=``, ``symmetric_strains()``,
  ``--strains`` / ``--max-strain``) instead of through ``strain.dat``.
* Symmetry checks tolerate the ~1e-7 numerical noise of relaxed structures (previously
  such cells could be reported as misoriented).
* Symmetry analysis uses the primitive cell (same orientation), so supercells keep their
  full point group – e.g. an fcc cell with x ∥ [1-10], y ∥ [11-2], z ∥ [111] has 3
  constants, not 13.
* ``auto_post`` also reports the stiffness matrix in the crystal's standard setting with the
  textbook constants (``frame="standard"`` returns it; ``--frame standard`` on the
  command line).
* ``elastic-vasp post --output`` writes a plain 6×6 matrix (no header line), readable
  directly by elastic-vis.
* This documentation.

0.0.9
-----

* CI: test Python 3.10, 3.12 and 3.13 (Python 3.9 has no build for the current
  ``ubuntu-latest``); ``python_requires >= 3.10``.

0.0.8
-----

* POSCAR reading and writing with ASE ≥ 3.23: any scale-factor form, Cartesian
  coordinates (atoms now move with the strained cell), selective dynamics, POTCAR-style
  species labels.
* Pre-processing no longer deletes a user file named ``output.txt``.

0.0.7
-----

* One table per crystal system drives both pre- and post-processing.
* Energies from OUTCAR ``energy(sigma->0)`` only; unfinished and unconverged runs are
  listed and post-processing stops.
* POSCAR and ``strain.dat`` are checked before any folder is touched; ``strain.dat`` may
  span several lines and is created with defaults if missing; at least 3 strains.
* Trigonal example is now α-Al\ :sub:`2`\ O\ :sub:`3`.

0.0.6
-----

* Fixed the trigonal II (C15) and monoclinic (C15, C46, C44_C66_C46_II pattern) constants.
* Orientation checks for monoclinic, hexagonal and trigonal structures.
* Existing deformation folders are not overwritten without ``overwrite=True``.
