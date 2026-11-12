Prepare the inputs
==================

POSCAR
------

The structure to deform. It should be **fully relaxed** (cell and atoms) with the same
settings you use for the strained runs, so that the reference is stress-free.

Accepted forms (read with ASE):

* VASP 5 format with a species line. POTCAR-style labels such as ``Fe_pv`` or ``O_s`` are
  fine. VASP 4 files without a species line work only when a ``POTCAR`` is next to them.
* **Direct or Cartesian** coordinates. Atoms always move with the strained cell.
* Any **scale factor**: a positive number, a negative number (target volume) or three
  per-axis values (VASP 6).
* Selective dynamics flags.

The strained POSCARs are always written with scale 1.0 and Direct coordinates. The title
line, the atom order and the selective-dynamics flags are kept, so your ``POTCAR`` and
``MAGMOM`` still match. Your own POSCAR is never modified.

.. note::

   For the named functions the crystal must be oriented the way the deformations expect,
   e.g. the hexagonal c axis along z. See :doc:`orientation`, or use :doc:`auto`.

Strain values
-------------

The strain values :math:`\delta` to apply, as fractions (``0.01`` = 1 %), can be passed
directly:

.. code-block:: python

   eg.cubic_pre(strains=[-0.02, -0.01, 0, 0.01, 0.02])
   eg.cubic_pre(strains=eg.symmetric_strains(0.03, 7))   # -0.03 ... 0.03 in 7 steps

or on the command line with ``--strains -0.02 -0.01 0 0.01 0.02`` or
``--max-strain 0.03 --points 7`` (see :doc:`cli`). Every ``*_pre`` function, including
``auto_pre``, takes ``strains=``.

strain.dat
~~~~~~~~~~

If no strains are given, they are read from ``strain.dat``:

.. code-block:: text

   # comments start with '#'
   -0.03, -0.02, -0.01, 0.00
   0.01 0.02 0.03

* Values may be separated by commas and/or spaces, over any number of lines.
* At least **3 distinct** values are needed for the quadratic fit; 5–7 values
  symmetric around zero are recommended.
* ``0``, ``0.0``, ``-0.00`` all give the folder ``strain_0.00``. Repeated values and
  typos stop pre-processing before any folder is touched.
* If the file is missing, it is created with
  ``-0.03,-0.02,-0.01,0.00,0.01,0.02,0.03`` and a note is printed.

Choosing the range: ±1 % to ±3 % is typical. Smaller strains suffer more from numerical
noise in the energies; larger strains include more anharmonic (higher-order) effects. With
the default ±3 % the anharmonic error is a few tenths of a percent for typical materials.

The zero-strain folder is optional. Its volume is used as the reference volume ``V0`` if
present; otherwise the volume of the unstrained ``POSCAR`` one level up is used.

INCAR, KPOINTS, POTCAR
----------------------

Not read by elastic_vasp. Copy them into every strain folder (see the ``copy.sh`` script
in the examples). Use ``ISIF = 2``; see :doc:`workflow` for the other settings.

.. admonition:: Next
   :class: tip

   :doc:`workflow` – write the strained structures and run VASP.
