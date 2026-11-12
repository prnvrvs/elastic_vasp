Quick start
===========

One complete calculation, using the command line. Every step has a Python equivalent,
shown below each block.

You need a folder with a **relaxed** structure and your VASP inputs:

.. code-block:: text

   POSCAR     relaxed structure (VASP 5 format)
   INCAR      with ISIF = 2
   KPOINTS
   POTCAR

1. Check
--------

.. code-block:: bash

   elastic-vasp check

.. code-block:: text

   Space group:        Fm-3m (225) – cubic (Laue class m-3m)
   Elastic constants:  3 independent (C11, C12, C44 in this orientation)
   Recommended:        cubic_pre() / cubic_post()  – works in this orientation

Python: ``from elastic_vasp import elastic_gen as eg; eg.check_symmetry()``

If the recommendation does **not** work in your orientation, use ``auto`` instead of
``cubic`` in the next steps.

2. Generate the strained structures
-----------------------------------

.. code-block:: bash

   elastic-vasp pre cubic --max-strain 0.03 --points 7

This creates ``C11_C12_I/``, ``C11_C12_II/`` and ``C44/``, each with seven ``strain_*``
folders (−3 % … +3 %) holding a ``POSCAR``. Python:
``eg.cubic_pre(strains=eg.symmetric_strains(0.03, 7))``

3. Run VASP
-----------

Copy ``INCAR``, ``KPOINTS`` and ``POTCAR`` into every ``strain_*`` folder and run VASP
there (see :doc:`guide/workflow` for a loop and the settings that matter).

4. Analyse
----------

.. code-block:: bash

   elastic-vasp post cubic --output C.txt

It prints the fitted curvature of every deformation, the constants (C11, C12, C44), the
6×6 stiffness matrix, whether the crystal is mechanically stable, and the Voigt, Reuss and
Hill bulk and shear moduli, Young's modulus and Poisson ratio.

Python: ``from elastic_vasp import constants; C = constants.cubic_post()``

``C.txt`` holds the 6×6 stiffness matrix in GPa.

.. admonition:: Next
   :class: tip

   * See a complete run with real VASP output: :doc:`worked_example`.
   * Before running VASP for real, read the settings advice in :doc:`guide/workflow`.
   * For your own material, follow the :doc:`user guide <guide/check>` step by step.
