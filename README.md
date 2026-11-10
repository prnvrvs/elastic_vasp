This package is based on python to create deformed structure and calculate **second-order elastic constants** from strain-energy relationship. This tool is only compatible with DFT package VASP.

Kindly **cite** following articles if you find this package useful. 

1. *Kumar, P., and Adlakha, I. (July 26, 2022). "Effect of Interstitial Hydrogen on Elastic Behavior of Metals: an Ab-Initio Study." ASME. J. Eng. Mater. Technol. doi: https://doi.org/10.1115/1.4055097* 
```
@article{Elastic_vasp,
    author = {Kumar, P. and Adlakha, I.},
    title = "{Effect of Interstitial Hydrogen on Elastic Behavior of Metals: An Ab-Initio Study}",
    journal = {Journal of Engineering Materials and Technology},
    volume = {145},
    number = {1},
    year = {2022},
    month = {08},
    issn = {0094-4289},
    doi = {10.1115/1.4055097},
}

```

**Installation**
```
pip install elastic-vasp
```

How to use this Package

The repository now keeps the original code line on the `legacy` branch and the current refactored code on `main`.

**Preprocessing**
```
from elastic_vasp import elastic_gen as eg

eg.cubic_pre()
```
The above code for cubic elastic constants creates three folders: **C11_C12_I**, **C11_C12_II**, and **C44**. Each folder contains one subfolder per strain value from `strain.dat`. Zero strain is normalized automatically, so `0`, `0.0`, and `0.00` all work.

Inside each strain folder, a new `POSCAR` is created for the deformed structure. After preprocessing, copy `KPOINTS`, `POTCAR`, and `INCAR` into every strain folder. You can use the bundled `copy.sh` script to automate that step.

After the VASP calculations finish, run post-processing to print the elastic constants.

**Post Processing**\
you can run post processing command to print elastic constants

```
from elastic_vasp import constants 

C = constants.cubic_post()
```
Available post-processing functions:
```
cubic_post()
hexagonal_post()
orthogonal_post()
trigonal_1_post()
trigonal_2_post()
monoclinic_post()
triclinic_post()
```

Available preprocessing functions:
```
cubic_pre()
orthogonal_pre() ------ for orthorhombic and Tetragonal
hexagonal_pre() 
trigonal_1_pre()
trigonal_2_pre()
monoclinic_pre()
triclinic_pre()
```

Things to keep in mind:
1. Use direct coordinate in initial POSCAR
2. Use ISIF=2 only for calculation
3. You should know the symmetry of the crystal before preprocessing. Use the example folders for tests. To determine symmetry, use spglib or phonopy.
4. Preprocessing refuses to run if the deformation folders already exist, so finished calculations are not deleted by accident. Pass `overwrite=True` (e.g. `cubic_pre(overwrite=True)`) to regenerate them.
5. Post-processing auto-detects the zero-strain folder from the `strain_*` names, so it no longer depends on `strain_0.00` specifically.
6. Post-processing reads the energy from `energy(sigma->0)` of the final ionic step in each OUTCAR (OSZICAR is not used). If any strain folder has no OUTCAR, an OUTCAR without that line, or a run that did not finish normally (no `General timing and accounting` footer, e.g. killed by the walltime), or an ionic relaxation (`NSW > 0`, `IBRION` = 1, 2 or 3) that ran out of steps without `reached required accuracy`, it lists every such folder and stops before fitting.
7. `strain.dat` may list strains separated by commas and/or spaces, over several lines (`#` starts a comment). At least 3 values are needed for the quadratic fit. If it is missing, pre-processing creates it with `-0.03,-0.02,-0.01,0.00,0.01,0.02,0.03` and says so. Pre-processing checks POSCAR and strain.dat before deleting or creating any folder.

| Crystal system       | Space-group          |No. of independent elastic constants  |
| ------------- |:-------------:| -----:|
|Triclinic|1-2|21|
|Monoclinic	     |3-15		|	13|
|Orthorhombic	    | 16-74|			9|
|Tetragonal I	    | 89-142|			6|
|Tetragonal II	  |   75-88	|		7|
|Trigonal I	    | 149-167		|	6|
|Trigonal II	  |   143-148	|		7|
|Hexagonal	   |  168-194			|5 |
|Cubic		     |195-230			|3 |


For tetragonal systems, it is usually simpler to calculate orthorhombic elastic constants first. That avoids confusion about which axis is longer and lets you compute all 9 elastic constants directly.

Similarly, for trigonal systems, calculate all 7 elastic constants.

**Note: use Trigonal II (`trigonal_2_pre()` / `trigonal_2_post()`) for all trigonal crystals unless you are sure of the in-plane orientation.** Trigonal I additionally assumes that x lies along a 2-fold axis (normal to a mirror plane). That holds for groups such as P-3m1, R-3m, P321 and R3m with a along x, but not for "31m"-type groups (P-31m, P-31c, P31m, P31c, P312, P3<sub>1</sub>12, P3<sub>2</sub>12). Trigonal II fits C14 and C15 together, so it is correct for any in-plane orientation at the cost of one extra deformation set.

### Required orientation

The strains are applied along the Cartesian x, y, z axes of the POSCAR lattice, so the crystal must be aligned accordingly:

| Function | Required orientation |
| --- | --- |
| `cubic_*` | cube axes along x, y, z (standard conventional or primitive fcc/bcc cells are fine) |
| `orthogonal_*` | a ∥ x, b ∥ y, c ∥ z |
| `hexagonal_*`, `trigonal_2_*` | c (6-fold/3-fold axis) ∥ z, any in-plane orientation |
| `trigonal_1_*` | c ∥ z and x along a 2-fold axis (see note above) |
| `monoclinic_*` | unique axis (b, the 2-fold axis) ∥ y; a and c in the xz plane |
| `triclinic_*` | any |

`hexagonal_pre()`, `trigonal_1_pre()`, `trigonal_2_pre()` and `monoclinic_pre()` check this automatically. This covers rhombohedral cells, e.g. a primitive cell with the 3-fold axis along [111], and c-unique monoclinic cells. If the axis is misaligned they print a warning saying where the axis actually is, write a rotated copy to `POSCAR.reoriented` (a pure rotation: fractional coordinates and atom order are unchanged, so POTCAR and MAGMOM ordering stay valid) and then continue with the current orientation. To use the rotated cell, copy `POSCAR.reoriented` over `POSCAR` and rerun the pre-processing with `overwrite=True`.


**For any contact, visit my homepage**\
https://sites.google.com/view/prnvkmr4
