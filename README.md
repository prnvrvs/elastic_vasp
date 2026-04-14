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
4. Post-processing auto-detects the zero-strain folder from the `strain_*` names, so it no longer depends on `strain_0.00` specifically.

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


**For any contact, visit my homepage**\
https://sites.google.com/view/prnvkmr4
