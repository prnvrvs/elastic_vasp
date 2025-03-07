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
2. *Mishra P, Kumar P, Neelakantan L, Adlakha I. First-principles prediction of electrochemical polarization and mechanical behavior in Mg based intermetallics. Computational Materials Science. 2022 Nov 1;214:111667. https://doi.org/10.1016/j.commatsci.2022.111667* 
```
@article{MISHRA2022111667,
title = {First-principles prediction of electrochemical polarization and mechanical behavior in Mg based intermetallics},
journal = {Computational Materials Science},
volume = {214},
pages = {111667},
year = {2022},
issn = {0927-0256},
doi = {https://doi.org/10.1016/j.commatsci.2022.111667},
author = {Pragyandipta Mishra and Pranav Kumar and Lakshman Neelakantan and Ilaksh Adlakha},
}
```


**Installation**
```
pip install elastic-vasp
```

How to use this Package

**Preprocessing**
```
from elastic_vasp import elastic_gen as eg 
from elastic_vasp import constants 
eg.cubic_pre() 
```
The above code for cubic elastic constants is to create three folder **C11_C12_I**, **C11_C12_II**, **C44** each contain bunch of folder corresponding to strain range. You need to set strain level in file "strain.dat". (Note: zero strain should be metioned as 0.00 to avoid-post processing error)

Inside each strain folder, New POSCAR file is available which is basically deformed structure of POSCAR. Once preprocessing is done. You need to copy KPOINTS, POTCAR, INCAR to all strain folder inside all three folder. You can use same bash script copy.sh to ease this process.

Once calculation of each folder is done, you can go for post processing

**Post Processing**\
you can run post processing command to print elastic constants

```
from elastic_vasp import elastic_gen as eg  # This will import module
from elastic_vasp import constants 
#======= post-processing =======
C=constants.cubic_post() 
```
Following option are available:
```
cubic_pre()
orthogonal_pre() ------ for orthorhombic and Tetragonal
hexagonal_pre() 
trigonal_1_pre()
trigonal_2_pre()
monoclinic_pre()
triclinic_pre()
```

Few thing need to keep in mind
1. Use direct coordinate in initial POSCAR
2. Use ISIF=2 only for calculation
3. You should know the symmetry of crystal and based on symmetry before preprocessing , for test suite use different example folder. For finding symmetry of crystal use spglib/phonopy 

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


I will suggest you to find orthorhombic elastic constants if you are dealing with tetragonal system. This will reduce confusion which side is longer and which elastic constants to calculate. Just directly calculate all 9 elastic constants.\
Similarly, for trigonal calculate all 7 elastic constants. 


**For any contact, visit my homepage**\
https://sites.google.com/view/prnvkmr4
