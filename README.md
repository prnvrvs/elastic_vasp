This package is based on python to create deformed structure and calculate **second-order elastic constants** from strain-energy relationship. This tool is only compatible with DFT package VASP.

Kindly **cite** following articles if you find this package useful. 

1. *Kumar P, Adlakha I. Effect of interstitial hydrogen on elastic behavior of metals: an ab-initio study, Journal of Engineering Materials and Technology.*  

2. *Mishra P, Kumar P, Neelakantan L, Adlakha I . First-principles prediction of electrochemical polarization and mechanical behavior in Mg based intermetallics, computational material science.*




How to use this Package

**Preprocessing**
```
from elastic_vasp import elastic_gen as eg 
from elastic_vasp import constants 
eg.cubic_pre() 
```
The above code for cubic elastic constants is to create three folder **C11_C12_I**, **C11_C12_II**, **C44** each contain bunch of folder corresponding to strain range. You need to set strain level in file strain.dat

Inside each strain folder, New POSCAR file is available which is basically deformed structure of POSCAR. Once preprocessing is done. You need to copy KPOINTS, POTCAR, INCAR to all strain folder inside all three folder. You can use same bash script copy.sh to ease this process.

Once calculation of each folder is done, you can go to post processing

**Post Processing**\
you can run post processing command to print elastic constants

```
from elastic_vasp import elastic_gen as eg  # This will import module
from elastic_vasp import constants 
#======= post-processing =======
C=constants.cubic_post() 
```
Few thing need to keep in mind
1. Use direct coordinate in initial POSCAR
2. Use ISIF=2 only for calculation
3. You should know the symmetry of crystal and based on symmetry before preprocessing , for test suite use different example folder. For finding symmetry of crystal use phonopy or vaspkit

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
Triclinic and Monolclinic is not added yet in this package. In case Future requirement occurs, we will add to this package

**For any contact, visit my homepage**\
https://sites.google.com/view/prnvkmr4


 
