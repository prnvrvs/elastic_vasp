from elastic_vasp import elastic_gen as eg  # This will import module
from elastic_vasp import constants 
# This will generate folder of calculation through vasp
#======= pre-processing ========
#eg.orthogonal_pre()   # uncoment during preprocessing and comment during post processing

# This will create folder for different distorion and range of strain folder inside as mentioned in strain.dat file
# copy INCAR POTCAR KPOINTS to each strain-folder inside each folder and do calculation for each strain-folder

#for copy as well as run you can use copy.sh file edit and run as "bash -i copy.sh"

#======= post-processing =======
#C=constants.orthogonal_post() # uncoment during postprocessing and comment during pre-processing

#===OUTPUT is printed===

