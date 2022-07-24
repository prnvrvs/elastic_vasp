#!/bin/bash
export current=`pwd`
for i in `ls -d */`
do
echo $i
cd $i
for j in `ls -d */`
do
cd $j
echo $j
ln -sf $current/INCAR INCAR
ln -sf $current/POTCAR POTCAR
ln -sf $current/KPOINTS KPOINTS
ln -sf $current/run.pbs run.pbs
## call your vasp run script for eg. #qsub run.pbs # mpiexec -np 24 vasp.xx or anything that work for you
qsub run.pbs
cd ..
done
cd ..
done
