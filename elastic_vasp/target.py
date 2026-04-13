"""
Created on Sat Feb 29 16:07:35 2020

@author: pranav
"""
import numpy as np
#
#d1=float(sys.argv[1])
#d2=float(sys.argv[2])
#d3=float(sys.argv[3])
#d4=float(sys.argv[4])
#d5=float(sys.argv[5])
#d6=float(sys.argv[6])
def calc_gen(d1, d2, d3, d4, d5, d6, poscar_path="POSCAR", output_path="output.txt"):
     with open(poscar_path, "r") as file_handle:
         lines = file_handle.readlines()

     lattice = []
     for i in range(2, 5):
          row = [float(k) for k in lines[i].split()]
          lattice.append(row)
     lattice = np.array(lattice)

     strain_matrix = np.array(
         [[1 + d1, d6 / 2, d5 / 2],
          [d6 / 2, 1 + d2, d4 / 2],
          [d5 / 2, d4 / 2, 1 + d3]]
     )
     out = np.matmul(lattice, strain_matrix)

     head = lines[0:2]
     tail = lines[5:]

     with open(output_path, "w") as my_file:
         my_file.writelines(head)
         np.savetxt(my_file, out, fmt="%10.12f")
         my_file.writelines(tail)
     return
