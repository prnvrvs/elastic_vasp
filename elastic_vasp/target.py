"""
Created on Sat Feb 29 16:07:35 2020

@author: pranav
"""
import numpy as np
import os 
import sys
#
#d1=float(sys.argv[1])
#d2=float(sys.argv[2])
#d3=float(sys.argv[3])
#d4=float(sys.argv[4])
#d5=float(sys.argv[5])
#d6=float(sys.argv[6])
def calc_gen(d1,d2,d3,d4,d5,d6):
     with open("POSCAR", 'r') as file_handle:
         # read file content into list
         lines = file_handle.readlines()
     # print list content
     a=[];
     for i in range (2,5):
          b=lines[i].split()
          c=[float(k) for k in b]
          a.append(c)
     a =np.array(a)
     e=np.array([[1+d1,d6/2,d5/2],[d6/2,1+d2,d4/2],[d5/2,d4/2,1+d3]])
     out=np.matmul(a,e)
     head=lines[0:2]
     tail=lines[5:]
     out2=list(out)
     
     MyFile=open('output.txt','a')
     MyFile.writelines(head)
     np.savetxt(MyFile, out,fmt='%10.12f')
     MyFile.writelines(tail)
     MyFile.close() 
     return
     return
