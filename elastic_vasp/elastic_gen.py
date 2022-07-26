"""
Created on Sat Jun  5 19:16:26 2021

@author: pranav
"""
import os
import shutil
import numpy as np
from elastic_vasp.target import calc_gen
def cubic_pre():
    #===========C11_C12_I===============
    con=['C11_C12_I','C11_C12_II','C44']
    no_const=len(con)
    for folder in con:
     file='output.txt'
     cwd=os.getcwd()
     li=os.listdir()
     if folder in li:
         shutil.rmtree(cwd+'/'+folder)
     if file in li:   
         os.remove('output.txt')
     os.mkdir(folder)
     with open('strain.dat') as f:
         st = f.readlines()
     f.close()
     st[0]= st[0].strip().rstrip("\n")
     strain=st[0].split(",")
     for i in range(0,len(strain)):
       num=float(strain[i])
       os.mkdir(folder+'/'+'strain_'+strain[i])
       if folder=='C11_C12_I':
           calc_gen(num,-1*num,0.0,0.0,0.0,0.0)
       elif folder=='C11_C12_II':
           calc_gen(num,num,num,0.0,0.0,0.0)
       elif folder=='C44':
           calc_gen(0,0,0,num,num,num)
    
       shutil.copy(cwd+"/output.txt",cwd+"/"+folder+"/"+"strain_"+strain[i]+"/POSCAR")
       os.remove('output.txt')
    return
def orthogonal_pre():
    #===========C11_C12_I===============
    con=['C11_C22_C12_I','C11_C22_C12_II','C22_C33_C23_I','C22_C33_C23_II','C11_C33_C13_I','C11_C33_C13_II','C44','C55','C66']
    no_const=len(con)
    for folder in con:
     file='output.txt'
     cwd=os.getcwd()
     li=os.listdir()
     if folder in li:
         shutil.rmtree(cwd+'/'+folder)
     if file in li:
         os.remove('output.txt')
     os.mkdir(folder)
     with open('strain.dat') as f:
         st = f.readlines()
     f.close()
     st[0]= st[0].strip().rstrip("\n")
     strain=st[0].split(",")
     for i in range(0,len(strain)):
       num=float(strain[i])
       num2=(num**2)/(1-num**2)
       num3=(num**2)/(4-num**2)


       os.mkdir(folder+'/'+'strain_'+strain[i])
       if folder=='C11_C22_C12_I':
           calc_gen(num,-1*num,num2,0.0,0.0,0.0)
       elif folder=='C11_C22_C12_II':
           calc_gen(num,num,0,0.0,0.0,0.0)
       if folder=='C22_C33_C23_I':
           calc_gen(num2,num,-1*num,0.0,0.0,0.0)
       elif folder=='C22_C33_C23_II':
           calc_gen(0,num,num,0.0,0.0,0.0)
       if folder=='C11_C33_C13_I':
           calc_gen(num,num2,-1*num,0.0,0.0,0.0)
       elif folder=='C11_C33_C13_II':
           calc_gen(num,0,num,0.0,0.0,0.0)
       elif folder=='C44':
           calc_gen(num3,0,0,num,0,0)
       elif folder=='C55':
           calc_gen(0,num3,0,0,num,0)
       elif folder=='C66':
           calc_gen(0,0,num3,0,0,num)

       shutil.copy(cwd+"/output.txt",cwd+"/"+folder+"/"+"strain_"+strain[i]+"/POSCAR")
       os.remove('output.txt')
    return

def hexagonal_pre():
    #===========C11_C12_I===============
    con=['C11_C12_I','C11_C12_II','C11_C33_C13_I','C44','C11_C12_C13_C33']
    no_const=len(con)
    for folder in con:
     file='output.txt'
     cwd=os.getcwd()
     li=os.listdir()
     if folder in li:
         shutil.rmtree(cwd+'/'+folder)
     if file in li:
         os.remove('output.txt')
     os.mkdir(folder)
     with open('strain.dat') as f:
         st = f.readlines()
     f.close()
     st[0]= st[0].strip().rstrip("\n")
     strain=st[0].split(",")
     for i in range(0,len(strain)):
       num=float(strain[i])
       num2=(num**2)/(1-num**2)
       os.mkdir(folder+'/'+'strain_'+strain[i])
       if folder=='C11_C12_I':
           calc_gen(num,num,0.0,0.0,0.0,0.0)
       elif folder=='C11_C12_II':
           calc_gen(num,-1*num,num2,0.0,0.0,0.0)
       elif folder=='C11_C33_C13_I':
           calc_gen(num,num2,-num,0.0,0.0,0.0)
       elif folder=='C44':
           calc_gen(0,0,0,num,num,0)
       elif folder=='C11_C12_C13_C33':
           calc_gen(num,num,num,0.0,0.0,0.0)

       shutil.copy(cwd+"/output.txt",cwd+"/"+folder+"/"+"strain_"+strain[i]+"/POSCAR")
       os.remove('output.txt')                                      
    return 
def trigonal_1_pre():
    
    con=['C11_C12_I','C11_C12_II','C11_C33_C13_I','C11_C33_C13_II','C11_C44_C14_I','C11_C44_C14_II']
    no_const=len(con)
    for folder in con:
     file='output.txt'
     cwd=os.getcwd()
     li=os.listdir()
     if folder in li:
         shutil.rmtree(cwd+'/'+folder)
     if file in li:
         os.remove('output.txt')
     os.mkdir(folder)
     with open('strain.dat') as f:
         st = f.readlines()
     f.close()
     st[0]= st[0].strip().rstrip("\n")
     strain=st[0].split(",")
     for i in range(0,len(strain)):
       num=float(strain[i])
       num2=(num**2)/(1-num**2)
       os.mkdir(folder+'/'+'strain_'+strain[i])
       if folder=='C11_C12_I':
           calc_gen(num,num,0.0,0.0,0.0,0.0)
       elif folder=='C11_C12_II':
           calc_gen(num,-1*num,num2,0.0,0.0,0.0)
       elif folder=='C11_C33_C13_I':
           calc_gen(num,0,num,0.0,0.0,0.0)
       elif folder=='C11_C33_C13_II':
           calc_gen(num,num2,-1*num,0.0,0.0,0.0)
       elif folder=='C11_C44_C14_I':
           calc_gen(num,0,0,num,0,0)
       elif folder=='C11_C44_C14_II':
           calc_gen(num,0,0,-1*num,0.0,0.0)

       shutil.copy(cwd+"/output.txt",cwd+"/"+folder+"/"+"strain_"+strain[i]+"/POSCAR")
       os.remove('output.txt')
    return
    
def trigonal_2_pre():
    
    con=['C11_C12_I','C11_C12_II','C11_C33_C13_I','C11_C33_C13_II','C11_C44_C14_I','C11_C44_C14_II','C11_C44_C15']
    no_const=len(con)
    for folder in con:
     file='output.txt'
     cwd=os.getcwd()
     li=os.listdir()
     if folder in li:
         shutil.rmtree(cwd+'/'+folder)
     if file in li:
         os.remove('output.txt')
     os.mkdir(folder)
     with open('strain.dat') as f:
         st = f.readlines()
     f.close()
     st[0]= st[0].strip().rstrip("\n")
     strain=st[0].split(",")
     for i in range(0,len(strain)):
       num=float(strain[i])
       num2=(num**2)/(1-num**2)
       os.mkdir(folder+'/'+'strain_'+strain[i])
       if folder=='C11_C12_I':
           calc_gen(num,num,0.0,0.0,0.0,0.0)
       elif folder=='C11_C12_II':
           calc_gen(num,-1*num,num2,0.0,0.0,0.0)
       elif folder=='C11_C33_C13_I':
           calc_gen(num,0,num,0.0,0.0,0.0)
       elif folder=='C11_C33_C13_II':
           calc_gen(num,num2,-1*num,0.0,0.0,0.0)
       elif folder=='C11_C44_C14_I':
           calc_gen(num,0,0,num,0,0)
       elif folder=='C11_C44_C14_II':
           calc_gen(num,0,0,-1*num,0,0)
       elif folder=='C11_C44_C15':
           calc_gen(num,0,0,0,-1*num,0)

       shutil.copy(cwd+"/output.txt",cwd+"/"+folder+"/"+"strain_"+strain[i]+"/POSCAR")
       os.remove('output.txt')
    return
    
def monoclinic_pre():
    #===========C11_C12_I===============
    con=['C11_C22_C12_I','C11_C22_C12_II','C22_C33_C23_I','C22_C33_C23_II','C11_C33_C13_I','C11_C33_C13_II','C44','C11_C55_C15_I','C11_C55_C15_II','C22_C55_C25','C33_C55_C35','C44_C66_C46_I','C44_C66_C46_II']
    no_const=len(con)
    for folder in con:
     file='output.txt'
     cwd=os.getcwd()
     li=os.listdir()
     if folder in li:
         shutil.rmtree(cwd+'/'+folder)
     if file in li:
         os.remove('output.txt')
     os.mkdir(folder)
     with open('strain.dat') as f:
         st = f.readlines()
     f.close()
     st[0]= st[0].strip().rstrip("\n")
     strain=st[0].split(",")
     for i in range(0,len(strain)):
       num=float(strain[i])
       num2=(num**2)/(1-num**2)
       num3=(num**2)/(4-num**2)


       os.mkdir(folder+'/'+'strain_'+strain[i])
       if folder=='C11_C22_C12_I':
           calc_gen(num,-1*num,num2,0.0,0.0,0.0)
       elif folder=='C11_C22_C12_II':
           calc_gen(num,num,0,0.0,0.0,0.0)
       elif folder=='C22_C33_C23_I':
           calc_gen(num2,num,-1*num,0.0,0.0,0.0)
       elif folder=='C22_C33_C23_II':
           calc_gen(0,num,num,0.0,0.0,0.0)
       elif folder=='C11_C33_C13_I':
           calc_gen(num,num2,-1*num,0.0,0.0,0.0)
       elif folder=='C11_C33_C13_II':
           calc_gen(num,0,num,0.0,0.0,0.0)
       elif folder=='C44':
           calc_gen(num3,0,0,num,0,0)
       elif folder=='C11_C55_C15_I':
           calc_gen(num,0,0,0,num,0)
       elif folder=='C11_C55_C15_II':
           calc_gen(num,0,0,0,-1*num,0)
       elif folder=='C22_C55_C25':
           calc_gen(0,num,0,0,num,0)
       elif folder=='C33_C55_C35':
           calc_gen(0,0,num,0,num,0)
       elif folder=='C44_C66_C46_I':
           calc_gen(0,0,0,num,0,num)
       elif folder=='C44_C66_C46_II':
           calc_gen(0,0,0,0,num,-1*num)  

       shutil.copy(cwd+"/output.txt",cwd+"/"+folder+"/"+"strain_"+strain[i]+"/POSCAR")
       os.remove('output.txt')
    return
        


def triclinic_pre():
    #===========C11_C12_I===============
    con=['C11','C11_C22_C12','C22','C22_C33_C23','C33','C11_C33_C13','C44','C55','C66','C11_C44_C14','C11_C55_C15','C11_C66_C16','C22_C44_C24','C22_C55_C25','C22_C66_C26','C33_C44_C34','C33_C55_C35','C33_C66_C36','C44_C55_C45','C44_C66_C46','C55_C66_C56']
    no_const=len(con)
    for folder in con:
     file='output.txt'
     cwd=os.getcwd()
     li=os.listdir()
     if folder in li:
         shutil.rmtree(cwd+'/'+folder)
     if file in li:
         os.remove('output.txt')
     os.mkdir(folder)
     with open('strain.dat') as f:
         st = f.readlines()
     f.close()
     st[0]= st[0].strip().rstrip("\n")
     strain=st[0].split(",")
     for i in range(0,len(strain)):
       num=float(strain[i])
       num2=(num**2)/(1-num**2)
       num3=(num**2)/(4-num**2)


       os.mkdir(folder+'/'+'strain_'+strain[i])
       if folder=='C11':
           calc_gen(num,0,0,0.0,0.0,0.0)
       elif folder=='C11_C22_C12':
           calc_gen(num,num,0,0.0,0.0,0.0)
       elif folder=='C22':
           calc_gen(0,num,0,0.0,0.0,0.0)
       elif folder=='C22_C33_C23':
           calc_gen(0,num,num,0.0,0.0,0.0)
       elif folder=='C33':
           calc_gen(0,0,num,0.0,0.0,0.0)
       elif folder=='C11_C33_C13':
           calc_gen(num,0,num,0.0,0.0,0.0)
       elif folder=='C44':
           calc_gen(0,0,0,num,0,0)
       elif folder=='C55':
           calc_gen(0,0,0,0,num,0)
       elif folder=='C66':
           calc_gen(0,0,0,0,0,num)
       elif folder=='C11_C44_C14':
           calc_gen(num,0,0,num,0,0)
       elif folder=='C11_C55_C15':
           calc_gen(num,0,0,0,num,0)
       elif folder=='C11_C66_C16':
           calc_gen(num,0,0,0,0,num)
       elif folder=='C22_C44_C24':
           calc_gen(0,num,0,num,0,0)
       elif folder=='C22_C55_C25':
           calc_gen(0,num,0,0,num,0)
       elif folder=='C22_C66_C26':
           calc_gen(0,num,0,0,0,num)
       elif folder=='C33_C44_C34':
           calc_gen(0,0,num,num,0,0)
       elif folder=='C33_C55_C35':
           calc_gen(0,0,num,0,num,0)
       elif folder=='C33_C66_C36':
           calc_gen(0,0,num,0,0,num)
       elif folder=='C44_C55_C45':
           calc_gen(0,0,0,num,num,0) 
       elif folder=='C44_C66_C46':
           calc_gen(0,0,0,num,0,num)
       elif folder=='C55_C66_C56':
           calc_gen(0,0,0,0,num,num)  

       shutil.copy(cwd+"/output.txt",cwd+"/"+folder+"/"+"strain_"+strain[i]+"/POSCAR")
       os.remove('output.txt')
    return



