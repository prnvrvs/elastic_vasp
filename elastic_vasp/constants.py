"""
Created on Sat Jun  5 19:16:26 2021

@author: pranav
"""
import os
import shutil
import numpy as np
from elastic_vasp.coefficient import coeff
from elastic_vasp import art
from elastic_vasp import average
def cubic_post():
    art.loading() 
    cwd=os.getcwd()
    os.chdir(cwd+"/"+"C11_C12_I/")
    c11_c12_I=coeff()
    os.chdir(cwd+"/"+"C11_C12_II/")
    c11_c12_II=coeff()
    os.chdir(cwd+"/"+"C44/")
    c44=coeff()
    bulk=c11_c12_II*2/9
    shear=c11_c12_I*0.5
    C11=(3*bulk+4*shear)/3
    C12=(3*bulk-2*shear)/3
    C44=2*c44/3
    
    print("second order coefficients are as follow:")
    print("C11_C12_I/=",c11_c12_I,"\nC12/=",c11_c12_II,"\nC44/=",c44)
    print('''
    ===== Methodology =====
    Bulk=(C11_C12_II/) x 2/9
    shear=(C11_C12_I/) x 1/2
    C11=(3xbulk+4xshear)/3
    C12=(3xbulk-2xshear)/3
    C44=2/3 x C44/
           ''')
    print("C11=",C11,"\nC12=",C12,"\nC44=",C44,"\nbulk_modulus=",bulk)
    print('''
    ==== stiffness matrix ====''')
    stiffness=np.array([[C11,C12,C12,0,0,0],[C12,C11,C12,0,0,0],[C12,C12,C11,0,0,0],[0,0,0,C44,0,0],[0,0,0,0,C44,0],[0,0,0,0,0,C44]])
    print(np.around(stiffness,decimals=2))

    print('''
    check stability criteria
    C11−C12 > 0 ; C11+2C12 > 0 ; C44 > 0
         ''')
    #checking stability
    a=C11-C12
    b=C11+2*C12
    c=C44
    if (a>0)and(b>0)and(c>0):
        print("\033[1m"+"\033[1;32m structure is mechanically stable"+"\033[0m")
    else:
        print("\033[1m"+"\033[1;32m structure is mechanically unstable"+"\033[0m")
    average.avg(stiffness)
    art.baginga()
    return stiffness
def hexagonal_post():
    art.loading() 
    cwd=os.getcwd()
    os.chdir(cwd+"/"+"C11_C12_I/")
    c11_c12_I=coeff()
    os.chdir(cwd+"/"+"C11_C12_II/")
    c11_c12_II=coeff()
    os.chdir(cwd+"/"+"C11_C33_C13_I/")
    c11_c33_c13_I=coeff()
    os.chdir(cwd+"/"+"C44/")
    c44=coeff()
    os.chdir(cwd+"/"+"C11_C12_C13_C33/")
    c11_c12_c13_c33=coeff()
    C11=(c11_c12_I+c11_c12_II)/2
    C12=(c11_c12_I-c11_c12_II)/2
    C44=c44
    eq1=2*c11_c33_c13_I-C11
    eq2=2*(c11_c12_c13_c33-c11_c12_I)
    C13=(eq2-eq1)/6
    C33=eq1+2*C13
    print(eq1,eq2)
    bulk=c11_c12_c13_c33*2/9
    print("second order coefficients are as follow:")
    print("C11_C12_I/=",c11_c12_I,"\nC11_C12_II/=",c11_c12_II,"C11_C33_C13_I/=",c11_c33_c13_I)
    print("C44/=",c44,"\nC11_C12_C13_C33/=",c11_c12_c13_c33)
    print('''
    Methodology
    C11=(C11_C12_I/+C11_C12_I/)/2
    C12=(C11_C12_I/-C11_C12_II/)/2
    eq1=2*C11_C33_C13_I/-C11
    eq2=2*(C11_C12_C13_C33/-C11_C12_I/)
    C13=(eq2-eq1)/6
    C33=eq1+2*C13
    C44=C44/
    bulk=c11_c12_c13_c33*2/9
           ''')
    print("C11=",C11,"\nC12=",C12,"\nC33",C33,"\nC13=",C13,"\nC44=",C44,"\nbulk_modulus=",bulk)
    print('''
            ==== stiffness matrix ====''')
    stiffness=np.array([[C11,C12,C13,0,0,0],[C12,C11,C13,0,0,0],[C13,C13,C33,0,0,0],[0,0,0,C44,0,0],[0,0,0,0,C44,0],[0,0,0,0,0,(C11-C12)/2]])
          
    print(np.around(stiffness,decimals=2))

    print('''
        check stability criteria
    C11 >|C12| ; 2C13^2 < C33(C11 + C12)
    C44 > 0 ; C66 > 0 
            ''')
    #checking stability
    
    a=2*C13**2
    b=C33*(C11+C12)
    if (C11>abs(C12))and(C44>0)and((C11-C12)>0)and(a<b):
        print("\033[1m"+"\033[1;32m structure is mechanically stable"+"\033[0m")
    else:
        print("\033[1m"+"\033[1;32m structure is mechanically unstable"+"\033[0m")
    average.avg(stiffness)
    art.baginga()
    return stiffness
def orthogonal_post():
    art.loading()
    cwd=os.getcwd()
    os.chdir(cwd+"/"+"C11_C22_C12_I/")
    c11_c22_c12_I=coeff()
    os.chdir(cwd+"/"+"C11_C22_C12_II/")
    c11_c22_c12_II=coeff()
    os.chdir(cwd+"/"+"C22_C33_C23_I/")
    c22_c33_c23_I=coeff()
    os.chdir(cwd+"/"+"C22_C33_C23_II/")
    c22_c33_c23_II=coeff()
    os.chdir(cwd+"/"+"C11_C33_C13_I/")
    c11_c33_c13_I=coeff()
    os.chdir(cwd+"/"+"C11_C33_C13_II/")
    c11_c33_c13_II=coeff()
    os.chdir(cwd+"/"+"C44/")
    c44=coeff()
    os.chdir(cwd+"/"+"C55/")
    c55=coeff()
    os.chdir(cwd+"/"+"C66/")
    c66=coeff()
    A=c11_c22_c12_I
    B=c11_c22_c12_II
    C=c22_c33_c23_I
    D=c22_c33_c23_II
    E=c11_c33_c13_I
    F=c11_c33_c13_II
    G=c44
    H=c55
    I=c66
    AA=(A+B)/2
    BB=(B-A)/2
    CC=(D+C)/2
    DD=(D-C)/2
    EE=(E+F)/2
    FF=(F-E)/2
    AAA=AA-CC+EE
    BBB=AA+CC-EE
    CCC=-AA+CC+EE
    GG=G*2
    HH=H*2
    II=I*2
    C11=AAA
    C22=BBB
    C33=CCC
    C12=BB
    C23=DD
    C13=FF
    C44=GG
    C55=HH
    C66=II
    print("second order coefficients are as follow:")
    print("C11_C22_C12_I/=",c11_c22_c12_I,"\nC11_C22_C12_II/=",c11_c22_c12_II)
    print("C22_C33_C23_I/=",c22_c33_c23_I,"\nC22_C33_C23_II=",c22_c33_c23_II)
    print("C11_C33_C13_I/=",c11_c33_c13_I,"\nC11_C33_C13_II/=",c11_c33_c13_II)
    print("C44/=",c44,"\nC55/=",c55,"\nC66/=",c66)
    print('''
    Methodology
    A=c11_c22_c12_I/
    B=c11_c22_c12_II/
    C=c22_c33_c23_I/
    D=c22_c33_c23_II/
    E=c11_c33_c13_I/
    F=c11_c33_c13_II
    G=c44/
    H=c55/
    I=c66/
    AA=(A+B)/2
    BB=(B-A)/2
    CC=(D+C)/2
    DD=(D-C)/2
    EE=(E+F)/2
    FF=(F-E)/2
    AAA=AA-CC+EE
    BBB=AA+CC-EE
    CCC=-AA+CC+EE
    GG=G*2
    HH=H*2
    II=I*2
    C11=AAA
    C22=BBB
    C33=CCC
    C12=BB
    C23=DD
    C13=FF
    C44=GG
    C55=HH
    C66=II
        ''')
    print("C11=",C11,"\nC22=",C22,"\nC33=",C33,"\nC12=",C12,"\nC13=",C13,"\nC23=",C23)
    print("C44=",C44,"\nC55=",C55,"\nC66=",C66)
    print('''
            ==== stiffness matrix ====''')
    stiffness=np.array([[C11,C12,C13,0,0,0],[C12,C22,C23,0,0,0],[C13,C23,C33,0,0,0],[0,0,0,C44,0,0],[0,0,0,0,C55,0],[0,0,0,0,0,C66]])

    print(np.around(stiffness,decimals=2))

    print('''
        check stability criteria
    stiffness matrix is positive definite
    all eigen value are positive
            ''')
    #checking stability

    if (np.all(np.linalg.eigvals(stiffness)) > 0):
        print("\033[1m"+"\033[1;32m structure is mechanically stable"+"\033[0m")
    else:
        print("\033[1m"+"\033[1;32m structure is mechanically unstable"+"\033[0m")
    average.avg(stiffness)
    art.baginga()
    return stiffness
def trigonal_1_post():
    art.loading()
    cwd=os.getcwd()
    os.chdir(cwd+"/"+"C11_C12_I/")
    c11_c12_I=coeff()
    os.chdir(cwd+"/"+"C11_C12_II/")
    c11_c12_II=coeff() 
    os.chdir(cwd+"/"+"C11_C33_C13_I/")
    c11_c33_c13_I=coeff()  
    os.chdir(cwd+"/"+"C11_C33_C13_II/")
    c11_c33_c13_II=coeff()
    os.chdir(cwd+"/"+"C11_C44_C14_I/")
    c11_c44_c14_I=coeff() 
    os.chdir(cwd+"/"+"C11_C44_C14_II/")
    c11_c44_c14_II=coeff()  
    A=c11_c12_I
    B=c11_c12_II
    C=c11_c33_c13_I
    D=c11_c33_c13_II
    E=c11_c44_c14_I
    F=c11_c44_c14_II
    C11=(A+B)/2
    C12=(A-B)/2
    C33=(C+D)-C11
    C13=(C-D)/2
    C44=(E+F)-C11
    C14=(E-F)/2 
    C66=(C11-C12)/2
    print("second order coefficients are as follow:")
    print("C11_C12_I/=",c11_c12_I,"\nC11_C12_II/=",c11_c12_II)
    print("C11_C33_C13_I/=",c11_c33_c13_I,"\nC11_C33_C13_II=",c11_c33_c13_II)
    print("C11_C44_C14_I/=",c11_c44_c14_I,"\nC11_C44_C14_II/=",c11_c44_c14_II)
    print('''
    Methodology
    A=c11_c12_I/
    B=c11_c12_II/
    C=c11_c33_c13_I/
    D=c11_c33_c13_II/
    E=c11_c44_c14_I/
    F=c11_c44_c14_II/
    C11=(A+B)/2
    C12=(A-B)/2
    C33=(C+D)-C11
    C13=(C-D)/2
    C44=(E+F)-C11
    C14=(E-F)/2
	''')
    print("C11=",C11,"\nC22=",C11,"\nC33=",C33,"\nC12=",C12,"\nC13=",C13,"\nC23=",C13)
    print("C44=",C44,"\nC55=",C44,"\nC14=",C14,"\nC66=",C66)
    print('''
            ==== stiffness matrix ====''')
    stiffness=np.array([[C11,C12,C13,C14,0,0],[C12,C11,C13,-1*C14,0,0],[C13,C13,C33,0,0,0],[C14,-1*C14,0,C44,0,0],[0,0,0,0,C44,C14],[0,0,0,0,C14,C66]])

    print(np.around(stiffness,decimals=2))

    print('''
        check stability criteria
    stiffness matrix is positive definite
    all eigen value are positive
            ''')
    #checking stability

    if (np.all(np.linalg.eigvals(stiffness)) > 0):
        print("\033[1m"+"\033[1;32m structure is mechanically stable"+"\033[0m")
    else:
        print("\033[1m"+"\033[1;32m structure is mechanically unstable"+"\033[0m")
    average.avg(stiffness)
    art.baginga()

    return stiffness
    

