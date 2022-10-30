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


def trigonal_2_post():
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
    os.chdir(cwd+"/"+"C11_C44_C15/")
    c11_c44_c15=coeff() 
    A=c11_c12_I
    B=c11_c12_II
    C=c11_c33_c13_I
    D=c11_c33_c13_II
    E=c11_c44_c14_I
    F=c11_c44_c14_II
    G=c11_c44_c15
    C11=(A+B)/2
    C12=(A-B)/2
    C33=(C+D)-C11
    C13=(C-D)/2
    C44=(E+F)-C11
    C14=(E-F)/2 
    C15=G-(C11/2)-(C44)/2
    C66=(C11-C12)/2
    print("second order coefficients are as follow:")

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


def monoclinic_post():
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
    os.chdir(cwd+"/"+"C11_C55_C15_I/")
    c11_c55_c15_I=coeff()
    os.chdir(cwd+"/"+"C11_C55_C15_II/")
    c11_c55_c15_II=coeff()
    os.chdir(cwd+"/"+"C22_C55_C25/")
    c22_c55_c25=coeff()
    os.chdir(cwd+"/"+"C33_C55_C35/")
    c33_c55_c35=coeff()
    os.chdir(cwd+"/"+"C44_C66_C46_I/")
    c44_c66_c46_I=coeff()
    os.chdir(cwd+"/"+"C44_C66_C46_II/")
    c44_c66_c46_II=coeff()
    
    A=c11_c22_c12_I
    B=c11_c22_c12_II
    C=c22_c33_c23_I
    D=c22_c33_c23_II
    E=c11_c33_c13_I
    F=c11_c33_c13_II
    
    G=c44

    H=c11_c55_c15_I
    I=c11_c55_c15_II
    J=c22_c55_c25
    K=c33_c55_c35
    L=c44_c66_c46_I
    M=c44_c66_c46_II
  
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

    C11=AAA
    C22=BBB
    C33=CCC
    C12=BB
    C23=DD
    C13=FF
    C44=GG
    C55=(H+I)-C11
    C15=(I-H)/2
    C25=J-(C22/2)-(C55/2)
    C35=K-(C33/2)-(C55/2)
    C66=(L+M)-C44
    C46=(M-L)/2


    print('''
            ==== stiffness matrix ====''')
    stiffness=np.array([[C11,C12,C13,0,C15,0],[C12,C22,C23,0,C25,0],[C13,C23,C33,0,C35,0],[0,0,0,C44,0,C46],[C15,C25,C25,C35,C55,0],[0,0,0,0,C46,C66]])

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



























    
def triclinic_post():
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
    os.chdir(cwd+"/"+"C44")
    c44=coeff()
    os.chdir(cwd+"/"+"C11_C44_C14/")
    c11_c44_c14=coeff()
    os.chdir(cwd+"/"+"C55/")
    c55=coeff()
    os.chdir(cwd+"/"+"C11_C55_C15/")
    c11_c55_c15=coeff()
    os.chdir(cwd+"/"+"C66")
    c66=coeff()
    os.chdir(cwd+"/"+"C11_C66_C16")
    c11_c66_c16=coeff()
    os.chdir(cwd+"/"+"C22_C44_C24/")
    c22_c44_c24=coeff()
    os.chdir(cwd+"/"+"C22_C55_C25/")
    c22_c55_c25=coeff()
    os.chdir(cwd+"/"+"C22_C66_C26/")
    c22_c66_c26=coeff()
    os.chdir(cwd+"/"+"C33_C44_C34/")
    c33_c44_c34=coeff()
    os.chdir(cwd+"/"+"C33_C55_C35/")
    c33_c55_c35=coeff()
    os.chdir(cwd+"/"+"C33_C66_C36/")
    c33_c66_c36=coeff()
    os.chdir(cwd+"/"+"C44_C55_C45/")
    c44_c55_c45=coeff()
    os.chdir(cwd+"/"+"C44_C66_C46/")
    c44_c66_c46=coeff()
    os.chdir(cwd+"/"+"C55_C66_C56/")
    c55_c66_c56=coeff()
    
    A=c11_c22_c12_I
    B=c11_c22_c12_II
    C=c22_c33_c23_I
    D=c22_c33_c23_II
    E=c11_c33_c13_I
    F=c11_c33_c13_II
    
    G=c44
    H=c55
    I=c66
    J=c11_c44_c14
    K=c11_c55_c15
    L=c11_c66_c16
  
    M=c22_c44_c24
    N=c22_c55_c25
    O=c22_c66_c26
    P=c33_c44_c34
    Q=c33_c55_c35
    R=c33_c66_c36
    
    S=c44_c55_c45
    T=c44_c66_c46
    U=c55_c66_c56

    print('c11_c22_c12_I=',c11_c22_c12_I)
    print('c11_c22_c12_II=',c11_c22_c12_II)
    print('c22_c33_c23_I=',c22_c33_c23_I)
    print('c22_c33_c23_II=',c22_c33_c23_II)
    print('c11_c33_c13_I=',c11_c33_c13_I)
    print('c11_c33_c13_II=',c11_c33_c13_I)
    print('c44=',c44)
    print('c55=',c55)
    print('c66=',c66)
    print('c11_c44_c14=',c11_c44_c14)
    print('c11_c55_c15=',c11_c55_c15)
    print('c11_c66_c16=',c11_c66_c16)
    print('c22_c44_c24=',c22_c44_c24)
    print('c22_c55_c25=',c22_c55_c25)
    print('c22_c66_c26=',c22_c66_c26)
    print('c33_c44_c34=',c33_c44_c34)
    print('c33_c55_c35=',c33_c55_c35)
    print('c33_c66_c36=',c33_c66_c36)
    print('c44_c55_c45=',c44_c55_c45)
    print('c44_c66_c46=',c44_c66_c46)
    print('c55_c66_c56=',c55_c66_c56) 






 
    AA=(A+B)/2
    BB=(B-A)/2
    CC=(D+C)/2
    DD=(D-C)/2
    EE=(E+F)/2
    FF=(F-E)/2
    AAA=AA-CC+EE
    BBB=AA+CC-EE
    CCC=-AA+CC+EE

    C11=AAA
    C22=BBB
    C33=CCC
    C12=BB
    C23=DD
    C13=FF
    C44=2*G
    C55=2*H
    C66=2*I
    C14=J-(C11+C44)/2
    C15=K-(C11+C55)/2
    C16=L-(C11+C66)/2
    C24=M-(C22+C44)/2
    C25=N-(C22+C55)/2
    C26=O-(C22+C66)/2
    C34=P-(C33+C44)/2
    C35=Q-(C33+C55)/2
    C36=R-(C33+C66)/2
    C45=S-(C44+C55)/2
    C46=T-(C44+C66)/2
    C56=U-(C55+C66)/2

    print('''
            ==== stiffness matrix ====''')
    stiffness=np.array([[C11,C12,C13,C14,C15,C16],[C12,C22,C23,C24,C25,C26],[C13,C23,C33,C34,C35,C36],[C14,C24,C34,C44,C45,C46],[C15,C25,C25,C35,C55,C56],[C16,C26,C36,C46,C56,C66]])

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

def triclinic1_post():
    art.loading()
    cwd=os.getcwd()
    os.chdir(cwd+"/"+"C11/")
    c11=coeff()
    os.chdir(cwd+"/"+"C11_C22_C12/")
    c11_c22_c12=coeff()
    os.chdir(cwd+"/"+"C22/")
    c22=coeff()
    os.chdir(cwd+"/"+"C22_C33_C23/")
    c22_c33_c23=coeff()
    os.chdir(cwd+"/"+"C33/")
    c33=coeff()
    os.chdir(cwd+"/"+"C11_C33_C13/")
    c11_c33_c13=coeff()
    os.chdir(cwd+"/"+"C44")
    c44=coeff()
    os.chdir(cwd+"/"+"C11_C44_C14/")
    c11_c44_c14=coeff()
    os.chdir(cwd+"/"+"C55/")
    c55=coeff()
    os.chdir(cwd+"/"+"C11_C55_C15/")
    c11_c55_c15=coeff()
    os.chdir(cwd+"/"+"C66")
    c66=coeff()
    os.chdir(cwd+"/"+"C11_C66_C16")
    c11_c66_c16=coeff()
    os.chdir(cwd+"/"+"C22_C44_C24/")
    c22_c44_c24=coeff()
    os.chdir(cwd+"/"+"C22_C55_C25/")
    c22_c55_c25=coeff()
    os.chdir(cwd+"/"+"C22_C66_C26/")
    c22_c66_c26=coeff()
    os.chdir(cwd+"/"+"C33_C44_C34/")
    c33_c44_c34=coeff()
    os.chdir(cwd+"/"+"C33_C55_C35/")
    c33_c55_c35=coeff()
    os.chdir(cwd+"/"+"C33_C66_C36/")
    c33_c66_c36=coeff()
    os.chdir(cwd+"/"+"C44_C55_C45/")
    c44_c55_c45=coeff()
    os.chdir(cwd+"/"+"C44_C66_C46/")
    c44_c66_c46=coeff()
    os.chdir(cwd+"/"+"C55_C66_C56/")
    c55_c66_c56=coeff()
    
    A=c11
    B=c22    
    C=c33
    D=c11_c22_c12
    E=c11_c33_c13
    F=c22_c33_c23
    
    G=c44
    H=c55
    I=c66
    J=c11_c44_c14
    K=c11_c55_c15
    L=c11_c66_c16
  
    M=c22_c44_c24
    N=c22_c55_c25
    O=c22_c66_c26
    P=c33_c44_c34
    Q=c33_c55_c35
    R=c33_c66_c36
    
    S=c44_c55_c45
    T=c44_c66_c46
    U=c55_c66_c56

    print('c11=',c11)
    print('c22=',c22)
    print('c33=',c33)
    print('c11_c22_c12=',c11_c22_c12)
    print('c11_c33_c13=',c11_c33_c13)
    print('c22_c33_c23=',c22_c33_c23)
    print('c44=',c44)
    print('c55=',c55)
    print('c66=',c66)
    print('c11_c44_c14=',c11_c44_c14)
    print('c11_c55_c15=',c11_c55_c15)
    print('c11_c66_c16=',c11_c66_c16)
    print('c22_c44_c24=',c22_c44_c24)
    print('c22_c55_c25=',c22_c55_c25)
    print('c22_c66_c26=',c22_c66_c26)
    print('c33_c44_c34=',c33_c44_c34)
    print('c33_c55_c35=',c33_c55_c35)
    print('c33_c66_c36=',c33_c66_c36)
    print('c44_c55_c45=',c44_c55_c45)
    print('c44_c66_c46=',c44_c66_c46)
    print('c55_c66_c56=',c55_c66_c56) 




    C11=2*A
    C22=2*B
    C33=2*C
    C12=D-(C11+C22)/2
    C13=E-(C11+C33)/2
    C23=F-(C22+C33)/2
    C44=2*G
    C55=2*H
    C66=2*I
    C14=J-(C11+C44)/2
    C15=K-(C11+C55)/2
    C16=L-(C11+C66)/2
    C24=M-(C22+C44)/2
    C25=N-(C22+C55)/2
    C26=O-(C22+C66)/2
    C34=P-(C33+C44)/2
    C35=Q-(C33+C55)/2
    C36=R-(C33+C66)/2
    C45=S-(C44+C55)/2
    C46=T-(C44+C66)/2
    C56=U-(C55+C66)/2

    print('''
            ==== stiffness matrix ====''')
    stiffness=np.array([[C11,C12,C13,C14,C15,C16],[C12,C22,C23,C24,C25,C26],[C13,C23,C33,C34,C35,C36],[C14,C24,C34,C44,C45,C46],[C15,C25,C25,C35,C55,C56],[C16,C26,C36,C46,C56,C66]])

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
