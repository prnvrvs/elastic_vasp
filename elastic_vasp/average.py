# calculate average voigt rauss Hill average
import numpy as np
def avg(C):
    C11=C[0,0]
    C22=C[1,1]
    C33=C[2,2]
    C44=C[3,3]
    C55=C[4,4]
    C66=C[5,5]
    C12=C[0,1]
    C13=C[0,2]
    C23=C[1,2]
    K_voigt=((C11+C22+C33)+2*(C12+C13+C23))/9
    G_voigt=((C11+C22+C33)-(C12+C13+C23)+3*(C44+C55+C66))/15
    S=np.linalg.inv(C)
    S11=S[0,0]
    S22=S[1,1]
    S33=S[2,2]
    S44=S[3,3]
    S55=S[4,4]
    S66=S[5,5]
    S12=S[0,1]
    S13=S[0,2]
    S23=S[1,2]

    K_reuss=1/((S11+S22+S33)+2*(S12+S13+S23))
    G_reuss=15/(4*(S11+S22+S33)-4*(S12+S13+S23)+3*(S44+S55+S66))
    K_VRH=(K_voigt+K_reuss)/2
    G_VRH=(G_voigt+G_reuss)/2
    E_VRH=(9*K_VRH*G_VRH)/(G_VRH+3*K_VRH)
    nu_VRH=(3*K_VRH-2*G_VRH)/(6*K_VRH+2*G_VRH)
    print("===== Polycrystalline Average =====")
    print("K_voigt=",'%.2f' %K_voigt,' GPa')
    print("K_reuss=",'%.2f' %K_reuss,' GPa')
    print("K_VRH=",'%.2f' %K_VRH,' GPa')

    print("\nG_voigt=",'%.2f' %G_voigt,' GPa')
    print("G_reuss=",'%.2f' %G_reuss,' GPa')
    print("G_VRH=",'%.2f' %G_VRH,' GPa')
    print("\nE_VRH=",'%.2f' %E_VRH,' GPa')
    print("nu_VRH=",'%.2f' %nu_VRH)
    print("Pugh ratio=",'%.2f' %(K_VRH/G_VRH))
    return 
