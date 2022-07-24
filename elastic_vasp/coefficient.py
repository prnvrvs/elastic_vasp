import os
import numpy as np
from glob import glob
import matplotlib.pyplot as plt
def coeff():
  d=glob("*/") 
  p=[]
  for i in range(0,len(d)):
      s=d[i]
      p.append(s[s.find('_')+1:s.find('/')])
#===============================================
      st="grep 'F' "+d[i]+"OSZICAR|tail -n 1 |awk '{print $5}'>>log"
      os.system(st)
  DataIn = np.loadtxt('log')
  os.system('rm log')
  data=DataIn-np.amin(DataIn)
#===============================================
  st2="grep 'volume of cell' strain_0.00/OUTCAR |tail -n 1|awk '{print $5}'>log"
  os.system(st2)
  Vol = np.loadtxt('log')
  os.system('rm log')

#===============================================
  pp=[float(k) for k in p]
  out=np.stack((np.array(pp), data), axis=-1)
  out2 = out[out[:,0].argsort()]
  np.savetxt('Result.txt', out2,fmt='%5.7f')
  x=out2[:,0]
  y=out2[:,1]/Vol
  z = np.polyfit(x,y,2)
#=============================================
  constant=z[0]*160.2176621
#==============================================
  poly=np.poly1d(z)
  new_x = np.linspace(x[0], x[-1])
  new_y = poly(new_x)
  valeur_T=x
  valeur_min=y
  yhat = poly(valeur_T)
  ybar = sum(valeur_min)/len(valeur_min)
  SST = sum((valeur_min - ybar)**2)
  SSreg = sum((yhat - ybar)**2)
  R2 = SSreg/SST
  image = 'plot.jpg'
  lege = 'R^2='+str(R2)
 #=================================
  plt.plot(x, y, "o", new_x, new_y)
  plt.legend([lege])
  plt.savefig(image)
  plt.close()
  return constant
