import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import scienceplots

plt.style.use(['science','ieee'])

## DIMENSIONS
roh = 7850
nu = 0.3
g = 9.81
n = 15000
omega = n*np.pi/30

b = 100/1000
a = 25/1000

r = np.linspace(a,b,100)
sigmar = 1e-6*(3+nu)*g*roh*omega**2/(8*g)*(a**2+b**2-a**2*b**2/r**2-r**2)
sigmat = 1e-6*(3+nu)*g*roh*omega**2/(8*g)*(a**2+b**2+a**2*b**2/r**2-(1+3*nu)/(3+nu)*r**2)

rp = r*1000

## FEM
dfL = pd.read_csv('dados/rotatingcylinder.csv')
xL = dfL[dfL.columns[16]]
sigmatFEM = dfL[dfL.columns[4]]
sigmarFEM = dfL[dfL.columns[3]]


plt.figure(1,figsize=(3,3.2))
plt.subplot(211)
plt.plot(rp,sigmat,'k',label='Analytical')
plt.plot(xL,sigmatFEM,'b-.',label='FEM')
plt.ylabel(r'$\sigma_t$ / MPa')
plt.ylim((0,200))
plt.xlim((25,100))
#plt.grid()
plt.xticks(np.arange(min(rp), max(rp)+1, 5.))
plt.legend(loc='best')
plt.subplot(212)
plt.plot(rp,sigmar,'k',label='Analytical')
plt.plot(xL,sigmarFEM,'b-.',label='FEM')
plt.xlim((25,100))
plt.ylim((0,60))
plt.ylabel(r'$\sigma_r$ / MPa')
plt.xlabel('r / mm')
plt.legend(loc='best')
#plt.grid()
plt.xticks(np.arange(min(rp), max(rp)+1, 5.))
plt.savefig('pdf/rotatingcylinderFEM.pdf', bbox_inches = 'tight',
    pad_inches = 0.1, transparent=True)