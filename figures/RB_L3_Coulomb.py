import math
import numpy as np
import matplotlib.pyplot as plt

## 6206 DEEP GROOVE BALL BEARING
D = 62      
d = 30
dm = (D+d)/2
C0 = 11300

def Coulomb(cof,d,F):
    Mt = cof*F*d/2
    return Mt

## INFLUENCE LOAD
Fr = 2000
Fa = 1000
cofmin = 0.001
cofmax = 0.0015
cof = (cofmin+cofmax)/2
cofv = np.linspace(0,cofmax,100)

## OPERATING CONDITIONS
n = 3000                    ## SPEED / rpm
F = np.sqrt(Fr**2+Fa**2)    ## LOAD F = (Fr^2 + Fa^2)^0.5
vec = np.linspace(0,5000,1000)    

MF = Coulomb(cofmin,d,vec)
MC = Coulomb(cofv,d,F)

M = Coulomb(cof,d,F)

## FIGURE
plt.figure(1)
plt.plot(vec,MF,'k',label=r'$M$')
plt.xlabel(r'Load / N ')
plt.ylabel(r'$M$ / Nmm')
plt.ylim(0,int(10*math.ceil(max(MF)/10)))
plt.grid()
plt.legend()
plt.savefig('pdf/CoulombLoad.pdf', bbox_inches = 'tight',
    pad_inches = 0.1, transparent=True)

plt.figure(2)
plt.plot(cofv,MC,'k',label=r'$M$')
plt.xlabel(r'COF')
plt.ylabel(r'$M$ / Nmm')
plt.ylim(0,int(10*math.ceil(max(MC)/10)))
plt.grid()
plt.legend()
plt.savefig('pdf/CoulombCOF.pdf', bbox_inches = 'tight',
    pad_inches = 0.1, transparent=True)