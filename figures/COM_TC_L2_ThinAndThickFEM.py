import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

## DIMENSIONS

a = 25 #mm
b = 100

r = np.linspace(a,b,1000)

## STRESSES THICK WALL
pi = 1 #MPa
sigmat = pi*a**2*(1+b**2/r**2)/(b**2-a**2)
sigmar = pi*a**2*(1-b**2/r**2)/(b**2-a**2)

## STRESSES THICK WALL
pi = 1 #MPa
po = -1 #MPa
sigmarpo = (a**2*pi-b**2*po)/(b**2-a**2)-(pi-po)*a**2*b**2/(r**2*(b**2-a**2))
sigmatpo = (a**2*pi-b**2*po)/(b**2-a**2)+(pi-po)*a**2*b**2/(r**2*(b**2-a**2))

## FEM
dfL = pd.read_csv('dados/thickwall.csv')
xL = dfL[dfL.columns[16]]
sigmatFEM = dfL[dfL.columns[4]]
sigmarFEM = dfL[dfL.columns[3]]

dfpo = pd.read_csv('dados/thickwallpo.csv')
xpo = dfpo[dfpo.columns[16]]
sigmatFEMpo = dfpo[dfpo.columns[4]]
sigmarFEMpo = dfpo[dfpo.columns[3]]

plt.figure(1,figsize=(4.5,6))
plt.subplot(211)
plt.plot(r,sigmat,'k-',label='Lamé')
plt.plot(xL,sigmatFEM,'b-.',label='FEM')
plt.ylabel(r'$\sigma_t$ / MPa')
plt.ylim((0,1.5))
plt.xlim((25,100))
plt.grid()
plt.xticks(np.arange(min(r), max(r)+1, 5.))
plt.legend(loc='best')
plt.subplot(212)
plt.plot(r,sigmar,'k-',label=r'Lamé')
plt.plot(xL,sigmarFEM,'b-.',label='FEM')
plt.xlim((25,100))
plt.ylabel(r'$\sigma_r$ / MPa')
plt.xlabel('r / mm')
plt.legend(loc='best')
plt.grid()
plt.xticks(np.arange(min(r), max(r)+1, 5.))
plt.savefig('pdf/thickwallFEM.pdf', bbox_inches = 'tight',
    pad_inches = 0.1, transparent=True)

plt.figure(2,figsize=(4.5,6))
plt.subplot(211)
plt.plot(r,sigmatpo,'k-',label='Lamé')
plt.plot(xpo,sigmatFEMpo,'b-.',label='FEM')
plt.ylabel(r'$\sigma_t$ / MPa')
plt.ylim((0,3.26))
plt.xlim((25,100))
plt.grid()
plt.xticks(np.arange(min(r), max(r)+1, 5.))
plt.legend(loc='best')
plt.subplot(212)
plt.plot(r,sigmarpo,'k-',label=r'Lamé')
plt.plot(xpo,sigmarFEMpo,'b-.',label='FEM')
plt.xlim((25,100))
plt.ylabel(r'$\sigma_r$ / MPa')
plt.xlabel('r / mm')
plt.legend(loc='best')
plt.grid()
plt.xticks(np.arange(min(r), max(r)+1, 5.))
plt.savefig('pdf/thickwallFEMpo.pdf', bbox_inches = 'tight',
    pad_inches = 0.1, transparent=True)


## STRESSES THIN WALL
a = 95
rs = np.linspace(a,b,1000)

sigmatt = pi*a/(b-a)*np.ones(len(rs))
sigmart = -pi/2*np.ones(len(rs))

r = np.linspace(a,b,1000)
sigmatTWT = pi*a**2*(1+b**2/r**2)/(b**2-a**2)
sigmarTWT = pi*a**2*(1-b**2/r**2)/(b**2-a**2)

## FEM
dfT = pd.read_csv('dados/thinwall.csv')
xT = dfT[dfT.columns[16]]
sigmatTFEM = dfT[dfT.columns[4]]
sigmarTFEM = dfT[dfT.columns[3]]

plt.figure(3,figsize=(4,6))
plt.subplot(211)
plt.plot(rs,sigmatt,'k:',label='Thin wall theory')
plt.plot(r,sigmatTWT,'k-',label='Lamé')
plt.plot(xT,sigmatTFEM,'b-.',label='FEM')
plt.ylim((18,20))
plt.xlim((95,100))
plt.ylabel(r'$\sigma_t$ / MPa')
plt.grid()
plt.legend(loc='best')
plt.subplot(212)
plt.plot(rs,sigmart,'k:',label='SNCTTI')
plt.plot(r,sigmarTWT,'k-',label='Lamé')
plt.plot(xT,sigmarTFEM,'b-.',label='FEM')
plt.ylim((-2,1))
plt.xlim((95,100))
plt.legend(loc='best')
plt.ylabel(r'$\sigma_r$ / MPa')
plt.xlabel('r / mm')
plt.grid()
plt.savefig('pdf/thinwallFEM.pdf', bbox_inches = 'tight',
    pad_inches = 0.1, transparent=True)

## ERROR ANALYSIS
et = (max(sigmat)-max(sigmatFEM))/max(sigmatFEM)*100
etmin = (min(sigmat)-min(sigmatFEM))/min(sigmatFEM)*100
er = (max(sigmar)-max(sigmarFEM))/max(sigmarFEM)*100
ermin = (min(sigmar)-min(sigmarFEM))/min(sigmarFEM)*100

etpo = (max(sigmatpo)-max(sigmatFEMpo))/max(sigmatFEMpo)*100
erpo = (max(sigmarpo)-max(sigmarFEMpo))/max(sigmarFEMpo)*100

