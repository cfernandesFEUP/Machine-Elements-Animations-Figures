import numpy as np
import matplotlib.pyplot as plt

n = 100
beta  = np.linspace(0.001,0.9999,n)


## FRICTION RADIUS
ri = 1
ro = ri/beta

Rp = (2/3*(ro**3-ri**3)/(ro**2-ri**2))/ro
Rw = ((ro+ri)/2)/ro


## DIMENSIONLESS TORQUE
Tw = (1+beta)/4
Tp = (1-beta**3)/(3*(1-beta**2))


plt.figure(1,figsize = (5, 4))
plt.plot(beta,Tw,'k-',label='uniform wear')
plt.plot(beta,Tp,'k-.',label='uniform pressure')
plt.grid()
plt.legend()
plt.xlabel(r'Radius ratio, $\beta$')
plt.ylabel(r'Dimensionless torque, $\overline{T}$')
plt.savefig('pdf/clutchTorque.pdf', bbox_inches = 'tight',
    pad_inches = 0, transparent=True)
plt.show()


plt.figure(2,figsize = (5, 4))
plt.plot(beta,Rw,'k-',label='uniform wear')
plt.plot(beta,Rp,'k-.',label='uniform pressure')
plt.grid()
plt.legend()
plt.xlabel(r'Radius ratio, $\beta$')
plt.ylabel('Friction radius')
plt.savefig('pdf/clutchFrictionRadius.pdf', bbox_inches = 'tight',
    pad_inches = 0, transparent=True)
plt.show()