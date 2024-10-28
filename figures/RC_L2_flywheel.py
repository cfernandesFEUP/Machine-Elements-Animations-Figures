import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate

CS = 0.1
omega = 250

thetai = np.linspace(0,705,48)

## TURNING MOMENT
Ti = np.array([0,2800,2090,2430,2160,1840,1590,1210,1066,803,532,184,0,-107,\
               -206,-280,-323,-310,-242,-126,-8,89,125,85,0,-85,-125,-89,8,\
                126,242,310,323,280,206,107,0,-107,-206,-292,-355,-371,-362,\
                    -312,-272,-274,-548,-760])

TiSI = 0.112984829*Ti
    
U = integrate.trapz(Ti, thetai*np.pi/180)
Tm = U/max(thetai*np.pi/180)
deltaT = Ti - Tm
deltaU = integrate.trapz(deltaT[0:12], thetai[0:12]*np.pi/180)
I = deltaU/(CS*omega**2)

USI = integrate.trapz(TiSI, thetai*np.pi/180)
TmSI = USI/max(thetai*np.pi/180)
deltaTSI = TiSI - TmSI
deltaUSI = integrate.trapz(deltaTSI[0:12], thetai[0:12]*np.pi/180)
ISI = deltaUSI/(CS*omega**2)

print(USI, TmSI, deltaUSI, ISI)

plt.figure(1)
plt.plot(thetai, TiSI,'k-o')
plt.grid()
plt.xlabel(r'Crank angle, $\theta$ / $^\circ$')
plt.ylabel('Crank torque, T / Nm')
plt.savefig('pdf/crankEX.pdf', bbox_inches = 'tight',
    pad_inches = 0.1, transparent=True)
