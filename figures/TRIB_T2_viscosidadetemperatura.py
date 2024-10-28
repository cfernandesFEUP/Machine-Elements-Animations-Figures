import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
## MEDIÇÕES VISCOSIDADE #######################################################
Texp = np.array([40,70,100])
viscexp = np.array([319.22,65.81,22.33])
xdata = Texp + 273.15
ydata = viscexp
Tlub = np.linspace(20,100,100)
TKelvin = Tlub + 273.15
roh = 902
alphaT = -5.8e-4
## LEI ASTM ###################################################################
def ASTM(T,c,m,n):
    return -c+10**(10**(m - n*np.log10(T)))
p01 = [0.7, 9., 3.5]
popt1, pcov1 = curve_fit(ASTM, xdata, ydata,p01)
niuASTM = ASTM(TKelvin,popt1[0],popt1[1],popt1[2])
## VOGEL ######################################################################
def VOGEL(T,a,b,c): 
    return a*np.exp(b/(T-c))
p02 = [0.1 , 1000., -100.]
popt2, pcov2 = curve_fit(VOGEL, xdata, ydata,p02)
niuVOGEL = VOGEL(TKelvin,popt2[0],popt2[1],popt2[2])
## PLOT #######################################################################
plt.figure(1)
plt.plot(Texp, viscexp,'ks', ms=6,label='experimental')
plt.plot(Tlub,niuASTM,'k',lw=1.5,label='ASTM')
plt.grid()
plt.legend(loc=1)
plt.xlabel(r'T / $^\circ$C')
plt.ylabel(r'$\nu~/~\mathrm{mm^2/s}$')
plt.savefig('pdf/viscosidade.pdf', bbox_inches = 'tight', transparent=True)

plt.figure(2)
plt.plot(Texp, viscexp,'ks', ms=6,label='experimental')
plt.plot(Tlub,niuASTM,'k',lw=1.5,label='ASTM')
plt.plot(Tlub,niuVOGEL,'k--',lw=1.5,label='VOGEL')
plt.legend(loc=1)
plt.grid()
plt.xlabel(r'T / $^\circ$C')
plt.ylabel(r'$\nu~/~\mathrm{mm^2/s}$')
plt.savefig('pdf/viscosidadelog.pdf', bbox_inches = 'tight',transparent=True)