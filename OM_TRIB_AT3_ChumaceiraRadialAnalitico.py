import matplotlib.pyplot as plt
import numpy as np
import scienceplots
plt.style.use('science')
## DADOS CHUMACEIRA ###########################################################
epsilon = 0.1
phi = np.linspace(0,2*np.pi,100)
phiG = np.linspace(0,np.pi,50)
phiG1 = np.linspace(np.pi,2*np.pi,50)
## SOMMERFELD #################################################################
PS = 6*epsilon*np.sin(phi)*(2+np.cos(phi))/((2+epsilon**2)*(1+epsilon*np.cos(phi))**2)
## HALF SOMMERFELD - GUMBEL ###################################################
PG0 = 6*epsilon*np.sin(phiG)*(2+np.cos(phiG))/((2+epsilon**2)*(1+epsilon*np.cos(phiG))**2)
PG1 = np.zeros(50)
## PLOT #######################################################################
plt.figure(1)
plt.plot(phi, PS,'k')
plt.xlabel(r'$\theta~/~rad$')
plt.ylabel(r'$P$')
pi = np.pi
plt.xticks(np.arange(0, pi+pi+pi/2, step=(pi/2)), ['0',r'$\pi$/2',r'$\pi$',r'3$\pi$/2',r'2$\pi$'])
plt.savefig('figures/pdf/sommerfeld.pdf', bbox_inches = 'tight',
    pad_inches = 0, transparent=True)

plt.figure(2)
plt.plot(phiG, PG0,'k')
plt.plot(phiG1, PG1,'k')
#plt.grid()
plt.xlabel(r'$\theta~/~rad$')
plt.ylabel(r'$P$')
pi = np.pi
plt.xticks(np.arange(0, pi+pi+pi/2, step=(pi/2)), ['0',r'$\pi$/2',r'$\pi$',r'3$\pi$/2',r'2$\pi$'])
plt.savefig('figures/pdf/gumbel.pdf', bbox_inches = 'tight',
    pad_inches = 0, transparent=True)
