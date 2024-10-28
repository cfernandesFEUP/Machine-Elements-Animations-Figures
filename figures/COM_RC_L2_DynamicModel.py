import numpy as np
import matplotlib.pyplot as plt
import scienceplots
plt.style.use(['science','ieee'])

Tl = 100

I = np.array([1,5,50])
F = ['-','-.','--']
j = 0
theta = np.linspace(0,2*np.pi,100)
plt.figure(1)
plt.plot(theta*180/np.pi,Tl*(1+np.sin(theta)),'k-')
plt.ylabel(r'$T_m$ / Nm')
#plt.grid()
plt.savefig('pdf/turningmoment.pdf', bbox_inches = 'tight',
    pad_inches = 0.1, transparent=True)
plt.figure(2)
for i in I:
    alpha = (Tl*np.sin(theta))/i
    plt.plot(theta*180/np.pi,alpha,'k'+F[j],label=str(i)+'I')
    plt.xlabel(r'$\theta$ / $^\circ$')
    plt.ylabel(r'$\ddot{\theta}$ / rad s$^{-2}$')
    plt.xlim([0,360])
    plt.legend()
    #plt.grid()
    j += 1
plt.savefig('pdf/flywheelDM.pdf', bbox_inches = 'tight', transparent=True)