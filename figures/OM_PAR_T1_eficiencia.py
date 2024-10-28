import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

cof = np.array([0.01, 0.02, 0.05, 0.1, 0.15, 0.2])
gama = np.linspace(0,np.pi/2,2000)
beta = 0*np.pi/180

mpl.rcParams["font.size"] = 12
plt.figure(1)
plt.xlim([0,90])
plt.ylim([0,100])
plt.xlabel(r'$\gamma$ / $^{\circ}$')
plt.ylabel(r'$\eta$ / $\%$')
# plt.axis('equal')
for i in cof:
    # eff = 100*(np.cos(beta)-i*np.tan(gama))/(np.cos(beta)+i*np.cos(gama))
    eff = 100*np.tan(gama)*(1-i*np.tan(gama)/np.cos(beta))/(i/np.cos(beta)+np.tan(gama))
    plt.plot(gama[eff>0]*180/np.pi,eff[eff>0],'k',label=str(i))

mpl.rcParams["font.size"] = 8
from labellines import labelLines 
labelLines(plt.gca().get_lines(),zorder=2.5)
plt.savefig('pdf/eficiencia0.pdf', bbox_inches = 'tight',
    pad_inches = 0.1, transparent=True)
plt.show()


gama = 40*np.pi/180
beta = np.linspace(0,60*np.pi/180,2000)

mpl.rcParams["font.size"] = 12
plt.figure(2)
plt.xlim([0,60])
plt.ylim([0,100])
plt.xlabel(r'$\alpha_n$ / $^{\circ}$')
plt.ylabel(r'$\eta$ / $\%$')
# plt.axis('equal')
for i in cof:
    # eff = 100*(np.cos(beta)-i*np.tan(gama))/(np.cos(beta)+i*np.cos(gama))
    eff = 100*np.tan(gama)*(1-i*np.tan(gama)/np.cos(beta))/(i/np.cos(beta)+np.tan(gama))
    plt.plot(beta[eff>0]*180/np.pi,eff[eff>0],'k',label=str(i))
mpl.rcParams["font.size"] = 8
from labellines import labelLines 
labelLines(plt.gca().get_lines(),zorder=2.5)
plt.savefig('pdf/eficienciagama40.pdf', bbox_inches = 'tight',
    pad_inches = 0.1, transparent=True)
plt.show()
