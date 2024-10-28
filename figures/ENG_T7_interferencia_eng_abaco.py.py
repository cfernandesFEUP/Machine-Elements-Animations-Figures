import numpy as np, matplotlib.pyplot as plt

z2 = np.linspace(6,1000, 10000)

alphaG = np.array([15., 20., 25., 30.])

alpha = alphaG*np.pi/180

form=['k--','k-','k-.','k:']

def func(alpha,z2):
    z1 = -z2 + np.sqrt(z2**2 + 4*(z2+1)/(np.sin(alpha)**2))
    z1 = z1[z2>z1]
    return z1




import matplotlib as mpl
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1) 
for i in range(len(alpha)):
    plt.semilogx(z2[z2>func(alpha[i],z2)[0]],func(alpha[i],z2), form[i], label=r'$\alpha=$'+str(int(alphaG[i]))+'\u00b0')
plt.semilogx(z2,z2,'k',lw=0.5)
plt.xlim(6,100)
plt.ylim(6,30)
plt.xlabel(r'$z_2$')
plt.ylabel(r'$z_1$')
plt.grid(b=True, which='major', color='gray', linestyle='-')
plt.grid(b=True, which='minor', color='gray', linestyle='-')
ax.xaxis.set_minor_formatter(mpl.ticker.StrMethodFormatter('{x:.0f}'))
ax.xaxis.set_major_formatter(mpl.ticker.StrMethodFormatter('{x:.0f}'))
plt.legend()

