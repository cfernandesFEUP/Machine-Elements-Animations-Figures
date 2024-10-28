import numpy as np, matplotlib.pyplot as plt

n = np.linspace(2,6, 100)


def SN(a,b,n):
    S = a*10**(b*n)
    return S

aF = 1.62
bF = -0.085091


aT = 1.88372
bT = -0.106935

aC = 1.93966
bC = -0.137554

import matplotlib as mpl
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1) 
plt.loglog(10**n, SN(aF, bF, n),'k-',label='Flexão')
plt.loglog([10**6,10**7], [SN(aF, bF, 6), SN(aF, bF, 6)],'k-')
plt.loglog(10**n, SN(aT, bT, n),'k-.',label='Tração')
plt.loglog([10**6,10**7], [SN(aT, bT, 6), SN(aT, bT, 6)],'k-.')
plt.loglog(10**n, SN(aC, bC, n),'k--',label='Corte')
plt.loglog([10**6, 10**7], [SN(aC, bC, 6), SN(aC, bC, 6)],'k--')
plt.xlim(1e3,1e7)
plt.ylim(0.2,1.)
plt.ylabel(r'Tensão Limite Fadiga / Tensão de Rotura')
plt.xlabel(r'Ciclos')
plt.grid(b=True, which='major', color='gray', linestyle='-')
plt.grid(b=True, which='minor', color='gray', linestyle='-')
ax.yaxis.set_minor_formatter(mpl.ticker.StrMethodFormatter('{x:.1f}'))
ax.yaxis.set_major_formatter(mpl.ticker.StrMethodFormatter('{x:.1f}'))
plt.legend()
plt.savefig('pdf/FAD_T3_SN.pdf', bbox_inches = 'tight',
    pad_inches = 0, transparent=True)

