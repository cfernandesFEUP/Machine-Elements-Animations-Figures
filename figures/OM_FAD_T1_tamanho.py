# import tikzplotlib as tik
import numpy as np
import matplotlib.pyplot as plt
n = 150
d1 = np.linspace(8,51,n)
d2 = np.linspace(51,254,n)
d0 = np.linspace(8,254,n)
Cb0 = np.linspace(1,1,n)
Cb1 = 1.24*d1**(-0.107)
Cb2 = 1.51*d2**(-0.157)
plt.figure(1)
plt.plot(d0,Cb0,'k',linewidth=5, label = 'Axial')
plt.plot(d1,Cb1,'k-.',linewidth=3, label = 'Flexão')
plt.plot(d2,Cb2,'k-.',linewidth=3)
plt.xlabel('diâmetro / mm')
plt.ylabel(r'$\mathtt{C_2}$')
plt.xlim([0,250])
plt.ylim([0,1.])
plt.grid(True)
plt.legend(loc=3, prop={'size': 12})
plt.savefig('pdf/FAD_T1_tamanho.pdf', bbox_inches = 'tight',
    pad_inches = 0, transparent=True)
# tik.save('tamanho.tex')
