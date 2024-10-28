import matplotlib.pyplot as plt
import numpy as np

z = np.linspace(10,300,1000)

y = 0.154 - 0.912/z
y14 = 0.124 - 0.684/z

YL = 1/(np.pi*y)
YL14 = 1/(np.pi*y14)

plt.figure(1)
plt.plot(z,YL,'k',label=r'Wallace $\alpha$=20$^\circ$')
plt.plot(z,YL14,'k--',label=r'Wallace $\alpha$=14.5$^\circ$')
plt.xlabel(r'$z$')
plt.ylabel(r'$Y_L$')
plt.grid()
plt.legend()
plt.savefig('pdf/LewisForm.pdf', bbox_inches = 'tight',
    pad_inches = 0.1, transparent=True)