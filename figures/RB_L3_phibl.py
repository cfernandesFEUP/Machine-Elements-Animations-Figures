import numpy as np
import matplotlib.pyplot as plt

## FACT = (n*v)**1.4*dm


FACT = np.linspace(1e5,1e9,2000)

phibl = 1/np.exp(2.6e-8*FACT)


## FIGURE
plt.figure(1)
plt.semilogx(FACT,phibl,'k')
plt.xlabel(r'$\left(n\cdot \nu\right)^{1.4}\cdot d_m$')
plt.ylabel(r'$\phi_{bl}$')
plt.grid()
plt.savefig('pdf/phi_bl.pdf', bbox_inches = 'tight',
    pad_inches = 0.1, transparent=True)