import numpy as np, matplotlib.pyplot as plt

# geometry
a = 0.1
c = 1.
b = np.linspace(a,c,100)

# internal pressure
pi = 1.
# Equivalent stress
Se = 2*pi*c**2*(1-1/(b**2/(b**2-a**2) + c**2/(c**2-b**2)))/(c**2-a**2)

bmin = (a*c)**0.5
Semin =  2*pi*c**2*(1-1/(bmin**2/(bmin**2-a**2) + c**2/(c**2-bmin**2)))/(c**2-a**2)

plt.figure(1)
plt.plot(b,Se)
plt.plot(bmin,Semin, 'ks')
plt.xlabel('b')
plt.ylabel(r'$\sigma_e$')