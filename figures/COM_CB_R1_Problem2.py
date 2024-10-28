import numpy as np
import matplotlib.pyplot as plt

## PROBLEM 1

phi = np.linspace(0,np.pi/2,100)
M_Pr = np.cos(phi)-1

plt.figure(1)
plt.plot(phi*180/np.pi,M_Pr)
plt.xlabel(r'$\phi$ / $^\circ$')
plt.ylabel(r'$\frac{M}{Pr}$')


## PROBLEM 2

plt.figure(2)
M_Pr0 = 0 + np.sin(phi)
M_Pr1 = 0.5 + np.sin(phi)
M_Pr2 = 1 + np.sin(phi)
M_Pr3 = 1.5 + np.sin(phi)
M_Pr4 = 2 + np.sin(phi)

plt.plot(phi*180/np.pi,M_Pr0,label='l/r = 0')
plt.plot(phi*180/np.pi,M_Pr1)
plt.plot(phi*180/np.pi,M_Pr2)
plt.plot(phi*180/np.pi,M_Pr3)
plt.plot(phi*180/np.pi,M_Pr4)
plt.xlabel(r'$\phi$ / $^\circ$')
plt.ylabel(r'$\frac{M}{Pr}$')