import numpy as np
import matplotlib.pyplot as plt
import scienceplots

plt.style.use(['science','ieee'])

roh = 7850
nu = 0.3
g = 9.81
w = 10000*np.pi/30
b = 0.1

C = b**2*g*roh*w**2*(3 + nu)/(8*g)
D = (1 + 3*nu)/(3 + nu)
def sigmar(alpha,x,C):
    return C*(1+alpha**2-x**2-alpha**2/x**2)

def sigmat(alpha,x,C,D):
    return C*(1+alpha**2-D*x**2+alpha**2/x**2)

alpha0 = 0.00
alpha01 = 0.01
alpha02 = 0.05
alpha03 = 0.10
alpha1 = 0.25
alpha2 = 0.50
alpha3 = 0.8

def x(alpha):
    return np.linspace(alpha,1,100)

plt.figure(1)
plt.plot(x(alpha0),1e-6*sigmar(alpha0,x(alpha0),C),label='solid')
plt.plot(x(alpha01),1e-6*sigmar(alpha01,x(alpha01),C),label='K='+str("%.0f" % alpha01**-1))
plt.plot(x(alpha02),1e-6*sigmar(alpha02,x(alpha02),C),label='K='+str("%.0f" % alpha02**-1))
plt.plot(x(alpha03),1e-6*sigmar(alpha03,x(alpha03),C),label='K='+str("%.0f" % alpha03**-1))
plt.plot(x(alpha1),1e-6*sigmar(alpha1,x(alpha1),C),label='K='+str("%.0f" % alpha1**-1))
plt.plot(x(alpha2),1e-6*sigmar(alpha2,x(alpha2),C),label='K='+str("%.0f" % alpha2**-1))
plt.plot(x(alpha3),1e-6*sigmar(alpha3,x(alpha3),C),label='K='+str("%.2f" % alpha3**-1))
plt.legend(loc='upper right')
plt.xlim(0,1)
# plt.ylim(0,100)
#plt.grid()
plt.xlabel(r'$x=\frac{r}{b}$')
plt.ylabel(r'$\sigma_r$ / MPa')
plt.savefig('pdf/rotatingcylinderSr.pdf', bbox_inches = 'tight',
    pad_inches = 0.1, transparent=True)


plt.figure(2)
plt.plot(x(alpha0),1e-6*sigmat(alpha0,x(alpha0),C,D),label='solid')
plt.plot(x(alpha01),1e-6*sigmat(alpha01,x(alpha01),C,D),label='K='+str("%.0f" % alpha01**-1))
plt.plot(x(alpha02),1e-6*sigmat(alpha02,x(alpha02),C,D),label='K='+str("%.0f" % alpha02**-1))
plt.plot(x(alpha03),1e-6*sigmat(alpha03,x(alpha03),C,D),label='K='+str("%.0f" % alpha03**-1))
plt.plot(x(alpha1),1e-6*sigmat(alpha1,x(alpha1),C,D),label='K='+str("%.0f" % alpha1**-1))
plt.plot(x(alpha2),1e-6*sigmat(alpha2,x(alpha2),C,D),label='K='+str("%.0f" % alpha2**-1))
plt.plot(x(alpha3),1e-6*sigmat(alpha3,x(alpha3),C,D),label='K='+str("%.2f" % alpha3**-1))
# plt.legend(loc='best')
plt.xlim(0,1)
# plt.ylim(0,200)
#plt.grid()
plt.xlabel(r'$x=\frac{r}{b}$')
plt.ylabel(r'$\sigma_t$ / MPa')
plt.savefig('pdf/rotatingcylinderSt.pdf', bbox_inches = 'tight',
    pad_inches = 0.1, transparent=True)