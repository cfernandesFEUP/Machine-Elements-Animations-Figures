import numpy as np
import matplotlib.pyplot as plt
import scienceplots

plt.style.use(['science','ieee'])

roh = 7850
nu = 0.3
g = 9.81
w = np.pi/30
b = 0.1

C = g*b**2*roh*w**2*(3 + nu)/(8*g)
D = (1 + 3*nu)/(3 + nu)
def sigmar(alpha,x,C,n):
    return n**2*C*(1+alpha**2-x**2-alpha**2/x**2)

def sigmat(alpha,x,C,D,n):
    return n**2*C*(1+alpha**2-D*x**2+alpha**2/x**2)

alpha = 0.5

x = np.linspace(alpha,1,100)

n = [1000,2000,5000,10000]

plt.figure(1)
plt.plot(x,1e-6*sigmar(alpha,x,C,n[0]),label='n='+str(n[0])+' rpm')
plt.plot(x,1e-6*sigmar(alpha,x,C,n[1]),label='n='+str(n[1])+' rpm')
plt.plot(x,1e-6*sigmar(alpha,x,C,n[2]),label='n='+str(n[2])+' rpm')
plt.plot(x,1e-6*sigmar(alpha,x,C,n[3]),label='n='+str(n[3])+' rpm')
plt.legend(loc='best')
plt.xlim(0.5,1)
plt.ylim(0,10)
#plt.grid()
plt.xlabel(r'$x=\frac{r}{b}$')
plt.ylabel(r'$\sigma_r$ / MPa')
plt.savefig('pdf/rotatingcylinderSrspeed.pdf', bbox_inches = 'tight',
    pad_inches = 0.1, transparent=True)

plt.figure(2)
plt.plot(x,1e-6*sigmat(alpha,x,C,D,n[0]),label='n='+str(n[0])+' rpm')
plt.plot(x,1e-6*sigmat(alpha,x,C,D,n[1]),label='n='+str(n[1])+' rpm')
plt.plot(x,1e-6*sigmat(alpha,x,C,D,n[2]),label='n='+str(n[2])+' rpm')
plt.plot(x,1e-6*sigmat(alpha,x,C,D,n[3]),label='n='+str(n[3])+' rpm')
plt.legend(loc='best')
plt.xlim(0.5,1)
plt.ylim(0,100)
#plt.grid()
plt.xlabel(r'$x=\frac{r}{b}$')
plt.ylabel(r'$\sigma_t$ / MPa')
plt.savefig('pdf/rotatingcylinderStspeed.pdf', bbox_inches = 'tight',
    pad_inches = 0.1, transparent=True)