import matplotlib.pyplot as plt
import numpy as np
from fractions import Fraction
## ADIMENSIONAL ###############################################################
def REYNOLDS(l,h0,sh,U,eta):
    l = l/1000
    h0 = h0/1000
    sh = sh/1000
    xvec = np.linspace(0,l,100)
    H0 = h0/sh
    X = xvec/l
    H = H0 + 1 -X
    Hm = 2*H0*(1 + H0)/(1 + 2*H0)
    P = 6*X*(1 - X)/((1+2*H0)*(H0 + 1 - X)**2)
    dPdX = 6*((H-Hm)/H**3)
    pdim = P*eta*U*l/(sh**2)
    cof = (2*sh*np.log(H0/(H0+1))+3*sh/(1+2*H0))/(3*l*np.log(H0/(H0+1))+6*l/(1+2*H0))
    return P, X, dPdX, H0, pdim, cof
H0 = np.linspace(0.11,2,1000)
Wy = 6*np.log((H0+1)/H0)-12/(1+2*H0)
Fb = 4*np.log(H0/(H0+1))+6/(1+2*H0)
Fa = 2*np.log(H0/(H0+1))+6/(1+2*H0)
sh = 0.1/1000
l = 0.1
cof = l*(2*sh*np.log(H0/(H0+1))+3*sh/(1+2*H0))/(3*l*np.log(H0/(H0+1))+6*l/(1+2*H0))/sh     
## 4 CONDIÇÕES ################################################################
P1, X1, dPdX1, H01, pdim1, cof1 = REYNOLDS(l=100,h0=0.1,sh=0.1,U=10,eta=0.04)
P2, X2, dPdX2, H02, pdim2, cof2 = REYNOLDS(l=100,h0=0.1,sh=0.2,U=10,eta=0.04)
P3, X3, dPdX3, H03, pdim3, cof3 = REYNOLDS(l=100,h0=0.1,sh=0.3,U=10,eta=0.04)
P4, X4, dPdX4, H04, pdim4, cof4 = REYNOLDS(l=100,h0=0.1,sh=0.4,U=10,eta=0.04)
## PLOT #######################################################################
# PRESSÃO
plt.figure(1)
plt.plot(X1,P1,'k',lw=1.5,label='H0='+str(Fraction(H01)))
plt.plot(X2,P2,'k--',lw=1.5,label='H0='+str(Fraction(H02).limit_denominator()))
plt.plot(X3,P3,'k-.',lw=1.5,label='H0='+str(Fraction(H03).limit_denominator()))
plt.plot(X4,P4,'k:',lw=1.5,label='H0='+str(Fraction(H04).limit_denominator()))
plt.legend(loc=2)
plt.xlabel(r'$X$')
plt.ylabel(r'$P$')
# CARGA
plt.figure(2)
plt.plot(H0, Wy,'k',lw=1.5)
plt.xlim([0,2])
plt.xlabel(r'$H0$')
plt.ylabel(r'$W_y$')
# CARGA COMPONENTES
plt.figure(3)
plt.plot(H0, Wy,'k',lw=1.5, label=r'$W_y$')
plt.plot(H0, Fa,'k--',lw=1.5, label=r'$F_a$')
plt.plot(H0, Fb,'k-.',lw=1.5, label=r'$F_b$')
plt.legend(loc=1)
plt.xlim([0,2])
plt.xlabel(r'$H0$')
plt.ylabel(r'$W_y, F_a, F_b$')
# ATRITO
plt.figure(4)
plt.plot(H0, cof,'k',lw=1.5)
plt.xlim([0,2])
plt.xlabel(r'$H0$')
plt.ylabel(r'$\mu l/s_h$')

# plt.figure(3)
# plt.plot(X1,dPdX1,'k',lw=1.5,label='H0='+str(Fraction(H01)))
# plt.plot(X2,dPdX2,'k--',lw=1.5,label='H0='+str(Fraction(H02).limit_denominator()))
# plt.plot(X3,dPdX3,'k-.',lw=1.5,label='H0='+str(Fraction(H03).limit_denominator()))
# plt.plot(X4,dPdX4,'k:',lw=1.5,label='H0='+str(Fraction(H04).limit_denominator()))
# plt.xlabel(r'$X$')
# plt.ylabel(r'$dP/dX$')

# plt.figure(2)
# plt.plot(X,W,'k',lw=1.5)
# plt.xlabel(r'$X$')
# plt.ylabel(r'$W$')

# plt.figure(3)
# plt.plot(xvec,p/1e6,'k',lw=1.5)
# plt.xlabel('x / mm')
# plt.ylabel('p / MPa')