import math
import numpy as np
import matplotlib.pyplot as plt

## 6206 DEEP GROOVE BALL BEARING
D = 62      
d = 30
dm = (D+d)/2
C0 = 11300

## MODEL COEFFICIENTS
f0 = 1.75    ## 1.5-2
f1 = 1.45

def Harris(f0,f1,dm,nu,n,F):
    ## COF DEEP GROOVE BALL BEARING
    cof = 0.002*(F/C0)**(1/2)
    ## NO-LOAD
    M0 = f0*1e-7*(nu*n)**(2/3)*dm**3
    ## LOAD-DEPENDENT
    M1 = cof*f1*F*dm/2
    ## TOTAL
    M = M0 + M1
    return M, M0, M1

## INFLUENCE LOAD
Fr = 2000
Fa = 1000


## OPERATING CONDITIONS
nu = 20                     ## KINEMATIC VISCOSITY / mm2/s
n = 3000                    ## SPEED / rpm
F = np.sqrt(Fr**2+Fa**2)    ## LOAD F = (Fr^2 + Fa^2)^0.5
vec = np.linspace(0,5000,1000)    

MF, M0F, M1F = Harris(f0,f1,dm,nu,n,vec)
MS, M0S, M1S = Harris(f0,f1,dm,nu,vec,F)

M,M0,M1 = Harris(f0,f1,dm,nu,n,F)

## FIGURE
plt.figure(1)
plt.plot(vec,MF,'k',label=r'$M$')
plt.plot([0,max(vec)],[M0F, M0F],'k--',label=r'$M_0$')
plt.plot(vec,M1F,'k-.',label=r'$M_1$')
plt.xlabel(r'Load / N ')
plt.ylabel(r'$M$ / Nmm')
plt.ylim(0,int(10*math.ceil(max(MF)/10)))
plt.grid()
plt.legend()
plt.savefig('pdf/HarrisLoad.pdf', bbox_inches = 'tight',
    pad_inches = 0.1, transparent=True)

# plt.figure(2)
# plt.plot([0,max(vec)],[M0F, M0F],'k--',label=r'$M_0$')
# plt.xlabel(r'Load / N ')
# plt.ylabel(r'$M_0$ / Nmm')
# plt.ylim(0,100)
# plt.grid()
# plt.savefig('pdf/HarrisLoadM0.pdf', bbox_inches = 'tight',
#     pad_inches = 0.1, transparent=True)

# plt.figure(3)
# plt.plot(vec,M1F,'k-.',label=r'$M_1$')
# plt.xlabel(r'Load / N ')
# plt.ylabel(r'$M_1$ / Nmm')
# plt.ylim(0,100)
# plt.grid()
# plt.legend()
# plt.savefig('pdf/HarrisLoadM1.pdf', bbox_inches = 'tight',
#     pad_inches = 0.1, transparent=True)

plt.figure(4)
plt.plot(vec,MS,'k',label=r'$M$')
plt.plot(vec,M0S,'k--',label=r'$M_0$')
plt.plot([0,max(vec)],[M1S, M1S],'k-.',label=r'$M_1$')
plt.xlabel(r'n / rpm ')
plt.ylabel(r'$M$ / Nmm')
plt.ylim(0,int(10*math.ceil(max(MS)/10)))
plt.grid()
plt.legend()
plt.savefig('pdf/HarrisSpeed.pdf', bbox_inches = 'tight',
    pad_inches = 0.1, transparent=True)