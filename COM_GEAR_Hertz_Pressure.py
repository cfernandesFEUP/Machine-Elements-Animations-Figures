import numpy as np
import matplotlib.pyplot as plt


z1 = 20
z2 = 20
m = 2
b = 20
Es = 210e3 # MPa
nu = 0.3
Fn = 5000 #N
E = (2*(1-nu)/Es)**-1
r1 = z1*m/2
r2 = z2*m/2
alpha = np.pi/9
rb1 = r1*np.cos(alpha)
rb2 = r2*np.cos(alpha)
ra1 = r1 + m
ra2 = r2 + m
pb = np.pi*m*np.cos(alpha)
T1T2 = (r1+r2)*np.sin(alpha)
T2A = np.sqrt(ra2**2 - rb2**2)
T1A = T1T2 - T2A
AB = np.sqrt(ra1**2 - rb1**2) + np.sqrt(ra2**2 - rb2**2) - T1T2
x = np.linspace(0, AB, 200)
R1 = T1A + x
R2 = T2A - x
R = (1/R1 + 1/R2)**-1 

T1B = np.sqrt(ra1**2- rb1**2)
TR1 = T1B - pb - T1A

xi = x/pb
pH = np.zeros(len(x))
for i in range(len(x)):
    if x[i]<TR1:
        pH[i] = np.sqrt(0.5*Fn*E/(np.pi*b*R[i]))
    elif x[i]>pb:
        pH[i] = np.sqrt(0.5*Fn*E/(np.pi*b*R[i]))
    else:
        pH[i] = np.sqrt(Fn*E/(np.pi*b*R[i]))

plt.figure()
plt.plot(xi, pH/pH.max(), 'k')
plt.xlabel('$\overline{x}$')
plt.ylabel('$\overline{p_0}$')
plt.ylim(0.5,1.2)
plt.grid()
plt.savefig('figures/pdf/Hertz_Pressure_Equal_Tooth.pdf')

