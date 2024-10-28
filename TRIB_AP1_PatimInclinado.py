import numpy as np

import matplotlib.pyplot as plt

## INPUT ######################################################################

L = float(input('Dimensão do patim segundo z, L / mm: '))/1000

B = float(input('Dimensão do patim segundo x, B / mm: '))/1000

h1 = float(input(r'Altura $h_1$ / mm: ' ))/1000

h2 = float(input(r'Altura $h_2$ / mm: '))/1000

U = float(input('Velocidade da placa móvel U / m/s: '))

eta = float(input(r'Viscosidade do lubrificante $\eta$ / mPa$\cdot$s: '))/1000

sh = h1-h2

xvec = np.linspace(0,B,100)

hvec = h1 + xvec*(h2-h1)/B

hm = 2*h1*h2/(h1+h2)

## CAMPO DE PRESSÃO ###########################################################
ph = 6*U*B*eta/(h2-h1)*((-1/hvec)+(1/hvec**2)*(h1*h2/(h1+h2))+1/(h1+h2))

dpdh = 6*U*B*eta/(h2-h1)*(hvec-hm)/(hvec**3)

W = 6*eta*U*L*(B/(h2-h1))**2*(np.log(h1/h2)-2*(h1-h2/(h1+h2)))

Ft = eta*U*L*B/(h2-h1)*(4*np.log(h2/h1)+6*(h1-h2)/(h1+h2))

cof = Ft/W

tyx = eta*U*(3*hm/hvec**2-4/hvec) 

Q = L*U*h1*h2/(h1+h2)

print(u'$h_m$=', "%.4f" % hm)

plt.figure(1)
plt.title('Campo de pressão')
plt.plot(xvec*1000, ph/1e6)
plt.xlabel('x / mm')
plt.ylabel('p / MPa')
plt.show()

plt.figure(2)
plt.plot(xvec*1000, dpdh/1e6)
plt.xlabel('x / mm')
plt.ylabel('dp/dx / MPa/m')
plt.show()

plt.figure(3)
plt.plot(xvec*1000, tyx)
plt.xlabel('x / mm')
plt.ylabel(u'$\tau_{yx}$ / MPa')
plt.show()

input('Prima uma tecla para continuar')