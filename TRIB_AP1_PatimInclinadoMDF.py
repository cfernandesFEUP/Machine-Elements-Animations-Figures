import numpy as np

## INPUT ######################################################################

X1 = float(input('Dimensão segundo X / mm: '))

Y1 = float(input('Dimensão segundo Y / mm: '))

# H1 = float(input('Altura h1 / mm: '))

# H2 = float(input('Altura h2 / mm: '))

# DECLIVE = (H2-H1)/X1

N = 100 #int(input('Número de nós segundo X: '))

M = 100 #int(input('Número de nós segundo Z: '))

C = 0.030

## MDF ########################################################################

deltaX = 1/N

deltaY = 1/M

C1 = X1*X1/(Y1*Y1)

ITER = 1000

p = np.zeros((N+2,M+2))

X = np.zeros(N)

Y = np.zeros(M)

SUM = np.zeros(ITER+1)

error = 1

for K in range(ITER): 
    SUMIJ = 0 
    for I in range(1,N+1):
        X[I-1] = 1/N*(I-1)
        h = 2/3*(2-X[I-1])
        hm = 2/3*(2-X[I-1]-0.5*deltaX)
        hp = 2/3*(2-X[I-1]+0.5*deltaX)
        hm1 = 2/3*(2-X[I-1]-deltaX)
        hp1 = 2/3*(2-X[I-1]+deltaX)
        cubh = h*h*h
        cubhm = hm*hm*hm
        cubhp = hp*hp*hp
        C2 = (cubhp + cubhm + 2*C1*cubh);
        A = C1*cubh/C2
        CA = cubhp/C2
        D = cubhm/C2
        E = -0.5*deltaX*(hp1 - hm1)*C/(C2*X1) 
        for J in range(1,M+1):
            Y[J-1] = 1/M*(J-1)
            p[I,J] = A*p[I,J+1] + A*p[I,J-1] + CA*p[I+1,J] + D*p[I-1,J] - E
            SUMIJ = SUMIJ + p[I,J]
    SUM[K+1] = SUMIJ
    while error < 0.1:
        error = abs(SUM[K+1]-SUM[K])/abs(SUM[K+1])
print('Número de iterações', K)
## PLOT #######################################################################

import matplotlib.pyplot as plt

from matplotlib import cm

x = np.linspace(0, X1, N+2)
y = np.linspace(0, Y1, M+2)
X, Y = np.meshgrid(x, y)

fig = plt.figure(figsize=(8,5))
ax = plt.axes(projection='3d')
ax.plot_surface(X, Y, p, cmap=cm.coolwarm, linewidth=0.5, color='black')

ax.set_xlabel('B / mm', fontsize=12)
ax.set_ylabel('L / mm', fontsize=12)
ax.set_zlabel('p', fontsize=12)
ax.set_xlim([0,X1])
ax.set_ylim([0,Y1])
ax.set_zlim([0,1.1*p.max()])
plt.show()