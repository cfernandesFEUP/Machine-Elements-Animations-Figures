import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

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

n = 100

alpha = np.linspace(0.01,1,n)

x = np.zeros([len(alpha),n])
Sr = np.zeros([len(alpha),n])
St = np.zeros([len(alpha),n])
for i in range(len(alpha)):    
    x[i] = np.linspace(alpha[i],1,n)
    Sr[i] = 1e-6*sigmar(alpha[i],x[i],C)
    St[i] = 1e-6*sigmat(alpha[i],x[i],C,D)

Y= x*b*1000
X = alpha
# X = np.meshgrid(alpha,Y)

fig = plt.figure(1)
ax = plt.axes(projection='3d')
ax.plot_surface(X, Y, Sr,cmap='jet', edgecolor='k')
ax.set_xlabel(r'$\alpha$')
ax.set_ylabel(r'r / mm')
ax.set_zlabel(r'$\sigma_r$ / MPa')
# ax.view_init(10, 15)
plt.show()

# Normalize to [0,1]
norm = plt.Normalize(Sr.min(), Sr.max())
colors = cm.jet(norm(Sr))
rcount, ccount, _ = colors.shape
fig = plt.figure()
ax = fig.gca(projection='3d')
surf = ax.plot_surface(X, Y, Sr, rcount=rcount, ccount=ccount,
                       facecolors=colors, shade=False)
ax.set_xlabel(r'$\alpha$')
ax.set_ylabel(r'r / mm')
ax.set_zlabel(r'$\sigma_r$ / MPa')
surf.set_facecolor((0,0,0,0))
plt.show()

fig = plt.figure(2)
ax = plt.axes(projection='3d')
surf = ax.plot_surface(X, Y, 1e-6*St,cmap='jet', edgecolor='k')
ax.set_xlabel(r'$\alpha$')
ax.set_ylabel(r'r / mm')
ax.set_zlabel(r'$\sigma_t$ / MPa')
# ax.view_init(10, 15)
plt.show()

# Normalize to [0,1]
norm = plt.Normalize(St.min(), St.max())
colors = cm.jet(norm(St))
rcount, ccount, _ = colors.shape
fig = plt.figure()
ax = fig.gca(projection='3d')
surf = ax.plot_surface(X, Y, St, rcount=rcount, ccount=ccount,
                       facecolors=colors, shade=False)
ax.set_xlabel(r'$\alpha$')
ax.set_ylabel(r'r / mm')
ax.set_zlabel(r'$\sigma_t$ / MPa')
surf.set_facecolor((0,0,0,0))
plt.show()

plt.plot(x,Sr)

# plt.plot(x,St)