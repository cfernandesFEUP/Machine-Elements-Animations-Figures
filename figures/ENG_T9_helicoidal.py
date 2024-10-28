## IMPORTAR LIBRARIAS ##
import matplotlib.pyplot as plt
import numpy as np
import functions.dimensioning as dm
from matplotlib import rcParams
# rcParams['font.family'] = 'sans-serif'
# rcParams['font.sans-serif'] = ['Fira Sans']
rcParams['font.size'] = 15
rcParams['legend.fontsize'] = 'medium'
params = {'mathtext.default': 'regular' }  
## DADOS ENGRENAGENS ##
m = 2.
z1 = 20.
z2 = 10.
alpha = np.pi/6
beta = np.pi/10
alphat = np.arctan(np.tan(alpha)/np.cos(beta))
betab = np.arctan(np.tan(beta)*np.cos(alphat))
x1 = 0.
x2 = 0.
hfP = 1.
haP = 1.25
roh = 0.38
b = 100
## CALCULAR GEOMETRIA ##
r1, r2 = m*z1/(2*np.cos(beta)), m*z2/(2*np.cos(beta))
rl1, rl2 = r1 + m*x1, r2 + m*x2
ra1, ra2 = r1 + m*(1+x1), r2 + m*(1+x2)
rb1, rb2 = r1*np.cos(alphat), r2*np.cos(alphat)
a = rl1+rl2
p = np.pi*m
pb = p*np.cos(alpha)
s = p/2
## POSIÇÃO DOS CENTROS ##
O1 = np.array([0,a])
O2 = np.array([0,0])
## CÁLCULO DAS EVOLVENTES ##
evsize = 100
rasize = 10
T1T2 = a*np.sin(alpha)
T2B = np.sqrt(ra2**2-rb2**2)
T1A = np.sqrt(ra1**2-rb1**2)
T2A = T1T2 - T1A
T1B = T1T2 - T2B
rA2 = np.sqrt(rb2**2+T2A**2)
rAN2 = np.sqrt(rb2**2+(T2A+pb)**2)
## DEFINIÇÃO DOS ÂNGULOS DE ROTAÇÃO ##
alphaA1 = np.arctan(T1A/rb1)
alphaA2 = np.arctan(T2B/rb2) 
deltaA1 = -0.1
deltaA2 = -0.015
def phiA(ang_inc):
    return np.linspace(0, ang_inc + dm.inv(ang_inc), evsize)
## COORDENADAS DAS EVOLVENTES ##
xye1 = dm.ev(rb1,phiA(-alphaA1 + deltaA1))
xye2 = dm.ev(rb2,phiA(-alphaA2 + deltaA2))
xyi1 = dm.ev(rb1,phiA(-alphaA1))
xyi2 = dm.ev(rb2,phiA(-alphaA2))
## COORDENADAS DA CABEÇA ##
sa1 = ra1*(s/rl1+2*(dm.inv(alpha)-dm.inv(alphaA1)))
sa2 = ra2*(s/rl2+2*(dm.inv(alpha)-dm.inv(alphaA2)))
te1, te2 = sa1/ra1, sa2/ra2
teta1 = np.linspace(dm.inv(-alphaA1),dm.inv(-alphaA1)-te1,rasize)
teta2 = np.linspace(dm.inv(-alphaA2),dm.inv(-alphaA2)-te2,rasize)
xya1 = dm.dxy(ra1, 0, 0, teta1)
xya2 = dm.dxy(ra2, 0, 0, teta2)
vi1, vi2 = np.zeros((2,evsize+rasize)), np.zeros((2,evsize+rasize))
## CONCATENAR EVOLVENTE E CABEÇA ##
for i in range(0,2): 
    vi1[i] = np.concatenate((xyi1[i],xya1[i]))
    vi2[i] = np.concatenate((xyi2[i],xya2[i]))
## DEFINIR POSIÇÃO DE CONTACTO ##
def cont(rM2):
    alphaM2 = np.arccos(rb2/rM2)
    angR2 = np.tan(alphaM2)-alpha
    T2M = np.sqrt(rM2**2-rb2**2)
    T1M = T1T2 - T2M
    rM1 = np.sqrt(T1M**2+rb1**2)
    alphaM1 = np.arccos(rb1/rM1)
    angR1 = np.pi + np.tan(alphaM1) - alpha
    return alphaM1, alphaM2, angR1, angR2
###############################################################################
## ENGRENAMENTO COM ra1 #######################################################
###############################################################################
ang2s, ang2f = 0, 360

pz = 2*np.pi*r2/np.tan(beta)
bs = pz/(2*np.pi)
ts = np.linspace(0,-2*np.pi,100)
xhb = rb2*np.sin(ts)
yhb = rb2*np.cos(ts)
zh = -bs*ts

xh = r2*np.sin(ts-dm.inv(alpha))
yh = r2*np.cos(ts-dm.inv(alpha))


phib = -b*np.tan(betab)/rb2
vf2 = np.add(dm.matrix(0,vi2), (O2*np.ones((evsize+rasize,2))).T)
vf2t = np.add(dm.matrix(phib,vi2), (O2*np.ones((evsize+rasize,2))).T)

def data_for_cylinder_along_z(center_x,center_y,radius,height_0,height_z):
    z = np.linspace(height_0, height_z, 50)
    theta = np.linspace(0, 2*np.pi, 50)
    theta_grid, z_grid=np.meshgrid(theta, z)
    x_grid = radius*np.cos(theta_grid) + center_x
    y_grid = radius*np.sin(theta_grid) + center_y
    return x_grid,y_grid,z_grid

zv = 'z'
fig = plt.figure(1)
ax = fig.gca(projection='3d')
plt.axis('off')
ax.plot(vf2[0], vf2[1], zs=0, zdir=zv, color='k')
ax.plot(vf2t[0], vf2t[1], zs=b, zdir=zv, color='k')

ax.plot(xhb, yhb, zh, zdir=zv, color='b',label='hélice de base')
ax.plot(xh, yh, zh, zdir=zv, color='r',label='hélice primitiva')

ax.plot([0,0], [0,0], [0,pz], zdir=zv,linestyle='-.',color='k',lw=0.5)

ax.plot(.5*rb1*np.cos(ts), .5*rb1*np.sin(ts), zs=b, zdir=zv,color='k',lw=0.5)
ax.plot(.5*r1*np.cos(ts), .5*r1*np.sin(ts), zs=b, zdir=zv,color='k',lw=0.5)

ax.plot(.5*rb1*np.cos(ts), .5*rb1*np.sin(ts), zs=0,zdir=zv,color='k',lw=0.5)
ax.plot(.5*r1*np.cos(ts), .5*r1*np.sin(ts), zs=0,zdir=zv,color='k',lw=0.5)

ax.plot(.5*rb1*np.cos(ts), .5*rb1*np.sin(ts), zs=pz, zdir=zv,linestyle=':',color='k',lw=0.5)
ax.plot(.5*r1*np.cos(ts), .5*r1*np.sin(ts), zs=pz, zdir=zv,linestyle='--',color='k',lw=0.5)


ax.plot([xhb[0],xhb[-1]], [yhb[0],yhb[-1]], [0,pz], zdir=zv,linestyle='-',color='k',lw=0.5)
ax.plot([xh[0],xh[-1]], [yh[0],yh[-1]], [0,pz], zdir=zv,linestyle='-',color='k',lw=0.5)
# Cylinder
Xc,Yc,Zc = data_for_cylinder_along_z(0.,0.,.5*r1,b,pz)
ax.plot_surface(Xc, Yc, Zc,color='k', alpha=0.05)

# Xc,Yc,Zc = data_for_cylinder_along_z(0.,0.,.5*rb1,pz)
# ax.plot_surface(Xc, Yc, Zc,color='k', alpha=0.15)

Xc,Yc,Zc = data_for_cylinder_along_z(0.,0.,.5*rb1,0,b)
ax.plot_surface(Xc, Yc, Zc,color='b', alpha=0.25)

Xc,Yc,Zc = data_for_cylinder_along_z(0.,0.,.5*r1,0,b)
ax.plot_surface(Xc, Yc, Zc,color='r', alpha=0.25)

# ax.annotate('',[xhb[0],yhb[0]],[xhb[-1],yhb[-1]],arrowprops=dict(arrowstyle="<|-,head_width=0.1,head_length=0.2",\
                            # color='k',shrinkA=0,shrinkB=0,lw=0.5))

ax.legend(loc='best')
# ax.set_xlim(-ra2, ra2)
# ax.set_ylim(0,pz)
# ax.set_zlim(-ra2, ra2)

# ax.set_proj_type('ortho')
ax.view_init(20., 90.)

# plt.show()
plt.savefig('pdf/helice3D.pdf', bbox_inches = 'tight',
    pad_inches = 0.1, transparent=True)