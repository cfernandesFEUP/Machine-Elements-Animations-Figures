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
z = 20.
alpha = np.pi/5
beta = np.pi/9
betab = beta*np.cos(alpha)
x = 0.
hfP = 1.
haP = 1.25
roh = 0.38
## CALCULAR GEOMETRIA ##
r = m*z/2
rl = r + m*x
ra = r + m*(1+x)
rb = r*np.cos(alpha)
p = np.pi*m
pb = p*np.cos(alpha)
s = p/2
O1 = np.array([0,0])
xd = 2*np.pi*r
xb = 2*np.pi*rb
pz = 2*np.pi*r/np.tan(beta)
###############################################################################
## HÉLICES ####################################################################
###############################################################################
plt.figure(1)
fig = plt.gcf()
ax = fig.gca(aspect='equal')
plt.axis('off')
plt.plot([0, xd],[0,0],'k-')
plt.plot([0, 0],[0,pz],'k-')
plt.plot([xd, 0],[0,pz],'k-')

# primitvo 
dm.diml(ax, [0, -.25*pz], [xd, -.25*pz], r'$2\pi\cdot r$')
plt.plot([0, 0], [0, -.3*pz],'k-', lw=0.5)
plt.plot([xd, xd], [0, -.3*pz],'k-', lw=0.5)
# circulo de base
dm.diml(ax, [xd-xb, -.15*pz], [xd, -.15*pz], r'$2\pi\cdot r_b$')
plt.plot([xd-xb, xd-xb], [0, -.2*pz],'k-', lw=0.5)
# passo axial
dm.diml(ax, [-.1*pz, 0], [-.1*pz, pz], r'$p_z$')
plt.plot([0, -.15*pz], [0, 0],'k-', lw=0.5)
plt.plot([0, -.15*pz], [pz, pz],'k-', lw=0.5)


plt.plot([xd, xd-xb],[0,pz],'k--')
plt.plot([xd-xb, xd-xb],[0,pz],'k--')

plt.plot([xd, xd],[0,pz],'k-.', lw=0.5)


plt.plot([xd-xb, 0],[pz,pz],'k--', lw=0.5)

dm.circlewx(ax, .5*pz, 90, 90+beta*180/np.pi, xd, 0, r'$\beta$')
dm.circlewx(ax, .75*pz, 90, 90+betab*180/np.pi, xd, 0, r'$\beta_b$')
plt.savefig('pdf/helices.pdf', bbox_inches = 'tight',
    pad_inches = 0, transparent=True)