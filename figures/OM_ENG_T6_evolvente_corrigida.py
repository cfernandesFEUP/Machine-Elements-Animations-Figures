## IMPORTAR LIBRARIAS ##
import matplotlib.pyplot as plt
import numpy as np
import functions.dimensioning as dm
from matplotlib import rcParams
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Fira Sans']
rcParams['font.size'] = 15
rcParams['legend.fontsize'] = 'medium'
params = {'mathtext.default': 'regular' }  
## DADOS ENGRENAGENS ##
m = 2.
z1 = 20.
z2 = 30.
alpha = np.pi/5
x1 = 0.
x2 = 0.5
n = 100
hfP = 1.
haP = 1.25
roh = 0.38
## CALCULAR GEOMETRIA ##
r1 = m*z1/2
r2 = m*z2/2
ra1 = r1 + m*(1+x1)
ra2 = r2 + m*(1+x2)
rb1 = r1*np.cos(alpha)
rb2 = r2*np.cos(alpha)
a = r1+r2
O1 = [0,a]
O2 = [0,0]
## PLOT FIGURES ##
phi = np.linspace(0,np.pi/2.2,n)
phiC = np.linspace(-np.pi/2,np.pi/2,100)
###############################################################################
## EVOLVENTE DE CÍRCULO #######################################################
###############################################################################
## EQUAÇÕES ##
vi =  dm.ev(rb2,phi)
angR1 = -np.pi/9
vf = dm.matrix(angR1,vi)
## PLOT ##
plt.figure(1)
fig = plt.gcf()
ax = fig.gca(aspect='equal')
plt.axis('off')
plt.plot(vf[0], vf[1],'r')
## raio de base
dm.circle(rb2, 15, 130, 0, 0,'k')
dm.circle(ra2, 60, 110, 0, 0,'k-.')
ra2l = ra2 + x2*m
dm.circle(ra2l, 60, 110, 0, 0,'k')
## cotagem
dm.dimlr(ax,O2, [rb2*np.cos(np.pi/6),rb2*np.sin(np.pi/6)],r'$\mathrm{r_b}$')
dm.dimlr(ax,O2, [ra2*np.cos(np.pi/2.1),ra2*np.sin(np.pi/2.1)],r'$\mathrm{r_a}$')
dm.dimlr(ax,O2, [ra2l*np.cos(np.pi/2.8),ra2l*np.sin(np.pi/2.8)],r"$\mathrm{r_a'}$")
dm.dimt(O2[0],O2[1],'O','left','top')
plt.savefig('pdf/evolvente_corrigida.pdf', bbox_inches = 'tight',
    pad_inches = 0, transparent=True)