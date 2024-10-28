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
alpha = np.pi/9
x1 = 0.
x2 = 0.
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
T2B = np.sqrt(ra2**2-rb2**2)
alphaA = T2B/rb2
phi = np.linspace(0,alphaA,n)
###############################################################################
## ESPESSURA DENTE ############################################################
###############################################################################
## EQUAÇÕES ##
vli =  dm.ev(rb2,phi)
p = np.pi*m
phil = -(p/(4*rb2)+dm.inv(alpha))
vfl = dm.matrix(phil,vli)
## Ponto M
M = int(n//1.2)
alphaM = phi[M]
xym = dm.ev(rb2,alphaM)
vfm = dm.matrix(phil,xym)
xyY = dm.ev(rb2,alpha+dm.inv(alpha))
vfY = dm.matrix(phil,xyY)
xa = np.linspace(vfl[0,-1],-vfl[0,-1],10)
ya = np.sqrt(ra2**2-xa**2)

invx = np.concatenate((vfl[0],xa,-vfl[0,::-1]))
invy = np.concatenate((vfl[1],ya,vfl[1,::-1]))
## PLOT ##
plt.figure(9)
fig = plt.gcf()
ax = fig.gca(aspect='equal')
plt.axis('off')
plt.plot(invx, invy,'r')
## raio de base
dm.circle(rb2, 60, 120, 0, 0,'k--')
dm.circle(r2, 70, 110, 0, 0,'k')
dm.circle(ra2, 80, 100, 0, 0,'k:')
## cotagem
dm.dimlr(ax,dm.cxy(0.5*rb2,0,0,np.pi/2+1.2*alpha),dm.cxy(rb2,0,0,np.pi/2+1.2*alpha),r'$\mathrm{r_b}$')
dm.dimlr(ax,dm.cxy(0.5*rb2,0,0,np.pi/2+.8*alpha),dm.cxy(r2,0,0,np.pi/2+0.8*alpha),r'$\mathrm{r}$')
dm.dimlr(ax,dm.cxy(0.5*rb2,0,0,np.pi/2+.5*alpha),dm.cxy(ra2,0,0,np.pi/2+0.5*alpha),r'$\mathrm{r_a}$')
# dm.dimln(ax,[xline[0],yline[0]],[xline[-1],yline[-1]],'','k-')
# dm.dimln(ax,[O2[0],O2[1]],[Tx,Ty],'','k-')
# plt.plot([O2[0],xi[0]],[O2[1],yi[0]],'k--',lw=0.5)
plt.plot([O2[0],vfm[0]],[O2[1],vfm[1]],'k--',lw=0.5)
# dm.dimt(O2[0],O2[1],r'$\mathrm{O}$','left','top')
dm.dimtd(vfl[0,0],vfl[1,0],r'$\mathrm{Q}$','right','top')
dm.dimtd(-vfl[0,0],vfl[1,0],r"$\mathrm{Q'}$",'left','top')
dm.dimtd(vfm[0],vfm[1],r'$\mathrm{M}$','right','bottom')
dm.dimtd(-vfm[0],vfm[1],r"$\mathrm{M'}$",'left','bottom')
dm.dimtd(vfY[0],vfY[1],r'$\mathrm{Y}$','right','top')
dm.dimtd(-vfY[0],vfY[1],r"$\mathrm{Y'}$",'left','top')
# dm.dimtt(xline[0],yline[0],r'$t_e$','right','bottom')
# th1 = np.arccos(xm/np.sqrt(xm**2+ym**2))
# th2 = np.arccos(Tx/rb2)
# xang1 = 0
# yang1 = 0.7*rb2*np.sin(90)
# xang2 = 0.55*rb2*np.cos(th1+(th2-th1)/2)
# yang2 = 0.55*rb2*np.sin(th1+(th2-th1)/2)
# dm.dimang(ax,0,0.6*rb2,0.6*rb2*np.cos(th1),0.6*rb2*np.sin(th1),"arc3,rad=0.1")
# dm.dimang(ax,0.5*rb2*np.cos(th1),0.5*rb2*np.sin(th1),0.5*rb2*np.cos(th2),0.5*rb2*np.sin(th2),"arc3,rad=0.25")
# dm.dimt(rb2*np.cos(th1),rb2*np.sin(th1),r'$\mathrm{X}$','left','top')
# dm.dimtt(xang1,yang1,r'inv$\alpha_M$','center','bottom')
# dm.dimtt(xang2,yang2,r'$\alpha_M$','left','center')
plt.savefig('pdf/espessura.pdf', bbox_inches = 'tight',
    pad_inches = 0)

