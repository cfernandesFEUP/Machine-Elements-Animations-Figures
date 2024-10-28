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
phi = np.linspace(0,np.pi/2.2,n)
phiC = np.linspace(-np.pi/2,np.pi/2,100)
###############################################################################
## RODAS DE ATRITO ############################################################
###############################################################################
## PLOT ##
plt.figure(1)
fig = plt.gcf()
ax = fig.gca(aspect='equal')
plt.axis('off')
## circulos primitivos
dm.circle(r1, 0, 360, 0, a, 'k')
dm.circle(r2, 0, 360, 0, 0, 'k')


w1i = dm.cxy(0.6*r1, 0, a, -np.pi/2+np.pi/3)
w1e = dm.cxy(0.6*r1, 0, a, -np.pi/2-np.pi/3)
F1 = 0.6*r1- np.sqrt((0.6*r1)**2-(w1i[0])**2)
w2i = dm.cxy(0.6*r1, 0, 0, np.pi/2+np.pi/3)
w2e = dm.cxy(0.6*r1, 0, 0, np.pi/2-np.pi/3)
# dm.dimangr(ax,w1i[0],w1i[1],w1e[0],w1e[1],'arc3,rad='+str(F1/w1i[0]))
# dm.dimangr(ax,w2i[0],w2i[1],w2e[0],w2e[1],'arc3,rad='+str(F1/w1i[0]))

## cotagem
dm.circlew(ax,0.5*r1, -120, -50, 0, a,r'$\omega_1$')
dm.circlew(ax,0.5*r2, 120, 50, 0, 0,r'$\omega_2$')
dm.dimt(0,r2,r'$\mathrm{I}$','left','top')
dm.dimt(O1[0], O1[1],r'$\mathrm{O_1}$','left','bottom')
dm.dimt(O2[0], O2[1],r'$\mathrm{O_2}$','left','top')
dm.dimlr(ax,O2,[-r2*np.sin(np.pi/3),r2*np.cos(np.pi/3)],r'$\mathrm{r_2}$')
dm.dimlr(ax,O1,[-r1*np.sin(np.pi/3),a+r1*np.cos(np.pi/3)],r'$\mathrm{r_1}$')
plt.savefig('pdf/rodas_atrito.pdf', bbox_inches = 'tight',
    pad_inches = 0, transparent=True)
###############################################################################
## EVOLVENTE DE CÍRCULO #######################################################
###############################################################################
## EQUAÇÕES ##
xi =  dm.ev(rb2,phi)[0]
yi =  dm.ev(rb2,phi)[1]
## derivada da curva involuta
dydx = np.diff(yi)/np.diff(xi)
## Ponto M
M = int(n//1.3)
alphaM = phi[M]
xm = dm.ev(rb2,alphaM)[0]
ym = dm.ev(rb2,alphaM)[1]
# declive da tangente em M
slope = dydx[M]
b = ym - slope*xm
## reta tangente
l = rb2
xline = np.linspace(xm-np.sqrt(l/slope), xm+np.sqrt(l/slope), 2)
yline = slope*xline + b
Tx = rb2*np.sin(alphaM)
Ty = rb2*np.cos(alphaM)
tcx = rb2*np.sin(-0.5*alphaM)
tcy = rb2*np.cos(-0.5*alphaM)
## PLOT ##
plt.figure(2)
fig = plt.gcf()
ax = fig.gca(aspect='equal')
plt.axis('off')
plt.plot(xi, yi,'r')
## raio de base
dm.circle(rb2, 0, 130, 0, 0,'k')
## cotagem
dm.dimlr(ax,[0,0],[tcx,tcy],'rb')
dm.dimln(ax,[xm,ym],[Tx,Ty],r'$t_r$','k-')
dm.dimln(ax,[xline[0],yline[0]],[xline[-1],yline[-1]],'','k-')
dm.dimln(ax,[O2[0],O2[1]],[Tx,Ty],'','k-')
plt.plot([O2[0],xi[0]],[O2[1],yi[0]],'k--',lw=0.5)
dm.dimt(O2[0],O2[1],'O','left','top')
dm.dimt(xi[0],yi[0],'Q','right','top')
dm.dimt(xm,ym,'M','center','bottom')
dm.dimt(Tx,Ty,'T','left','bottom')
dm.dimtt(xline[0],yline[0],r'$t_e$','right','bottom')
plt.text(xi[-1], yi[-1],'evolvente',fontsize=16, horizontalalignment='left',
      verticalalignment='center')
plt.savefig('pdf/evolvente.pdf', bbox_inches = 'tight',
    pad_inches = 0, transparent=True)
###############################################################################
## ÂNGULO DE INCIDÊNCIA #######################################################
###############################################################################
## PLOT ##
plt.figure(3)
fig = plt.gcf()
ax = fig.gca(aspect='equal')
plt.axis('off')
plt.plot(xi, yi,'r')
## raio de base
dm.circle(rb2, 0, 130, 0, 0,'k')
## cotagem
dm.dimlr(ax,[0,0],[tcx,tcy],r'$\mathrm{r_b}$')
dm.dimln(ax,[xm,ym],[Tx,Ty],r'$t_r$','k-')
dm.dimln(ax,[xline[0],yline[0]],[xline[-1],yline[-1]],'','k-')
dm.dimln(ax,[O2[0],O2[1]],[Tx,Ty],'','k-')
plt.plot([O2[0],xi[0]],[O2[1],yi[0]],'k--',lw=0.5)
plt.plot([O2[0],xm],[O2[1],ym],'k--',lw=0.5)
dm.dimt(O2[0],O2[1],r'$\mathrm{O}$','left','top')
dm.dimt(xi[0],yi[0],r'$\mathrm{Q}$','right','top')

dm.dimt(xm,ym,r'$\mathrm{M}$','center','bottom')
dm.dimt(Tx,Ty,r'$\mathrm{T}$','left','bottom')
dm.dimtt(xline[0],yline[0],r'$t_e$','right','bottom')
th1 = np.arccos(xm/np.sqrt(xm**2+ym**2))
th2 = np.arccos(Tx/rb2)
xang1 = 0
yang1 = 0.7*rb2*np.sin(90)
xang2 = 0.55*rb2*np.cos(th1+(th2-th1)/2)
yang2 = 0.55*rb2*np.sin(th1+(th2-th1)/2)
dm.dimang(ax,0,0.6*rb2,0.6*rb2*np.cos(th1),0.6*rb2*np.sin(th1),"arc3,rad=0.1")
dm.dimang(ax,0.5*rb2*np.cos(th1),0.5*rb2*np.sin(th1),0.5*rb2*np.cos(th2),0.5*rb2*np.sin(th2),"arc3,rad=0.25")
dm.dimt(rb2*np.cos(th1),rb2*np.sin(th1),r'$\mathrm{X}$','left','top')
dm.dimtt(xang1,yang1,r'inv$\alpha_M$','center','bottom')
dm.dimtt(xang2,yang2,r'$\alpha_M$','left','center')
plt.savefig('pdf/incidencia.pdf', bbox_inches = 'tight',
    pad_inches = 0, transparent=True)
###############################################################################
## EQUAÇÕES PARAMÉTRICAS ######################################################
###############################################################################
xxi =  dm.ev(rb2,phi)[1]
yyi =  dm.ev(rb2,phi)[0]
## derivada da curva involuta
dyydxx = np.diff(yyi)/np.diff(xxi)
## Ponto M
xxm = dm.ev(rb2,alphaM)[1]
yym = dm.ev(rb2,alphaM)[0]
# declive da tangente em M
slopep = dyydxx[M]
bp = yym - slopep*xxm
## reta tangente
xlinep = np.linspace(xxm-np.sqrt(l/slopep), xxm+np.sqrt(l/slopep), 2)
ylinep = slopep*xlinep + bp
Txx = rb2*np.sin(np.pi/2-alphaM)
Tyy = rb2*np.cos(np.pi/2-alphaM)
tcxx = rb2*np.sin(alphaM)
tcyy = rb2*np.cos(alphaM)
## PLOT ##
plt.figure(4)
fig = plt.gcf()
ax = fig.gca(aspect='equal')
plt.axis('off')
plt.plot(xxi, yyi,'r')
## raio de base
dm.circle(rb2, 0, 90, 0, 0,'k')
## cotagem
dm.dimlr(ax,[0,0],[1.5*rb2,0],'')
dm.dimlr(ax,O2,[Txx,Tyy],r'$\mathrm{r_b}$')
dm.dimtt(1.5*rb2,0,r'$x$','left','top')
dm.dimlr(ax,[0,0],[0,1.25*rb2],'')
dm.dimtt(0,1.26*rb2,r'$y$','right','top')

dm.dimln(ax,[xxm,yym],[Txx,Tyy],'','k-')
# plt.plot([O2[0],Txx],[O2[1],Tyy],'k-')
plt.plot([O2[0],xxm],[O2[1],yym],'k--',lw=0.5)
# dm.dimt(O2[0],O2[1],'O','left','top')
dm.dimt(xxm,yym,'M','left','top')
dm.dimt(Txx,Tyy,'T','left','bottom')
th1 = np.arccos(xxm/np.sqrt(xxm**2+yym**2))
th2 = np.arccos(Txx/rb2)
xang1 = 0.55*rb2*np.cos(th1/2)
yang1 = 0.55*rb2*np.sin(th1/2)
xang2 = 0.55*rb2*np.cos(th1+(th2-th1)/2)
yang2 = 0.55*rb2*np.sin(th1+(th2-th1)/2)
dm.dimang(ax,0.5*rb2,0,0.5*rb2*np.cos(th1),0.5*rb2*np.sin(th1),"arc3,rad=-0.1")
dm.dimang(ax,0.5*rb2*np.cos(th1),0.5*rb2*np.sin(th1),0.5*rb2*np.cos(th2),0.5*rb2*np.sin(th2),"arc3,rad=-0.25")

dm.dimtt(xang1,yang1,r'inv $\alpha_M$','left','center')
dm.dimtt(xang2,yang2,r'$\alpha_M$','left','center')
# ax.add_patch(Arc((0, 0), .9*rb2, .9*rb2,
#                  theta1=0, theta2=th1, edgecolor='k', lw=0.5))
# ax.add_patch(Arc((0, 0), .9*rb2, .9*rb2,
#                   theta1=th1, theta2=th2, edgecolor='k', lw=0.5))
plt.savefig('pdf/evol_parametrica.pdf', bbox_inches = 'tight',
    pad_inches = 0, transparent=True)
###############################################################################
## DUPLA EVOLVENTE DE CÍRCULO #################################################
###############################################################################
## EQUAÇÕES ##
phi1 = np.linspace(0,np.pi/3,n)
xi =  dm.ev(rb2,phi1)[0]
yi =  dm.ev(rb2,phi1)[1]
## derivada da curva involuta
dydx = np.diff(yi)/np.diff(xi)
## Ponto M
M = int(n//1.4)
alphaM = phi1[M]
xm = dm.ev(rb2,alphaM)[0]
ym = dm.ev(rb2,alphaM)[1]
# declive da tangente em M
slope = dydx[M]
b = ym - slope*xm
## reta tangente
l = rb2
xline = np.linspace(xm-np.sqrt(l/slope), xm+np.sqrt(l/slope), 2)
yline = slope*xline + b
Tx = rb2*np.sin(alphaM)
Ty = rb2*np.cos(alphaM)
tcx = rb2*np.sin(-0.5*alphaM)
tcy = rb2*np.cos(-0.5*alphaM)
phi2 = np.linspace(0,np.pi/3.5,n)
vi = dm.ev(rb2,phi2)
pb = np.pi*m*np.cos(alpha)
angPb =  pb/rb2
vf = dm.matrix(angPb,vi)
## PLOT ##
plt.figure(5)
fig = plt.gcf()
ax = fig.gca(aspect='equal')
plt.axis('off')
plt.plot(xi, yi,'r')
plt.plot(vf[0], vf[1],'r')
## raio de base
dm.circle(rb2, 20, 120, 0, 0,'k')
## cotagem
dm.dimlr(ax,[0,0],[tcx,tcy],'rb')
dm.dimln(ax,[xm,ym],[Tx,Ty],'','k-')
dm.dimln(ax,[xline[0],yline[0]],[xline[-1],yline[-1]],'','k-')
dm.dimln(ax,[O2[0],O2[1]],[Tx,Ty],'','k-')
dm.dimt(O2[0],O2[1],'O','left','top')
dm.dimt(xi[0],yi[0],r'$Q_2$','right','top')
dm.dimt(vf[0,0],vf[1,0],r'$Q_1$','right','top')
dm.dimt(xm,ym,r'$M_2$','center','bottom')
## M1 ##
rM2 = np.sqrt(xm**2+ym**2)
TM2 = np.sqrt(rM2**2-rb2**2)
rM1 = np.sqrt((TM2-pb)**2 + rb2**2)
alpha2 = np.arccos(rb2/rM2)
alphaT = np.pi/2-dm.inv(alpha2)- alpha2
alpha1 = alphaT + np.arccos(rb2/rM1)
dm.dimt(rM1*np.cos(alpha1),rM1*np.sin(alpha1),r'$M_1$','center','bottom')
dm.dimt(Tx,Ty,'T','left','bottom')
dm.dimtt(xline[0],yline[0],'','right','bottom')
plt.savefig('pdf/dupla_evolvente.pdf', bbox_inches = 'tight',
    pad_inches = 0, transparent=True)