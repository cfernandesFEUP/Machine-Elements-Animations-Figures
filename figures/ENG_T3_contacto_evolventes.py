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
O1 = [-a,0]
O2 = [0,0]
## PLOT FIGURES ##
phi = np.linspace(0,np.pi/2.2,n)
###############################################################################
## EVOLVENTE DE CÍRCULO #######################################################
###############################################################################
## Ponto M
M = int(n//1.3)
alphaM = phi[M]
## EQUAÇÕES ##
xy1 =  dm.ev(rb1,phi)
phi1 = np.pi-alphaM
xi1, yi1 = dm.matrix(phi1,xy1)[0], dm.matrix(phi1,xy1)[1]
## derivada da curva involuta
dydx1 = np.diff(yi1)/np.diff(xi1)

xym1 = dm.ev(rb1,alphaM)
xm1, ym1 = dm.matrix(phi1,xym1)[0], dm.matrix(phi1,xym1)[1]
# declive da tangente em M
slope1 = dydx1[M]
b1 = ym1 - slope1*xm1
## reta tangente
l1 = 0.05*rb1
xline1 = np.linspace(xm1-np.sqrt(l1/abs(slope1)), xm1+np.sqrt(l1/abs(slope1)), 2)
yline1 = slope1*xline1 + b1
Tx1 = -a+rb1*np.sin(phi1+alphaM)
Ty1 = rb1*np.cos(phi1+alphaM)
tcx1 = -a+rb1*np.sin(-0.5*alphaM)
tcy1 = rb1*np.cos(-0.5*alphaM)
## PLOT ##
plt.figure(1)
fig = plt.gcf()
ax = fig.gca(aspect='equal')
plt.axis('off')
plt.plot(-a+xi1, yi1,'b')
## raio de base
dm.circle(rb1, 180, 360, -a, 0,'k')
## cotagem
dm.dimln(ax,[-a+xm1,ym1],[Tx1,Ty1],'','k-')
dm.dimln(ax,[-a+xline1[0],yline1[0]],[-a+xline1[-1],yline1[-1]],'','k-')
dm.dimln(ax,[O1[0],O1[1]],[Tx1,Ty1],'','k-')
plt.plot([O1[0],-a+xi1[0]],[O1[1],yi1[0]],'k--',lw=0.5)
dm.dimt(O1[0],O1[1],r'$O_1$','right','bottom')
dm.dimt(-a+xi1[0],yi1[0],r'$Q_1$','right','top')
dm.dimt(-a+xm1,ym1,'X','right','bottom')
dm.dimt(Tx1,Ty1,r'$T_1$','left','bottom')
dm.dimtt(-a+xline1[0],yline1[0],r'$t_1$','right','bottom')
plt.savefig('pdf/evolvente1.pdf', bbox_inches = 'tight',
    pad_inches = 0, transparent=True)



xy2 =  dm.ev(rb2,phi)
phi2 = -alphaM
xi2, yi2 = dm.matrix(phi2,xy2)[0], dm.matrix(phi2,xy2)[1]

## derivada da curva involuta
dydx2 = np.diff(yi2)/np.diff(xi2)
## Ponto M
M = int(n//1.3)
alphaM2 = phi[M]
xym2 = dm.ev(rb2,alphaM2)

xm2, ym2 = dm.matrix(phi2,xym2)[0], dm.matrix(phi2,xym2)[1]
# declive da tangente em M
slope2 = dydx2[M]
b2 = ym2 - slope2*xm2
## reta tangente
l2 = 0.05*rb1
xline2 = np.linspace(xm2-np.sqrt(l2/abs(slope2)), xm2+np.sqrt(l2/abs(slope2)), 2)
yline2 = slope2*xline2 + b2
Tx2 = rb2*np.sin(phi2+alphaM2)
Ty2 = rb2*np.cos(phi2+alphaM2)
tcx2 = rb2*np.sin(-0.5*alphaM2)
tcy2 = rb2*np.cos(-0.5*alphaM2)
## PLOT ##
plt.figure(2)
fig = plt.gcf()
ax = fig.gca(aspect='equal')
plt.axis('off')
plt.plot(xi2, yi2,'r')
## raio de base
dm.circle(rb2, 0, 180, 0, 0,'k')
## cotagem
dm.dimln(ax,[xm2,ym2],[Tx2,Ty2],'','k-')
dm.dimln(ax,[xline2[0],yline2[0]],[xline2[-1],yline2[-1]],'','k-')
dm.dimln(ax,[O2[0],O2[1]],[Tx2,Ty2],'','k-')
plt.plot([O2[0],xi2[0]],[O2[1],yi2[0]],'k--',lw=0.5)
dm.dimt(O2[0],O2[1],r'$O_2$','left','top')
dm.dimt(xi2[0],yi2[0],r'$Q_2$','right','top')
dm.dimt(xm2,ym2,'Y','left','bottom')
dm.dimt(Tx2,Ty2,r'$T_2$','left','bottom')
dm.dimtt(xline2[0],yline2[0],r'$t_2$','left','bottom')
plt.savefig('pdf/evolvente2.pdf', bbox_inches = 'tight',
    pad_inches = 0, transparent=True)

## PLOT ##
plt.figure(3)
fig = plt.gcf()
ax = fig.gca(aspect='equal')
plt.axis('off')
plt.plot(xi2, yi2,'r')
## raio de base
dm.circle(rb2, 0, 180, 0, 0,'k')
## cotagem
dm.dimln(ax,[xm2,ym2],[Tx2,Ty2],'','k-')
dm.dimln(ax,[xline2[0],yline2[0]],[xline2[-1],yline2[-1]],'','k-')
dm.dimln(ax,[O2[0],O2[1]],[Tx2,Ty2],'','k-')
plt.plot([O2[0],xi2[0]],[O2[1],yi2[0]],'k--',lw=0.5)
dm.dimt(O2[0],O2[1],r'$O_2$','left','top')
dm.dimt(xi2[0],yi2[0],r'$Q_2$','right','top')
dm.dimt(xm2,ym2,'Y','left','bottom')
dm.dimt(Tx2,Ty2,r'$T_2$','left','bottom')
dm.dimtt(xline2[0],yline2[0],r'$t_2$','left','bottom')

al = rb1+rb2
at = a+(xm1-xm2-a)
plt.plot(-at+xi1, al+yi1,'b')

## raio de base
dm.circle(rb1, 180, 360, -at, al,'k')
## cotagem
dm.dimln(ax,[-at+xm1,al+ym1],[Tx1-(xm1-xm2-a),al+Ty1],'','k-')
dm.dimln(ax,[-at+xline1[0],al+yline1[0]],[-at+xline1[-1],al+yline1[-1]],'','k-')
dm.dimln(ax,[O1[0]-(xm1-xm2-a),al+O1[1]],[Tx1-(xm1-xm2-a),al+Ty1],'','k-')
plt.plot([O1[0]-(xm1-xm2-a),-at+xi1[0]],[al+O1[1],al+yi1[0]],'k--',lw=0.5)
dm.dimt(O1[0]-(xm1-xm2-a),al+O1[1],r'$O_1$','right','bottom')
dm.dimt(-at+xi1[0],al+yi1[0],r'$Q_1$','right','top')
dm.dimt(-at+xm1,al+ym1,'X','right','bottom')
dm.dimt(Tx1-(xm1-xm2-a),al+Ty1,r'$T_1$','right','bottom')
dm.dimtt(-at+xline1[0],al+yline1[0],r'$t_1$','right','bottom')
plt.savefig('pdf/evolvente12.pdf', bbox_inches = 'tight',
    pad_inches = 0, transparent=True)