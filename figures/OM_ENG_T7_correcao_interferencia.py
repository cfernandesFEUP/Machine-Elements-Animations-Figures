## IMPORTAR LIBRARIAS ##
import numpy as np
import functions.dimensioning as dm
import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Fira Sans']
rcParams['font.size'] = 15
rcParams['legend.fontsize'] = 'medium'
## DADOS ENGRENAGENS ##
m = 2.
z1 = 10.
z2 = 30.
alpha = np.pi/9
zl = 2/(np.sin(alpha)**2)
x1 = (zl-z1)/zl
x2 = 0.
n = 100
hfP = 1.
haP = 1.25
roh = 0.38
## CALCULAR GEOMETRIA ##
r1 = m*z1/2
r2 = m*z2/2
rl1 = r1 + m*x1
rl2 = r2 + m*x2
ra1 = r1 + m*(1+x1)
ra2 = r2 + m*(1+x2)
rb1 = r1*np.cos(alpha)
rb2 = r2*np.cos(alpha)
a = r1+r2
O1 = [r1,r1]
O2 = [0,0]
###############################################################################
## CREMALHEIRA 1 ##############################################################
###############################################################################
vi, vp = dm.rack(m,z2,rb2,x2,hfP,haP,roh,alpha)
Tx = r1+rb1*np.sin(alpha)
Ty = r1-rb1*np.cos(alpha)
p = np.pi*m
xra = hfP*m*np.tan(alpha)
yra = hfP*m
yta = haP*m
## PLOT ##
end = -250
plt.figure(1)
fig = plt.gcf()
ax = fig.gca(aspect='equal')
plt.axis('off')
plt.plot(vi[0,:end], -vi[1,:end]-x1*m,'k')
plt.plot([-r1,1.5*r1], [0,0],'k',lw=0.75)
plt.plot([-r1,1.5*r1], [-x1*m,-x1*m],'k-.',lw=0.75)
plt.plot([-r1,1.5*r1], [m-x1*m,m-x1*m],'r',lw=0.75)
ax.fill_between(vi[0,:end],-2*yra, -vi[1,:end]-x1*m, color='gray')
dm.secant(Tx,r1,Ty,0,-0.05*r1,1.5*r1)
dm.secant(-p/4+p-2*x1*np.tan(alpha),p-p/4+xra-2*x1*np.tan(alpha),0,-yra,0.5*p-2*x1*np.tan(alpha),1.15*p)
dm.secant(-p/4-2*x1*np.tan(alpha),-p/4+xra-2*x1*np.tan(alpha),0,-yra,-0.5*p-2*x1*np.tan(alpha),0.27*p)
## raio de base
dm.circle(rb1, -25, -155, r1, r1,'k--')
dm.circle(r1, -25, -155, r1, r1,'k')
dm.circle(ra1, -25, -155, r1, r1,'k:')
## cotagem
dm.dimlr(ax,O1,[O1[0]+rb1*np.sin(np.pi/3),O1[1]-rb1*np.cos(np.pi/3)],r'$\mathrm{r_b}$')
dm.dimlr(ax,O1,[O1[0]+r1*np.sin(np.pi/4),O1[1]-r1*np.cos(np.pi/4)],r'$\mathrm{r}$')
dm.dimlr(ax,O1,[O1[0]+ra1*np.sin(np.pi/6),O1[1]-ra1*np.cos(np.pi/6)],r'$\mathrm{r_a}$')
dm.dimt(1.35*rl1,-x1*m,'reta de referência\n da cremalheira','left','top')
dm.dimln(ax,O1,[r1,0],'','k-')
dm.dimln(ax,O1, [Tx,Ty],'','k-')
dm.dimln(ax,[-5*p/4-2*x1*np.tan(alpha),0], [-5*p/4-2*x1*np.tan(alpha),-3.2*yra],'','k-')
dm.dimln(ax,[-p/4-2*x1*np.tan(alpha),0], [-p/4-2*x1*np.tan(alpha),-3.2*yra],'','k-')
dm.diml(ax,[-5*p/4-2*x1*np.tan(alpha),-3*yra], [-p/4-2*x1*np.tan(alpha),-3*yra],'p')

dm.dimln(ax,[-p/4-2*x1*np.tan(alpha),0], [-p/4-2*x1*np.tan(alpha),2.2*m],'','k:')
dm.circlec(ax,2*m, 90, 110, -p/4-2*x1*np.tan(alpha), 0,r'$\alpha$')

dm.circlew(ax,2*m, 120, 260, r1, r1,r'$\omega$')

dm.dimlr(ax,[m,-5*yra],[4*m,-5*yra],'')
dm.dimtt(4*m,-5*yra,r'$v=\omega r$','center','bottom')
# plt.plot([O1[0],xi[0]],[O1[1],yi[0]],'k--')
dm.dimt(O1[0],O1[1],'O','left','bottom')
# dm.dimt(xi[0],yi[0],'Q','right','top')
# dm.dimt(xm,ym,'M','left','top')
dm.dimt(Tx,Ty,r'$\mathrm{T}$','left','bottom')
dm.dimt(r1,0,r'$\mathrm{I}$','left','bottom')
plt.savefig('pdf/interferenciaC.pdf', bbox_inches = 'tight',
    pad_inches = 0, transparent=True)
