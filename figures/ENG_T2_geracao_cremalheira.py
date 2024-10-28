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
plt.plot(vi[0,:end], -vi[1,:end],'k')
plt.plot([-0.5*r1,1.5*r1], [0,0],'k',lw=0.75)
ax.fill_between(vi[0,:end],-2*yra, -vi[1,:end], color='gray')
dm.secant(Tx,r1,Ty,0,-0.05*r1,1.5*r1)
dm.secant(-p/4+p,p-p/4+xra,0,-yra,0.5*p,1.15*p)
dm.secant(-p/4,-p/4+xra,0,-yra,-0.5*p,0.27*p)
## raio de base
dm.circle(rb1, -25, -155, r1, r1,'k--')
dm.circle(r1, -25, -155, r1, r1,'k')
dm.circle(ra1, -25, -155, r1, r1,'k:')
## cotagem
dm.dimlr(ax,O1,[O1[0]+rb1*np.sin(np.pi/3),O1[1]-rb1*np.cos(np.pi/3)],r'$\mathrm{r_b}$')
dm.dimlr(ax,O1,[O1[0]+r1*np.sin(np.pi/4),O1[1]-r1*np.cos(np.pi/4)],r'$\mathrm{r}$')
dm.dimlr(ax,O1,[O1[0]+ra1*np.sin(np.pi/6),O1[1]-ra1*np.cos(np.pi/6)],r'$\mathrm{r_a}$')
dm.dimt(1.35*r1,0,'reta de referência\n da cremalheira','left','top')
dm.dimln(ax,O1,[r1,0],'','k-')
dm.dimln(ax,O1, [Tx,Ty],'','k-')
dm.dimln(ax,[-5*p/4,0], [-5*p/4,-3.2*yra],'','k-')
dm.dimln(ax,[-p/4,0], [-p/4,-3.2*yra],'','k-')
dm.diml(ax,[-5*p/4,-3*yra], [-p/4,-3*yra],'p')

dm.dimln(ax,[-p/4,0], [-p/4,2.2*m],'','k:')
dm.circlec(ax,2*m, 90, 110, -p/4, 0,r'$\alpha$')

dm.circlew(ax,2*m, 120, 260, r1, r1,r'$\omega$')

dm.dimlr(ax,[m,-5*yra],[4*m,-5*yra],'')
dm.dimtt(4*m,-5*yra,r'$v=\omega r$','center','bottom')
# plt.plot([O1[0],xi[0]],[O1[1],yi[0]],'k--')
dm.dimt(O1[0],O1[1],'O','left','bottom')
# dm.dimt(xi[0],yi[0],'Q','right','top')
# dm.dimt(xm,ym,'M','left','top')
dm.dimt(Tx,Ty,r'$\mathrm{T}$','left','bottom')
dm.dimt(r1,0,r'$\mathrm{I}$','left','bottom')
plt.savefig('pdf/cremalheiraA.pdf', bbox_inches = 'tight',
    pad_inches = 0, transparent=True)
###############################################################################
## CREMALHEIRA 2 ##############################################################
###############################################################################
## evolvente
T1B = np.sqrt(ra1**2-rb1**2)
alphaA = T1B/rb1
phiE1 = np.linspace(0,-alphaA,n)
xi1 =  dm.ev(rb1,phiE1)[0]
yi1 =  dm.ev(rb1,phiE1)[1]
vi1 =  np.array([xi1,yi1])
T1M = r1*np.sin(alpha) + (r1-p+p/4)*np.cos(alpha)
alpha1 = np.arctan(T1M/rb1)
phi1 = np.pi -alpha + alpha1 + dm.inv(alpha1)
vf1 = dm.matrix(phi1,vi1)
vM =  dm.ev(rb1,-(alpha1+dm.inv(alpha1)))
vfM = dm.matrix(phi1,vM)

phiEc = np.linspace(-alphaA,-np.pi/2.6,n)
xi1c =  dm.ev(rb1,phiEc)[0]
yi1c =  dm.ev(rb1,phiEc)[1]
vi1c =  np.array([xi1c,yi1c])
vfc1 = dm.matrix(phi1,vi1c)


begin = 250
## PLOT ##
plt.figure(2)
fig = plt.gcf()
ax = fig.gca(aspect='equal')
plt.axis('off')
plt.plot(vi[0,begin:end], -vi[1,begin:end],'k')
plt.plot([-0.5*r1,1.5*r1], [0,0],'k',lw=0.75)
ax.fill_between(vi[0,begin:end],-2*m*hfP, -vi[1,begin:end], color='gray')
plt.plot(vf1[0]+r1,vf1[1]+r1,'r')
plt.plot(vfc1[0]+r1,vfc1[1]+r1,'k', lw=0.5)
## secantes
dm.secant(Tx,r1,Ty,0,0,1.5*r1)
dm.secant(-p/4,-p/4+xra,0,-yra,-0.5*p,0.27*p)
dm.secant(3*p/4,3*p/4+xra,0,-yra,0.5*p,1.15*p)
## raio de base
dm.circle(rb1, -25, -155, r1, r1,'k--')
dm.circle(r1, -25, -155, r1, r1,'k')
dm.circle(ra1, -25, -155, r1, r1,'k:')
## cotagem
dm.dimlr(ax,O1,[O1[0]+rb1*np.sin(np.pi/3),O1[1]-rb1*np.cos(np.pi/3)],r'$\mathrm{r_b}$')
dm.dimlr(ax,O1,[O1[0]+r1*np.sin(np.pi/4),O1[1]-r1*np.cos(np.pi/4)],r'$\mathrm{r}$')
dm.dimlr(ax,O1,[O1[0]+ra1*np.sin(np.pi/6),O1[1]-ra1*np.cos(np.pi/6)],r'$\mathrm{r_a}$')
dm.dimln(ax,O1,[r1,0],'','k-')
dm.dimln(ax,O1, [Tx,Ty],'','k-')
dm.dimln(ax,[-p/4,0], [-p/4,2.2*m],'','k:')
dm.circlec(ax,2*m, 90, 110, -p/4, 0,r'$\alpha$')
dm.circlew(ax,2*m, 120, 260, r1, r1,r'$\omega$')
dm.dimlr(ax,[m,-5*yra],[4*m,-5*yra],'')
dm.dimtt(4*m,-5*yra,r'$v=\omega r$','center','bottom')
# plt.plot([O1[0],xi[0]],[O1[1],yi[0]],'k--')
dm.dimt(O1[0],O1[1],r'$\mathrm{O}$','left','bottom')
dm.dimt(Tx,Ty,r'$\mathrm{T}$','left','bottom')
dm.dimt(r1,0,r'$\mathrm{I}$','left','bottom')
dm.dimt(vf1[0,0]+r1,vf1[1,0]+r1,r'$\mathrm{Q}$','left','bottom')
dm.dimt(vfM[0]+r1,vfM[1]+r1,r'$\mathrm{M}$','left','bottom')
plt.savefig('pdf/cremalheiraB.pdf', bbox_inches = 'tight',
    pad_inches = 0, transparent=True)
###############################################################################
## CREMALHEIRA 3 ##############################################################
###############################################################################
## evolvente
phiE2 = np.linspace(0,-alphaA,n)
xi2 =  dm.ev(rb1,phiE2)[0]
yi2 =  dm.ev(rb1,phiE2)[1]
vi2 =  np.array([xi2,yi2])
T1M2 = r1*np.sin(alpha) + (r1+p/4)*np.cos(alpha)
alpha12 = np.arctan(T1M2/rb1)
phi12 = np.pi -alpha + alpha12 + dm.inv(alpha12)
vf12 = dm.matrix(phi12,vi2)
vM1 =  dm.ev(rb1,-(alpha12+dm.inv(alpha12)))
vfM1 = dm.matrix(phi12,vM1)

phiEc2 = np.linspace(-alphaA,-np.pi/2.1,n)
xi2c =  dm.ev(rb1,phiEc2)[0]
yi2c =  dm.ev(rb1,phiEc2)[1]
vi2c =  np.array([xi2c,yi2c])

vfc12 = dm.matrix(phi12,vi2c)

## PLOT ##
plt.figure(3)
fig = plt.gcf()
ax = fig.gca(aspect='equal')
plt.axis('off')
plt.plot(vi[0,begin:end], -vi[1,begin:end],'k')
plt.plot([-0.5*r1,1.5*r1], [0,0],'k',lw=0.75)
ax.fill_between(vi[0,begin:end],-2*m*hfP, -vi[1,begin:end], color='gray')
plt.plot(vf1[0]+r1,vf1[1]+r1,'r')
plt.plot(vfc1[0]+r1,vfc1[1]+r1,'k', lw=0.5)
plt.plot(vf12[0]+r1,vf12[1]+r1,'r')
plt.plot(vfc12[0]+r1,vfc12[1]+r1,'k', lw=0.5)
## secantes
dm.secant(Tx,r1,Ty,0,0,1.5*r1)
dm.secant(-p/4,-p/4+xra,0,-yra,-0.5*p,0.27*p)
dm.secant(3*p/4,3*p/4+xra,0,-yra,0.5*p,1.15*p)
## raio de base
dm.circle(rb1, -25, -155, r1, r1,'k--')
dm.circle(r1, -25, -155, r1, r1,'k')
dm.circle(ra1, -25, -155, r1, r1,'k:')
dm.circleb(rb1, (np.pi/2-phi1)*180/np.pi, (np.pi/2-phi12)*180/np.pi, r1, r1,'b')
plt.plot([vfM[0]+r1,vfM1[0]+r1],[vfM[1]+r1,vfM1[1]+r1],'b')
## cotagem
dm.dimlr(ax,O1,[O1[0]+rb1*np.sin(np.pi/3),O1[1]-rb1*np.cos(np.pi/3)],r'$\mathrm{r_b}$')
dm.dimlr(ax,O1,[O1[0]+r1*np.sin(np.pi/4),O1[1]-r1*np.cos(np.pi/4)],r'$\mathrm{r}$')
dm.dimlr(ax,O1,[O1[0]+ra1*np.sin(np.pi/6),O1[1]-ra1*np.cos(np.pi/6)],r'$\mathrm{r_a}$')
dm.dimln(ax,O1,[r1,0],'','k-')
dm.dimln(ax,O1, [Tx,Ty],'','k-')

## passo p
dm.dimln(ax,[-3*p/4,0], [-3*p/4,-3*yra],'','k-')
dm.dimln(ax,[p/4,0], [p/4,-3*yra],'','k-')
dm.diml(ax,[-3*p/4,-2.5*yra], [p/4,-2.5*yra],'p')

dm.circlew(ax,2*m, 120, 260, r1, r1,r'$\omega$')
dm.dimlr(ax,[m,-5*yra],[4*m,-5*yra],'')
dm.dimtt(4*m,-5*yra,r'$v=\omega r$','center','bottom')
# plt.plot([O1[0],xi[0]],[O1[1],yi[0]],'k--')
dm.dimt(O1[0],O1[1],r'$\mathrm{O}$','left','bottom')
dm.dimt(Tx,Ty,r'$\mathrm{T}$','left','bottom')
dm.dimt(r1,0,r'$\mathrm{I}$','left','bottom')
dm.dimt(vf1[0,0]+r1,vf1[1,0]+r1,r'$\mathrm{Q}$','left','bottom')
dm.dimt(vf12[0,0]+r1,vf12[1,0]+r1,r'$\mathrm{Q_1}$','left','bottom')
dm.dimt(vfM[0]+r1,vfM[1]+r1,r'$\mathrm{M}$','left','bottom')
dm.dimt(vfM1[0]+r1,vfM1[1]+r1,r'$\mathrm{M_1}$','right','top')
plt.savefig('pdf/cremalheiraC.pdf', bbox_inches = 'tight',
    pad_inches = 0, transparent=True)