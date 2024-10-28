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
z2 = 30.
alpha = np.pi/6
x1 = 0.
x2 = 0.
hfP = 1.
haP = 1.25
roh = 0.38
## CALCULAR GEOMETRIA ##
r1, r2 = m*z1/2, m*z2/2
rl1, rl2 = r1 + m*x1, r2 + m*x2
ra1, ra2 = r1 + m*(1+x1), r2 + m*(1+x2)
rb1, rb2 = r1*np.cos(alpha), r2*np.cos(alpha)
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
## RAZÃO DE CONDUÇÃO  #########################################################
###############################################################################
# ang2s, ang2f = 30, 150
# ang1s, ang1f = ang2s + 180, ang2f +180
ang2s, ang2f = 50, 130
ang1s, ang1f = ang2s + 180, ang2f +180
## PLOT ##
plt.figure()
fig = plt.gcf()
ax = fig.gca()
ax.set_aspect("equal")
plt.axis('off')
# raio de base
dm.circle(rb1, ang1s, ang1f, 0, a,'k--')
dm.circle(rb2, ang2s, ang2f, 0, 0,'k--')
## raio de addendum
dm.circle(ra1, ang1s, ang1f, 0, a,'k')
dm.circle(ra2, ang2s, ang2f, 0, 0,'k')
## reta de engrenamento
angA =  -np.pi/2-alphaA1+alpha
angB =  np.pi/2+alpha-alphaA2
T1 = dm.cxy(rb1,0,a,alpha-np.pi/2)
T2 = dm.cxy(rb2,0,0,alpha+np.pi/2)
dm.dimln(ax,T1,T2,'','b-')
dm.dimln(ax,dm.cxy(ra1, O1[0], O1[1], angA),dm.cxy(ra2, O2[0], O2[1], angB),'','r-')
dm.dimln(ax,O2,T2,'','k-')
dm.dimln(ax,O1,T1,'','k-')
dm.dimln(ax,O1,O2,'','k-')
## cotagem
dm.circlewx(ax,0.5*r2, 90+180*alpha/np.pi, 90, 0, 0,r'$\alpha$')
dm.circlewx(ax,0.5*r1, -90+180*alpha/np.pi, -90, 0, a,r'$\alpha$')
dm.circlew(ax,0.15*r1, 150, 290, 0, a,r'$\omega_1$')
dm.circlew(ax,0.15*r1, -150, -290, 0, 0,r'$\omega_2$')
dm.dimt(T1[0],T1[1],r'$\mathrm{T_1}$','left','bottom')
dm.dimt(T2[0],T2[1],r'$\mathrm{T_2}$','right','top')
dm.dimt(O1[0],O1[1],r'$\mathrm{O_1}$','left','bottom')
dm.dimt(O2[0],O2[1],r'$\mathrm{O_2}$','left','top')
dm.dimt(0,r2,r'$\mathrm{I}$','left','bottom')
dm.dimlr(ax,O1,dm.cxy(ra1, O1[0], O1[1], -np.pi/2-alphaA1),r'$\mathrm{r_{a1}}$')
dm.dimlr(ax,O2,dm.cxy(ra2, O2[0], O2[1], 1.5*alphaA2),r'$\mathrm{r_{a2}}$')
dm.dimtarray(dm.cxy(ra1, O1[0], O1[1], angA),r'$\mathrm{A}$','right','top')
dm.dimtarray(dm.cxy(ra2, O2[0], O2[1], angB),r'$\mathrm{B}$','left','bottom')
plt.savefig('pdf/conducao.pdf', bbox_inches = 'tight',
    pad_inches = 0, transparent=True)
###############################################################################
## ENGRENAMENTO ###############################################################
###############################################################################
## ROTAÇÃO ##
rMa2 = (rl2+rb2)/1.95
alphaM1, alphaM2, angR1, angR2 = cont(rMa2)
vf1 = np.add(dm.matrix(angR1,xye1), (O1*np.ones((evsize,2))).T)
vf2 = np.add(dm.matrix(angR2,xye2), (O2*np.ones((evsize,2))).T)
## PLOT ##
plt.figure()
fig = plt.gcf()
ax = fig.gca()
ax.set_aspect("equal")
plt.axis('off')
plt.plot(vf1[0],vf1[1],'b')
plt.plot(vf2[0],vf2[1],'r')
# raio de base
dm.circle(rb1, ang1s, ang1f, 0, a,'k--')
dm.circle(rb2, ang2s, ang2f, 0, 0,'k--')
## reta de engrenamento 
T1 = dm.cxy(rb1,0,a,alpha-np.pi/2)
T2 = dm.cxy(rb2,0,0,alpha+np.pi/2)
dm.dimln(ax,T1,T2,'','k-')
dm.dimln(ax,O2,T2,'','k-')
dm.dimln(ax,O1,T1,'','k-')
dm.dimln(ax,O1,O2,'','k-')
## cotagem
plt.plot([0,.4*r2],[r2,r2],'k',lw=0.5)
dm.circlewxh(ax,0.32*r2, 0, alpha*180/np.pi, 0, r2,r"$\alpha$")
dm.circlewx(ax,0.5*r2, 90+180*alpha/np.pi, 90, 0, 0,r'$\alpha$')
dm.circlewx(ax,0.5*r1, -90+180*alpha/np.pi, -90, 0, a,r'$\alpha$')
dm.circlew(ax,0.15*r1, 150, 290, 0, a,r'$\omega_1$')
dm.circlew(ax,0.15*r1, -150, -290, 0, 0,r'$\omega_2$')
dm.dimlr(ax,O1,[-T1[0],T1[1]],r'$\mathrm{r_{b1}}$')
dm.dimlr(ax,O2,[-T2[0],T2[1]],r'$\mathrm{r_{b2}}$')
dm.dimt(T1[0],T1[1],r'$\mathrm{T_1}$','left','bottom')
dm.dimt(T2[0],T2[1],r'$\mathrm{T_2}$','right','top')
dm.dimt(O1[0],O1[1],r'$\mathrm{O_1}$','left','bottom')
dm.dimt(O2[0],O2[1],r'$\mathrm{O_2}$','left','top')
dm.dimt(0,r2,r'$\mathrm{I}$','left','bottom')
dm.dimt(-rMa2*np.sin(alpha-alphaM2),rMa2*np.cos(alpha-alphaM2),r'$\mathrm{M}$','right','top')
dm.dimt(vf1[0,0],vf1[1,0],r'$\mathrm{Q_1}$','left','bottom')
dm.dimt(vf2[0,0],vf2[1,0],r'$\mathrm{Q_2}$','left','top')
plt.savefig('pdf/engrenamento.pdf', bbox_inches = 'tight',
    pad_inches = 0, transparent=True)
###############################################################################
## ENGRENAMENTO antes de ra1 ##################################################
###############################################################################
## ROTAÇÃO ##
alphaM1, alphaM2, angR1, angR2 = cont(0.97*rA2)
vf1 = np.add(dm.matrix(angR1,vi1), (O1*np.ones((evsize+rasize,2))).T)
vf2 = np.add(dm.matrix(angR2,vi2), (O2*np.ones((evsize+rasize,2))).T)
# ## PLOT ##
plt.figure()
fig = plt.gcf()
ax = fig.gca()
ax.set_aspect("equal")
plt.axis('off')
plt.plot(vf1[0],vf1[1],'b')
plt.plot(vf2[0],vf2[1],'r')
# raio de base
dm.circle(rb1, ang1s, ang1f, 0, a,'k--')
dm.circle(rb2, ang2s, ang2f, 0, 0,'k--')
## raio de addendum
dm.circle(ra1, ang1s, ang1f, 0, a,'k:')
dm.circle(ra2, ang2s, ang2f, 0, 0,'k:')
## reta de engrenamento
T1 = dm.cxy(rb1,0,a,alpha-np.pi/2)
T2 = dm.cxy(rb2,0,0,alpha+np.pi/2)
dm.dimln(ax,T1,T2,'','k-')
dm.dimln(ax,O2,T2,'','k-')
dm.dimln(ax,O1,T1,'','k-')
dm.dimln(ax,O1,O2,'','k-')
## cotagem
dm.circlewx(ax,0.5*r2, 90+180*alpha/np.pi, 90, 0, 0,r'$\alpha$')
dm.circlewx(ax,0.5*r1, -90+180*alpha/np.pi, -90, 0, a,r'$\alpha$')
dm.circlew(ax,0.15*r1, 150, 290, 0, a,r'$\omega_1$')
dm.circlew(ax,0.15*r1, -150, -290, 0, 0,r'$\omega_2$')
dm.dimt(T1[0],T1[1],r'$\mathrm{T_1}$','left','bottom')
dm.dimt(T2[0],T2[1],r'$\mathrm{T_2}$','right','top')
dm.dimt(O1[0],O1[1],r'$\mathrm{O_1}$','left','bottom')
dm.dimt(O2[0],O2[1],r'$\mathrm{O_2}$','left','top')
dm.dimt(0,r2,r'$\mathrm{I}$','right','bottom')
dm.dimt(vf1[0,evsize],vf1[1,evsize],'','left','bottom')
plt.savefig('pdf/engrenamentoAb.pdf', bbox_inches = 'tight',
    pad_inches = 0, transparent=True)
###############################################################################
## ENGRENAMENTO COM ra1 #######################################################
###############################################################################
## ROTAÇÃO ##
alphaM1, alphaM2, angR1, angR2 = cont(rA2)
vf1 = np.add(dm.matrix(angR1,vi1), (O1*np.ones((evsize+rasize,2))).T)
vf2 = np.add(dm.matrix(angR2,vi2), (O2*np.ones((evsize+rasize,2))).T)
# ## PLOT ##
plt.figure()
fig = plt.gcf()
ax = fig.gca()
ax.set_aspect("equal")
plt.axis('off')
plt.plot(vf1[0],vf1[1],'b')
plt.plot(vf2[0],vf2[1],'r')
# raio de base
dm.circle(rb1, ang1s, ang1f, 0, a,'k--')
dm.circle(rb2, ang2s, ang2f, 0, 0,'k--')
## raio de addendum
dm.circle(ra1, ang1s, ang1f, 0, a,'k:')
dm.circle(ra2, ang2s, ang2f, 0, 0,'k:')
## reta de engrenamento
T1 = dm.cxy(rb1,0,a,alpha-np.pi/2)
T2 = dm.cxy(rb2,0,0,alpha+np.pi/2)
dm.dimln(ax,T1,T2,'','k-')
dm.dimln(ax,O2,T2,'','k-')
dm.dimln(ax,O1,T1,'','k-')
dm.dimln(ax,O1,O2,'','k-')
## cotagem
dm.circlewx(ax,0.5*r2, 90+180*alpha/np.pi, 90, 0, 0,r'$\alpha$')
dm.circlewx(ax,0.5*r1, -90+180*alpha/np.pi, -90, 0, a,r'$\alpha$')
dm.circlew(ax,0.15*r1, 150, 290, 0, a,r'$\omega_1$')
dm.circlew(ax,0.15*r1, -150, -290, 0, 0,r'$\omega_2$')
dm.dimt(T1[0],T1[1],r'$\mathrm{T_1}$','left','bottom')
dm.dimt(T2[0],T2[1],r'$\mathrm{T_2}$','right','top')
dm.dimt(O1[0],O1[1],r'$\mathrm{O_1}$','left','bottom')
dm.dimt(O2[0],O2[1],r'$\mathrm{O_2}$','left','top')
dm.dimt(0,r2,r'$\mathrm{I}$','right','bottom')
dm.dimt(vf1[0,evsize],vf1[1,evsize],r'$A$','left','bottom')
plt.savefig('pdf/engrenamentoA.pdf', bbox_inches = 'tight',
    pad_inches = 0, transparent=True)
###############################################################################
## ENGRENAMENTO COM ra2 #######################################################
###############################################################################
## ROTAÇÃO ##
alphaM1, alphaM2, angR1, angR2 = cont(ra2)
vf1 = np.add(dm.matrix(angR1,vi1), (O1*np.ones((evsize+rasize,2))).T)
vf2 = np.add(dm.matrix(angR2,vi2), (O2*np.ones((evsize+rasize,2))).T)
# ## PLOT ##
plt.figure()
fig = plt.gcf()
ax = fig.gca()
ax.set_aspect("equal")
plt.axis('off')
plt.plot(vf1[0],vf1[1],'b')
plt.plot(vf2[0],vf2[1],'r')
# raio de base
dm.circle(rb1, ang1s, ang1f, 0, a,'k--')
dm.circle(rb2, ang2s, ang2f, 0, 0,'k--')
## raio de addendum
dm.circle(ra1, ang1s, ang1f, 0, a,'k:')
dm.circle(ra2, ang2s, ang2f, 0, 0,'k:')
## reta de engrenamento
T1 = dm.cxy(rb1,0,a,alpha-np.pi/2)
T2 = dm.cxy(rb2,0,0,alpha+np.pi/2)
dm.dimln(ax,T1,T2,'','k-')
dm.dimln(ax,O2,T2,'','k-')
dm.dimln(ax,O1,T1,'','k-')
dm.dimln(ax,O1,O2,'','k-')
## cotagem
dm.circlewx(ax,0.5*r2, 90+180*alpha/np.pi, 90, 0, 0,r'$\alpha$')
dm.circlewx(ax,0.5*r1, -90+180*alpha/np.pi, -90, 0, a,r'$\alpha$')
dm.circlew(ax,0.15*r1, 150, 290, 0, a,r'$\omega_1$')
dm.circlew(ax,0.15*r1, -150, -290, 0, 0,r'$\omega_2$')
dm.dimt(T1[0],T1[1],r'$\mathrm{T_1}$','left','bottom')
dm.dimt(T2[0],T2[1],r'$\mathrm{T_2}$','right','top')
dm.dimt(O1[0],O1[1],r'$\mathrm{O_1}$','left','bottom')
dm.dimt(O2[0],O2[1],r'$\mathrm{O_2}$','left','top')
dm.dimt(0,r2,r'$\mathrm{I}$','right','bottom')
dm.dimt(vf2[0,evsize],vf2[1,evsize],r'$B$','left','bottom')
plt.savefig('pdf/engrenamentoB.pdf', bbox_inches = 'tight',
    pad_inches = 0, transparent=True)
###############################################################################
## ENGRENAMENTO depois de ra2 #################################################
###############################################################################
## ROTAÇÃO ##
alphaM1, alphaM2, angR1, angR2 = cont(1.03*ra2)
vf1 = np.add(dm.matrix(angR1,vi1), (O1*np.ones((evsize+rasize,2))).T)
vf2 = np.add(dm.matrix(angR2,vi2), (O2*np.ones((evsize+rasize,2))).T)
# ## PLOT ##
plt.figure()
fig = plt.gcf()
ax = fig.gca()
ax.set_aspect("equal")
plt.axis('off')
plt.plot(vf1[0],vf1[1],'b')
plt.plot(vf2[0],vf2[1],'r')
# raio de base
dm.circle(rb1, ang1s, ang1f, 0, a,'k--')
dm.circle(rb2, ang2s, ang2f, 0, 0,'k--')
## raio de addendum
dm.circle(ra1, ang1s, ang1f, 0, a,'k:')
dm.circle(ra2, ang2s, ang2f, 0, 0,'k:')
## reta de engrenamento
T1 = dm.cxy(rb1,0,a,alpha-np.pi/2)
T2 = dm.cxy(rb2,0,0,alpha+np.pi/2)
dm.dimln(ax,T1,T2,'','k-')
dm.dimln(ax,O2,T2,'','k-')
dm.dimln(ax,O1,T1,'','k-')
dm.dimln(ax,O1,O2,'','k-')
## cotagem
dm.circlewx(ax,0.5*r2, 90+180*alpha/np.pi, 90, 0, 0,r'$\alpha$')
dm.circlewx(ax,0.5*r1, -90+180*alpha/np.pi, -90, 0, a,r'$\alpha$')
dm.circlew(ax,0.15*r1, 150, 290, 0, a,r'$\omega_1$')
dm.circlew(ax,0.15*r1, -150, -290, 0, 0,r'$\omega_2$')
dm.dimt(T1[0],T1[1],r'$\mathrm{T_1}$','left','bottom')
dm.dimt(T2[0],T2[1],r'$\mathrm{T_2}$','right','top')
dm.dimt(O1[0],O1[1],r'$\mathrm{O_1}$','left','bottom')
dm.dimt(O2[0],O2[1],r'$\mathrm{O_2}$','left','top')
dm.dimt(0,r2,r'$\mathrm{I}$','right','bottom')
dm.dimt(vf2[0,evsize],vf2[1,evsize],'','left','bottom')
plt.savefig('pdf/engrenamentoBb.pdf', bbox_inches = 'tight',
    pad_inches = 0, transparent=True)
###############################################################################
## ENGRENAMENTO COM ra1 + pb ##################################################
###############################################################################
## ROTAÇÃO ##
alphaM1, alphaM2, angR1, angR2 = cont(rA2)
pb0 = angR2
vf1 = np.add(dm.matrix(angR1,vi1), (O1*np.ones((evsize+rasize,2))).T)
vf2 = np.add(dm.matrix(angR2,vi2), (O2*np.ones((evsize+rasize,2))).T)
alphaM1, alphaM2, angR1, angR2 = cont(rAN2)
vf11 = np.add(dm.matrix(angR1,vi1), (O1*np.ones((evsize+rasize,2))).T)
vf21 = np.add(dm.matrix(angR2,vi2), (O2*np.ones((evsize+rasize,2))).T)
pb1 = angR2
plt.figure()
fig = plt.gcf()
ax = fig.gca()
ax.set_aspect("equal")
plt.axis('off')
plt.plot(vf1[0],vf1[1],'b')
plt.plot(vf2[0],vf2[1],'r')
plt.plot(vf11[0],vf11[1],'b')
plt.plot(vf21[0],vf21[1],'r')
# raio de base
dm.circle(rb1, ang1s, ang1f, 0, a,'k--')
dm.circle(rb2, ang2s, ang2f, 0, 0,'k--')
## raio de addendum
dm.circle(ra1, ang1s, ang1f, 0, a,'k:')
dm.circle(ra2, ang2s, ang2f, 0, 0,'k:')
## reta de engrenamento
T1 = dm.cxy(rb1,0,a,alpha-np.pi/2)
T2 = dm.cxy(rb2,0,0,alpha+np.pi/2)
dm.dimln(ax,T1,T2,'','k-')
dm.dimln(ax,O2,T2,'','k-')
dm.dimln(ax,O1,T1,'','k-')
dm.dimln(ax,O1,O2,'','k-')
## cotagem
dm.circlewxb(ax,rb2, 90-180*pb1/np.pi, 90-180*pb0/np.pi, 0, 0,r'$p_b$')
pA = dm.cxy(ra1, O1[0], O1[1], angA)
T2Apb = T2A + pb
angApb = np.arctan(T2Apb/rb2)
rM2 = np.sqrt(T2Apb**2+rb2**2)
pB = dm.cxy(rM2, O2[0], O2[1], np.pi/2+alpha-angApb)
dm.diml(ax, pA, pB,r'$p_b$')
dm.circlewx(ax,0.5*r2, 90+180*alpha/np.pi, 90, 0, 0,r'$\alpha$')
dm.circlewx(ax,0.5*r1, -90+180*alpha/np.pi, -90, 0, a,r'$\alpha$')
dm.circlew(ax,0.15*r1, 150, 290, 0, a,r'$\omega_1$')
dm.circlew(ax,0.15*r1, -150, -290, 0, 0,r'$\omega_2$')
dm.dimt(T1[0],T1[1],r'$\mathrm{T_1}$','left','bottom')
dm.dimt(T2[0],T2[1],r'$\mathrm{T_2}$','right','top')
dm.dimt(O1[0],O1[1],r'$\mathrm{O_1}$','left','bottom')
dm.dimt(O2[0],O2[1],r'$\mathrm{O_2}$','left','top')
dm.dimt(0,r2,r'$\mathrm{I}$','right','bottom')
dm.dimtarray(dm.cxy(ra1, O1[0], O1[1], angA),r'$\mathrm{A}$','right','top')
dm.dimtarray(dm.cxy(ra2, O2[0], O2[1], angB),r'$\mathrm{B}$','left','bottom')
plt.savefig('pdf/engrenamentoC.pdf', bbox_inches = 'tight',
    pad_inches = 0, transparent=True)
###############################################################################
## ENGRENAMENTO COM ra1 > pb ##################################################
###############################################################################
## ROTAÇÃO ##
alphaM1, alphaM2, angR1, angR2 = cont(ra2)
vf1 = np.add(dm.matrix(angR1,vi1), (O1*np.ones((evsize+rasize,2))).T)
vf2 = np.add(dm.matrix(angR2,vi2), (O2*np.ones((evsize+rasize,2))).T)
pb1 = angR2
alphaM1, alphaM2, angR1, angR2 = cont(0.98*rA2)
vf11 = np.add(dm.matrix(angR1,vi1), (O1*np.ones((evsize+rasize,2))).T)
vf21 = np.add(dm.matrix(angR2,vi2), (O2*np.ones((evsize+rasize,2))).T)
pb0 = angR2
plt.figure()
fig = plt.gcf()
ax = fig.gca()
ax.set_aspect("equal")
plt.axis('off')
plt.plot(vf1[0],vf1[1],'b')
plt.plot(vf2[0],vf2[1],'r')
plt.plot(vf11[0],vf11[1],'b:')
plt.plot(vf21[0],vf21[1],'r:')
# raio de base
dm.circle(rb1, ang1s, ang1f, 0, a,'k--')
dm.circle(rb2, ang2s, ang2f, 0, 0,'k--')
## raio de addendum
dm.circle(ra1, ang1s, ang1f, 0, a,'k:')
dm.circle(ra2, ang2s, ang2f, 0, 0,'k:')
## reta de engrenamento
T1 = dm.cxy(rb1,0,a,alpha-np.pi/2)
T2 = dm.cxy(rb2,0,0,alpha+np.pi/2)
dm.dimln(ax,T1,T2,'','k-')
dm.dimln(ax,O2,T2,'','k-')
dm.dimln(ax,O1,T1,'','k-')
dm.dimln(ax,O1,O2,'','k-')
## cotagem
dm.circlewxb(ax,rb2, 90-180*pb1/np.pi, 90-180*pb0/np.pi, 0, 0,r'$p_b$')
pA = dm.cxy(0.98*rA2, O2[0], O2[1], np.pi/2+alpha-np.arccos(rb2/(0.98*rA2)))
dm.diml(ax, pA, dm.cxy(ra2, O2[0], O2[1], angB),r'$p_b$')
dm.circlewx(ax,0.5*r2, 90+180*alpha/np.pi, 90, 0, 0,r'$\alpha$')
dm.circlewx(ax,0.5*r1, -90+180*alpha/np.pi, -90, 0, a,r'$\alpha$')
dm.circlew(ax,0.15*r1, 150, 290, 0, a,r'$\omega_1$')
dm.circlew(ax,0.15*r1, -150, -290, 0, 0,r'$\omega_2$')
dm.dimt(T1[0],T1[1],r'$\mathrm{T_1}$','left','bottom')
dm.dimt(T2[0],T2[1],r'$\mathrm{T_2}$','right','top')
dm.dimt(O1[0],O1[1],r'$\mathrm{O_1}$','left','bottom')
dm.dimt(O2[0],O2[1],r'$\mathrm{O_2}$','left','top')
dm.dimt(0,r2,r'$\mathrm{I}$','right','bottom')
# dm.dimt(vf2[0,evsize],vf2[1,evsize],r'$A$','left','bottom')
dm.dimtarray(dm.cxy(ra1, O1[0], O1[1], angA),r'$\mathrm{A}$','left','center')
dm.dimtarray(dm.cxy(ra2, O2[0], O2[1], angB),r'$\mathrm{B}$','left','bottom')
plt.savefig('pdf/engrenamentoD.pdf', bbox_inches = 'tight',
    pad_inches = 0, transparent=True)
###############################################################################
## RAZÃO DE CONDUÇÃO  #########################################################
###############################################################################
## PLOT ##
plt.figure()
fig = plt.gcf()
ax = fig.gca()
ax.set_aspect("equal")
plt.axis('off')
# raio de base
dm.circle(rb1, ang1s, ang1f, 0, a,'k--')
dm.circle(rb2, ang2s, ang2f, 0, 0,'k--')
## raio de addendum
dm.circle(ra1, ang1s, ang1f, 0, a,'k')
dm.circle(ra2, ang2s, ang2f, 0, 0,'k')
## reta de engrenamento
angA =  -np.pi/2-alphaA1+alpha
angB =  np.pi/2+alpha-alphaA2
T1 = dm.cxy(rb1,0,a,alpha-np.pi/2)
T2 = dm.cxy(rb2,0,0,alpha+np.pi/2)
dm.dimln(ax,T1,T2,'','b-')
dm.dimln(ax,dm.cxy(ra1, O1[0], O1[1], angA),dm.cxy(ra2, O2[0], O2[1], angB),'','r-')
dm.dimln(ax,O2,T2,'','k-')
dm.dimln(ax,O1,T1,'','k-')
dm.dimln(ax,O1,O2,'','k-')
## cotagem
dm.circlewx(ax,0.5*r2, 90+180*alpha/np.pi, 90, 0, 0,r'$\alpha$')
dm.circlewx(ax,0.5*r1, -90+180*alpha/np.pi, -90, 0, a,r'$\alpha$')
dm.circlew(ax,0.15*r1, 150, 290, 0, a,r'$\omega_1$')
dm.circlew(ax,0.15*r1, -150, -290, 0, 0,r'$\omega_2$')
dm.dimt(T1[0],T1[1],r'$\mathrm{T_1}$','left','bottom')
dm.dimt(T2[0],T2[1],r'$\mathrm{T_2}$','right','top')
dm.dimt(O1[0],O1[1],r'$\mathrm{O_1}$','left','bottom')
dm.dimt(O2[0],O2[1],r'$\mathrm{O_2}$','left','top')
dm.dimt(0,r2,r'$\mathrm{I}$','left','bottom')
dm.dimlr(ax,O1,dm.cxy(ra1, O1[0], O1[1], -np.pi/2-alphaA1),r'$\mathrm{r_{a1}}$')
dm.dimlr(ax,O2,dm.cxy(ra2, O2[0], O2[1], 1.5*alphaA2),r'$\mathrm{r_{a2}}$')
dm.dimtarray(dm.cxy(ra1, O1[0], O1[1], angA),r'$\mathrm{A}$','right','top')
dm.dimtarray(dm.cxy(ra2, O2[0], O2[1], angB),r'$\mathrm{B}$','left','bottom')
plt.savefig('pdf/conducaoS.pdf', bbox_inches = 'tight',
    pad_inches = 0, transparent=True)
###############################################################################
## TP1  #######################################################################
###############################################################################
ang2s, ang2f = 50, 130
ang1s, ang1f = ang2s + 180, ang2f +180
## PLOT ##
plt.figure()
fig = plt.gcf()
ax = fig.gca()
ax.set_aspect("equal")
plt.axis('off')
# raio de base
dm.circle(rb1, ang1s, ang1f, 0, a,'k--')
dm.circle(rb2, ang2s, ang2f, 0, 0,'k--')
## raio de addendum
dm.circle(r1, ang1s, ang1f, 0, a,'k')
dm.circle(r2, ang2s, ang2f, 0, 0,'k')
## reta de engrenamento
T1 = dm.cxy(rb1,0,a,alpha-np.pi/2)
T2 = dm.cxy(rb2,0,0,alpha+np.pi/2)
dm.dimln(ax,T1,T2,'','k-')
dm.dimln(ax,O2,T2,'','k-')
dm.dimln(ax,O1,T1,'','k-')
dm.dimln(ax,O1,O2,'','k-')
## cotagem
dm.circlew(ax,0.15*r1, 150, 290, 0, a,r'$\omega_1$')
dm.circlew(ax,0.15*r1, -150, -290, 0, 0,r'$\omega_2$')
dm.dimt(T1[0],T1[1],r'$\mathrm{T_1}$','left','bottom')
dm.dimt(T2[0],T2[1],r'$\mathrm{T_2}$','right','top')
dm.dimt(O1[0],O1[1],r'$\mathrm{O_1}$','left','bottom')
dm.dimt(O2[0],O2[1],r'$\mathrm{O_2}$','left','top')
dm.dimt(0,r2,r'$\mathrm{I}$','left','bottom')
dm.dimlr(ax,O1,dm.cxy(r1, O1[0], O1[1], -np.pi/2-alphaA1),r'$\mathrm{r_{1}}$')
dm.dimlr(ax,O2,dm.cxy(r2, O2[0], O2[1], 1.5*alphaA2),r'$\mathrm{r_{2}}$')
dm.dimlr(ax,O1,dm.cxy(rb1, O1[0], O1[1], -np.pi/2-.4*alphaA1),r'$\mathrm{r_{b1}}$')
dm.dimlr(ax,O2,dm.cxy(rb2, O2[0], O2[1], 2*alphaA2),r'$\mathrm{r_{b2}}$')
plt.savefig('pdf/TP1_primitivos.pdf', bbox_inches = 'tight',
    pad_inches = 0, transparent=True)
###############################################################################
## RAIO ATIVO DE PÉ ###########################################################
###############################################################################
## ROTAÇÃO ##
alphaM1, alphaM2, angR1, angR2 = cont(rA2)
vf1 = np.add(dm.matrix(angR1,vi1), (O1*np.ones((evsize+rasize,2))).T)
vf2 = np.add(dm.matrix(angR2,vi2), (O2*np.ones((evsize+rasize,2))).T)
# ## PLOT ##
plt.figure()
fig = plt.gcf()
ax = fig.gca()
ax.set_aspect("equal")
plt.axis('off')
# plt.plot(vf1[0],vf1[1],'b')
plt.plot(vf2[0],vf2[1],'r')
# raio de base
dm.circle(rb1, ang1s, ang1f, 0, a,'k--')
dm.circle(rb2, ang2s, ang2f, 0, 0,'k--')
## raio de addendum
dm.circle(ra1, ang1s, ang1f, 0, a,'k:')
dm.circle(ra2, ang2s, ang2f, 0, 0,'k:')
## reta de engrenamento
T1 = dm.cxy(rb1,0,a,alpha-np.pi/2)
T2 = dm.cxy(rb2,0,0,alpha+np.pi/2)
dm.dimln(ax,T1,T2,'','k-')
dm.dimln(ax,O2,T2,'','k-')
dm.dimln(ax,O1,T1,'','k-')
dm.dimln(ax,O1,O2,'','k-')
dm.dimlr(ax,O2,[vf1[0,evsize],vf1[1,evsize]],r'$r_A$')
## cotagem
dm.circlewx(ax,0.25*r2, 90+180*alpha/np.pi, 90, 0, 0,r'$\alpha$')
dm.circlewx(ax,0.5*r1, -90+180*alpha/np.pi, -90, 0, a,r'$\alpha$')
dm.circlew(ax,0.15*r1, 150, 290, 0, a,r'$\omega_1$')
dm.circlew(ax,0.15*r1, -150, -290, 0, 0,r'$\omega_2$')
dm.dimt(T1[0],T1[1],r'$\mathrm{T_1}$','left','bottom')
dm.dimt(T2[0],T2[1],r'$\mathrm{T_2}$','right','top')
dm.dimt(O1[0],O1[1],r'$\mathrm{O_1}$','left','bottom')
dm.dimt(O2[0],O2[1],r'$\mathrm{O_2}$','left','top')
dm.dimt(0,r2,r'$\mathrm{I}$','right','bottom')
dm.dimt(vf1[0,evsize],vf1[1,evsize],r'$A$','left','bottom')
plt.savefig('pdf/raioAtivo.pdf', bbox_inches = 'tight',
    pad_inches = 0, transparent=True)
###############################################################################
## ENGRENAMENTO: AÇÃO CONJUGADA ###############################################
###############################################################################
## ROTAÇÃO ##
rMa2 = (rl2+rb2)/1.95
alphaM1, alphaM2, angR1, angR2 = cont(rMa2)
vf1 = np.add(dm.matrix(angR1,xye1), (O1*np.ones((evsize,2))).T)
vf2 = np.add(dm.matrix(angR2,xye2), (O2*np.ones((evsize,2))).T)
## PLOT ##
plt.figure()
fig = plt.gcf()
ax = fig.gca()
ax.set_aspect("equal")
plt.axis('off')
plt.plot(vf1[0],vf1[1],'b')
plt.plot(vf2[0],vf2[1],'r')
# raio de base
dm.circle(rb1, ang1s, ang1f, 0, a,'k--')
dm.circle(rb2, ang2s, ang2f, 0, 0,'k--')
dm.circle(r1, ang1s, ang1f, 0, a,'k-.')
dm.circle(r2, 1.2*ang2s, .9*ang2f, 0, 0,'k-.')
## reta de engrenamento 
T1 = dm.cxy(rb1,0,a,alpha-np.pi/2)
T2 = dm.cxy(rb2,0,0,alpha+np.pi/2)
dm.dimln(ax,T1,T2,'','k-')
dm.dimln(ax,O2,T2,'','k-')
dm.dimln(ax,O1,T1,'','k-')
dm.dimln(ax,O1,O2,'','k-')
## cotagem
dm.circlew(ax,0.15*r1, 150, 290, 0, a,r'$\omega_1$')
dm.circlew(ax,0.15*r1, -150, -290, 0, 0,r'$\omega_2$')
plt.text(r1*np.cos(.97*ang1f),a+r1*np.sin(0.97*ang1f),r'$\mathrm{C_1}$')
plt.text(r2*np.cos(1.2*ang2f),-r2*np.sin(1.2*ang2f),r'$\mathrm{C_2}$')
dm.dimt(T1[0],T1[1],r'$\mathrm{T_1}$','left','bottom')
dm.dimt(T2[0],T2[1],r'$\mathrm{T_2}$','right','top')
dm.dimt(O1[0],O1[1],r'$\mathrm{O_1}$','left','bottom')
dm.dimt(O2[0],O2[1],r'$\mathrm{O_2}$','left','top')
dm.dimt(0,r2,r'$\mathrm{I}$','left','bottom')
dm.dimt(-rMa2*np.sin(alpha-alphaM2),rMa2*np.cos(alpha-alphaM2),r'$\mathrm{M}$','right','top')
plt.savefig('pdf/acao_conjugada.pdf', bbox_inches = 'tight',
    pad_inches = 0, transparent=True)