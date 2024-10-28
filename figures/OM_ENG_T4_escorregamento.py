## IMPORTAR LIBRARIAS ##
import matplotlib.pyplot as plt
import numpy as np
import functions.dimensioning as dm
from matplotlib import rcParams
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Fira Sans']
rcParams['font.size'] = 12
rcParams['legend.fontsize'] = 'medium'
params = {'mathtext.default': 'regular' }  
## DADOS ENGRENAGENS ##
m = 2.
z1 = 20.
z2 = 30.
alpha = np.pi/7
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

ang2s, ang2f = 40, 140
ang1s, ang1f = ang2s + 180, ang2f +180
###############################################################################
## VELOCIDADES PRIMITIVO ######################################################
###############################################################################
## ROTAÇÃO ##
rM2 = rl2
alphaM1, alphaM2, angR1, angR2 = cont(rM2)
vf1 = np.add(dm.matrix(angR1,vi1), (O1*np.ones((evsize+rasize,2))).T)
vf2 = np.add(dm.matrix(angR2,vi2), (O2*np.ones((evsize+rasize,2))).T)
## VELOCIDADE 
wm = -0.75
xm = rM2*np.sin(alphaM2-alpha)
ym = rM2*np.cos(alphaM2-alpha)
O2M = np.array([xm, ym, 0])
O1M = np.array([0, a, 0]) - O2M
w2 = np.array([0, 0, wm])
w1 = w2*z2/z1
vm1 = np.cross(w1,O1M)
vm2 = np.cross(w2,O2M) 
# ## PLOT ##
plt.figure(1)
fig = plt.gcf()
ax = fig.gca(aspect='equal')
plt.axis('off')
## reta de engrenamento
T1 = dm.cxy(rb1,0,a,alpha-np.pi/2)
T2 = dm.cxy(rb2,0,0,alpha+np.pi/2)

## Perfis
plt.plot(vf1[0],vf1[1],'b')
plt.plot(vf2[0],vf2[1],'r')
# raio de base
dm.circle(rb1, ang1s, ang1f, 0, a,'k--')
dm.circle(rb2, ang2s, ang2f, 0, 0,'k--')

xl=20
dm.dimln(ax,[T1[0]+xl*np.cos(alpha),T1[1]+xl*np.sin(alpha)], T2,'','k-')
dm.dimln(ax,O2,T2,'','k-')
dm.dimln(ax,O1,T1,'','k-')
dm.dimln(ax,O1,O2,'','k-')
dm.dimln(ax,O2,[xm,ym],'','k-')
dm.dimln(ax,O1,[xm,ym],'','k-')
## Speed
dm.dimll(ax,[xm,ym],[xm+vm1[0],ym+vm1[1]],r'$\mathrm{v_{I1}}=\mathrm{v_{I2}}$','left','center','b',1.2)
dm.dimll(ax,[xm,ym],[xm+vm2[0],ym+vm2[1]],'','right','center','r',1)


V1x = np.sqrt(vm1[0]**2+vm1[1]**2)*np.cos(alphaM1)
V1y = np.sqrt(vm1[0]**2+vm1[1]**2)*np.sin(alphaM1)
V2x = np.sqrt(vm2[0]**2+vm2[1]**2)*np.cos(alphaM2)
V2y = np.sqrt(vm2[0]**2+vm2[1]**2)*np.sin(alphaM2)

dm.dimll(ax,[xm,ym],[xm+V1x*np.cos(alpha),ym+V1x*np.sin(alpha)],r'$\mathrm{v_{b1}}=\mathrm{v_{b2}}$','left','center','b',1.2)
dm.dimll(ax,[xm,ym],[xm+V2x*np.cos(alpha),ym+V2x*np.sin(alpha)],'','center','bottom','r',1)
dm.dimll(ax,[xm,ym],[xm+V2y*np.cos(alpha-np.pi/2),ym+V2y*np.sin(alpha-np.pi/2)],r'$\mathrm{v_{r1}}=\mathrm{v_{r2}}$','left','center','r',1)
dm.dimll(ax,[xm,ym],[xm+V1y*np.cos(alpha-np.pi/2),ym+V1y*np.sin(alpha-np.pi/2)],'','center','bottom','b',1.2)

## cotagem
dm.circlew(ax,0.15*r1, 150, 290, 0, a,r'$\omega_1$')
dm.circlew(ax,0.15*r1, -150, -290, 0, 0,r'$\omega_2$')
alpha1 = np.pi/13
# dm.circlewx(ax,0.75*r2, 90+180*alpha/np.pi, 90-180*(alphaM2-alpha)/np.pi, 0, 0,r'$\alpha_2$')
dm.circlewx(ax,0.5*r2, 90+180*alpha/np.pi, 90, 0, 0,r'$\alpha$')

# dm.circlewx(ax,0.75*r1, -90+180*alpha/np.pi, -90+180*(alpha-alphaM1)/np.pi, 0, a,r'$\alpha_1$')
dm.circlewx(ax,0.5*r1, -90+180*alpha/np.pi, -90, 0, a,r'$\alpha$')
dm.dimt(T1[0],T1[1],r'$\mathrm{T_1}$','left','bottom')
dm.dimt(T2[0],T2[1],r'$\mathrm{T_2}$','right','top')
dm.dimt(O1[0],O1[1],r'$\mathrm{O_1}$','left','bottom')
dm.dimt(O2[0],O2[1],r'$\mathrm{O_2}$','left','top')
dm.dimt(0,r2,r'$\mathrm{I}$','right','bottom')
# dm.dimt(xm,ym,r'$\mathrm{M}$','right','bottom')
plt.savefig('pdf/escorregamentoI.pdf', bbox_inches = 'tight',
    pad_inches = 0, transparent=True)

###############################################################################
## ENGRENAMENTO COM ra1 #######################################################
###############################################################################
## ROTAÇÃO ##
# rM2 = rl2 + m/2 
alphaMA1, alphaMA2, angRA1, angRA2 = cont(rA2)
vfA1 = np.add(dm.matrix(angRA1,vi1), (O1*np.ones((evsize+rasize,2))).T)
vfA2 = np.add(dm.matrix(angRA2,vi2), (O2*np.ones((evsize+rasize,2))).T)
# ## PLOT ##
plt.figure(2)
fig = plt.gcf()
ax = fig.gca(aspect='equal')
plt.axis('off')
plt.plot(vf1[0],vf1[1],'b',lw=0.5)
plt.plot(vf2[0],vf2[1],'r',lw=0.5)

#ponto A
plt.plot(vfA1[0],vfA1[1],'b',lw=0.5)
plt.plot(vfA2[0],vfA2[1],'r',lw=0.5)

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
dm.dimt(T1[0],T1[1],r'$\mathrm{T_1}$','left','bottom')
dm.dimt(T2[0],T2[1],r'$\mathrm{T_2}$','right','top')
dm.dimt(O1[0],O1[1],r'$\mathrm{O_1}$','left','bottom')
dm.dimt(O2[0],O2[1],r'$\mathrm{O_2}$','left','top')
dm.dimt(0,r2,r'$\mathrm{I}$','right','bottom')
angA =  -np.pi/2-alphaA1+alpha
angB =  np.pi/2+alpha-alphaA2
dm.dimtarray(dm.cxy(ra1, O1[0], O1[1], angA),r'$\mathrm{A}$','right','top')
dm.dimtarray(dm.cxy(ra2, O2[0], O2[1], angB),r'$\mathrm{B}$','left','bottom')
plt.savefig('pdf/escorregamento.pdf', bbox_inches = 'tight',
    pad_inches = 0, transparent=True)



###############################################################################
## ENGRENAMENTO COM ra1 #######################################################
###############################################################################
## ROTAÇÃO ##
# rM2 = rl2 + m/2 
alphaMA1, alphaMA2, angRA1, angRA2 = cont(ra2)
vfA1 = np.add(dm.matrix(angRA1,vi1), (O1*np.ones((evsize+rasize,2))).T)
vfA2 = np.add(dm.matrix(angRA2,vi2), (O2*np.ones((evsize+rasize,2))).T)
# ## PLOT ##
plt.figure(3)
fig = plt.gcf()
ax = fig.gca(aspect='equal')
plt.axis('off')
plt.plot(vf1[0],vf1[1],'b',lw=0.5)
plt.plot(vf2[0],vf2[1],'r',lw=0.5)

#ponto A
plt.plot(vfA1[0],vfA1[1],'b',lw=0.5)
plt.plot(vfA2[0],vfA2[1],'r',lw=0.5)

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
dm.dimt(T1[0],T1[1],r'$\mathrm{T_1}$','left','bottom')
dm.dimt(T2[0],T2[1],r'$\mathrm{T_2}$','right','top')
dm.dimt(O1[0],O1[1],r'$\mathrm{O_1}$','left','bottom')
dm.dimt(O2[0],O2[1],r'$\mathrm{O_2}$','left','top')
dm.dimt(0,r2,r'$\mathrm{I}$','right','bottom')
angA =  -np.pi/2-alphaA1+alpha
angB =  np.pi/2+alpha-alphaA2
dm.dimtarray(dm.cxy(ra1, O1[0], O1[1], angA),r'$\mathrm{A}$','right','top')
dm.dimtarray(dm.cxy(ra2, O2[0], O2[1], angB),r'$\mathrm{B}$','left','bottom')
plt.savefig('pdf/escorregamento.pdf', bbox_inches = 'tight',
    pad_inches = 0, transparent=True)


###############################################################################
## VELOCIDADES LINEARES #######################################################
###############################################################################
## ROTAÇÃO ##
rM2 = rl2 + m/2 
alphaM1, alphaM2, angR1, angR2 = cont(rM2)
vf1 = np.add(dm.matrix(angR1,vi1), (O1*np.ones((evsize+rasize,2))).T)
vf2 = np.add(dm.matrix(angR2,vi2), (O2*np.ones((evsize+rasize,2))).T)
## VELOCIDADE 
xm = rM2*np.sin(alphaM2-alpha)
ym = rM2*np.cos(alphaM2-alpha)
O2M = np.array([xm, ym, 0])
O1M = np.array([0, a, 0]) - O2M
w2 = np.array([0, 0, wm])
w1 = w2*z2/z1
vm1 = np.cross(w1,O1M)
vm2 = np.cross(w2,O2M) 
# ## PLOT ##
plt.figure(4)
fig = plt.gcf()
ax = fig.gca(aspect='equal')
plt.axis('off')
## reta de engrenamento
T1 = dm.cxy(rb1,0,a,alpha-np.pi/2)
T2 = dm.cxy(rb2,0,0,alpha+np.pi/2)


## Perfis
plt.plot(vf1[0],vf1[1],'b')
plt.plot(vf2[0],vf2[1],'r')
# raio de base
dm.circle(rb1, ang1s, ang1f, 0, a,'k--')
dm.circle(rb2, ang2s, ang2f, 0, 0,'k--')

xl=20
dm.dimln(ax,[T1[0]+xl*np.cos(alpha),T1[1]+xl*np.sin(alpha)], T2,'','k-')
dm.dimln(ax,O2,T2,'','k-')
dm.dimln(ax,O1,T1,'','k-')
dm.dimln(ax,O1,O2,'','k-')
dm.dimln(ax,O2,[xm,ym],'','k-')
dm.dimln(ax,O1,[xm,ym],'','k-')
## Speed
dm.dimll(ax,[xm,ym],[xm+vm1[0],ym+vm1[1]],r'$\mathrm{v_{M1}}$','left','center','b',1.2)
dm.dimll(ax,[xm,ym],[xm+vm2[0],ym+vm2[1]],r'$\mathrm{v_{M2}}$','left','center','r',1)


V1x = np.sqrt(vm1[0]**2+vm1[1]**2)*np.cos(alphaM1)
V1y = np.sqrt(vm1[0]**2+vm1[1]**2)*np.sin(alphaM1)
V2x = np.sqrt(vm2[0]**2+vm2[1]**2)*np.cos(alphaM2)
V2y = np.sqrt(vm2[0]**2+vm2[1]**2)*np.sin(alphaM2)

plt.plot([xm+V2x*np.cos(alpha),xm+vm2[0]], [ym+V2x*np.sin(alpha),ym+vm2[1]],'k--',lw=0.5)
# plt.plot([xm+V2y*np.cos(alpha-np.pi/2),xm+vm2[0]], [ym+V2y*np.sin(alpha-np.pi/2),ym+vm2[1]],'k--',lw=0.5)

## cotagem
dm.circlew(ax,0.15*r1, 150, 290, 0, a,r'$\omega_1$')
dm.circlew(ax,0.15*r1, -150, -290, 0, 0,r'$\omega_2$')
alpha1 = np.pi/13
dm.circlewx(ax,0.75*r2, 90+180*alpha/np.pi, 90-180*(alphaM2-alpha)/np.pi, 0, 0,r'$\alpha_2$')
dm.circlewx(ax,0.5*r2, 90+180*alpha/np.pi, 90, 0, 0,r'$\alpha$')

dm.circlewx(ax,0.75*r1, -90+180*alpha/np.pi, -90+180*(alpha-alphaM1)/np.pi, 0, a,r'$\alpha_1$')
dm.circlewx(ax,0.5*r1, -90+180*alpha/np.pi, -90, 0, a,r'$\alpha$')
dm.dimt(T1[0],T1[1],r'$\mathrm{T_1}$','left','bottom')
dm.dimt(T2[0],T2[1],r'$\mathrm{T_2}$','right','top')
dm.dimt(O1[0],O1[1],r'$\mathrm{O_1}$','left','bottom')
dm.dimt(O2[0],O2[1],r'$\mathrm{O_2}$','left','top')
dm.dimt(xm,ym,r'$\mathrm{M}$','right','center')
plt.savefig('pdf/vM.pdf', bbox_inches = 'tight',
    pad_inches = 0, transparent=True)
###############################################################################
## VELOCIDADES TANGENCIAIS BASE ###############################################
###############################################################################
## PLOT ##
plt.figure(5)
fig = plt.gcf()
ax = fig.gca(aspect='equal')
plt.axis('off')
## reta de engrenamento
T1 = dm.cxy(rb1,0,a,alpha-np.pi/2)
T2 = dm.cxy(rb2,0,0,alpha+np.pi/2)

xl=20
dm.dimln(ax,[T1[0]+xl*np.cos(alpha),T1[1]+xl*np.sin(alpha)], T2,'','k-')
dm.dimln(ax,O2,T2,'','k-')
dm.dimln(ax,O1,T1,'','k-')
dm.dimln(ax,O1,O2,'','k-')
dm.dimln(ax,O2,[xm,ym],'','k-')
dm.dimln(ax,O1,[xm,ym],'','k-')
## Speed
dm.dimll(ax,[xm,ym],[xm+vm1[0],ym+vm1[1]],r'$\mathrm{v_{M1}}$','left','center','b',1.2)
dm.dimll(ax,[xm,ym],[xm+vm2[0],ym+vm2[1]],r'$\mathrm{v_{M2}}$','left','center','r',1)


V1x = np.sqrt(vm1[0]**2+vm1[1]**2)*np.cos(alphaM1)
V1y = np.sqrt(vm1[0]**2+vm1[1]**2)*np.sin(alphaM1)
V2x = np.sqrt(vm2[0]**2+vm2[1]**2)*np.cos(alphaM2)
V2y = np.sqrt(vm2[0]**2+vm2[1]**2)*np.sin(alphaM2)

dm.dimll(ax,[xm,ym],[xm+V1x*np.cos(alpha),ym+V1x*np.sin(alpha)],r'$\mathrm{v_{b1}}=\mathrm{v_{b2}}$','left','center','b',1.2)
dm.dimll(ax,[xm,ym],[xm+V2x*np.cos(alpha),ym+V2x*np.sin(alpha)],'','center','bottom','r',1)

plt.plot([xm+V2x*np.cos(alpha),xm+vm2[0]], [ym+V2x*np.sin(alpha),ym+vm2[1]],'k--',lw=0.5)
## cotagem
dm.circlew(ax,0.15*r1, 150, 290, 0, a,r'$\omega_1$')
dm.circlew(ax,0.15*r1, -150, -290, 0, 0,r'$\omega_2$')
alpha1 = np.pi/13
dm.circlewx(ax,0.75*r2, 90+180*alpha/np.pi, 90-180*(alphaM2-alpha)/np.pi, 0, 0,r'$\alpha_2$')
dm.circlewx(ax,0.5*r2, 90+180*alpha/np.pi, 90, 0, 0,r'$\alpha$')

dm.circlewx(ax,0.75*r1, -90+180*alpha/np.pi, -90+180*(alpha-alphaM1)/np.pi, 0, a,r'$\alpha_1$')
dm.circlewx(ax,0.5*r1, -90+180*alpha/np.pi, -90, 0, a,r'$\alpha$')
dm.dimt(T1[0],T1[1],r'$\mathrm{T_1}$','left','bottom')
dm.dimt(T2[0],T2[1],r'$\mathrm{T_2}$','right','top')
dm.dimt(O1[0],O1[1],r'$\mathrm{O_1}$','left','bottom')
dm.dimt(O2[0],O2[1],r'$\mathrm{O_2}$','left','top')
dm.dimt(xm,ym,r'$\mathrm{M}$','right','center')
plt.savefig('pdf/vbM.pdf', bbox_inches = 'tight',
    pad_inches = 0, transparent=True)

###############################################################################
## VELOCIDADES ROLAMENTO ######################################################
###############################################################################
## PLOT ##
plt.figure(6)
fig = plt.gcf()
ax = fig.gca(aspect='equal')
plt.axis('off')
## reta de engrenamento
T1 = dm.cxy(rb1,0,a,alpha-np.pi/2)
T2 = dm.cxy(rb2,0,0,alpha+np.pi/2)


xl=20
dm.dimln(ax,[T1[0]+xl*np.cos(alpha),T1[1]+xl*np.sin(alpha)], T2,'','k-')
dm.dimln(ax,O2,T2,'','k-')
dm.dimln(ax,O1,T1,'','k-')
dm.dimln(ax,O1,O2,'','k-')
dm.dimln(ax,O2,[xm,ym],'','k-')
dm.dimln(ax,O1,[xm,ym],'','k-')
## Speed
dm.dimll(ax,[xm,ym],[xm+vm1[0],ym+vm1[1]],r'$\mathrm{v_{M1}}$','left','center','b',1.2)
dm.dimll(ax,[xm,ym],[xm+vm2[0],ym+vm2[1]],r'$\mathrm{v_{M2}}$','left','center','r',1)


V1x = np.sqrt(vm1[0]**2+vm1[1]**2)*np.cos(alphaM1)
V1y = np.sqrt(vm1[0]**2+vm1[1]**2)*np.sin(alphaM1)
V2x = np.sqrt(vm2[0]**2+vm2[1]**2)*np.cos(alphaM2)
V2y = np.sqrt(vm2[0]**2+vm2[1]**2)*np.sin(alphaM2)

dm.dimll(ax,[xm,ym],[xm+V2y*np.cos(alpha-np.pi/2),ym+V2y*np.sin(alpha-np.pi/2)],r'$\mathrm{v_{r2}}$','center','bottom','r',1)
dm.dimll(ax,[xm,ym],[xm+V1y*np.cos(alpha-np.pi/2),ym+V1y*np.sin(alpha-np.pi/2)],r'$\mathrm{v_{r1}}$','center','bottom','b',1.2)


plt.plot([xm+V2x*np.cos(alpha),xm+vm2[0]], [ym+V2x*np.sin(alpha),ym+vm2[1]],'k--',lw=0.5)
plt.plot([xm+V2y*np.cos(alpha-np.pi/2),xm+vm2[0]], [ym+V2y*np.sin(alpha-np.pi/2),ym+vm2[1]],'k--',lw=0.5)
## cotagem
dm.circlew(ax,0.15*r1, 150, 290, 0, a,r'$\omega_1$')
dm.circlew(ax,0.15*r1, -150, -290, 0, 0,r'$\omega_2$')
alpha1 = np.pi/13
dm.circlewx(ax,0.75*r2, 90+180*alpha/np.pi, 90-180*(alphaM2-alpha)/np.pi, 0, 0,r'$\alpha_2$')
dm.circlewx(ax,0.5*r2, 90+180*alpha/np.pi, 90, 0, 0,r'$\alpha$')

dm.circlewx(ax,0.75*r1, -90+180*alpha/np.pi, -90+180*(alpha-alphaM1)/np.pi, 0, a,r'$\alpha_1$')
dm.circlewx(ax,0.5*r1, -90+180*alpha/np.pi, -90, 0, a,r'$\alpha$')
dm.dimt(T1[0],T1[1],r'$\mathrm{T_1}$','left','bottom')
dm.dimt(T2[0],T2[1],r'$\mathrm{T_2}$','right','top')
dm.dimt(O1[0],O1[1],r'$\mathrm{O_1}$','left','bottom')
dm.dimt(O2[0],O2[1],r'$\mathrm{O_2}$','left','top')
dm.dimt(xm,ym,r'$\mathrm{M}$','right','center')
plt.savefig('pdf/vrM.pdf', bbox_inches = 'tight',
    pad_inches = 0, transparent=True)

###############################################################################
## VELOCIDADES TODAS ##########################################################
###############################################################################
## PLOT ##
plt.figure(7)
fig = plt.gcf()
ax = fig.gca(aspect='equal')
plt.axis('off')
## reta de engrenamento
T1 = dm.cxy(rb1,0,a,alpha-np.pi/2)
T2 = dm.cxy(rb2,0,0,alpha+np.pi/2)
dm.dimt(xm,ym,r'$\mathrm{M}$','right','center')

xl=20
dm.dimln(ax,[T1[0]+xl*np.cos(alpha),T1[1]+xl*np.sin(alpha)], T2,'','k-')
dm.dimln(ax,O2,T2,'','k-')
dm.dimln(ax,O1,T1,'','k-')
dm.dimln(ax,O1,O2,'','k-')
dm.dimln(ax,O2,[xm,ym],'','k-')
dm.dimln(ax,O1,[xm,ym],'','k-')
## Speed
dm.dimll(ax,[xm,ym],[xm+vm1[0],ym+vm1[1]],r'$\mathrm{v_{M1}}$','left','center','b',1.2)
dm.dimll(ax,[xm,ym],[xm+vm2[0],ym+vm2[1]],r'$\mathrm{v_{M2}}$','left','center','r',1)


V1x = np.sqrt(vm1[0]**2+vm1[1]**2)*np.cos(alphaM1)
V1y = np.sqrt(vm1[0]**2+vm1[1]**2)*np.sin(alphaM1)
V2x = np.sqrt(vm2[0]**2+vm2[1]**2)*np.cos(alphaM2)
V2y = np.sqrt(vm2[0]**2+vm2[1]**2)*np.sin(alphaM2)

dm.dimll(ax,[xm,ym],[xm+V1x*np.cos(alpha),ym+V1x*np.sin(alpha)],r'$\mathrm{v_{b1}}=\mathrm{v_{b2}}$','left','center','b',1.2)
dm.dimll(ax,[xm,ym],[xm+V2x*np.cos(alpha),ym+V2x*np.sin(alpha)],'','center','bottom','r',1)
dm.dimll(ax,[xm,ym],[xm+V2y*np.cos(alpha-np.pi/2),ym+V2y*np.sin(alpha-np.pi/2)],r'$\mathrm{v_{r2}}$','center','bottom','r',1)
dm.dimll(ax,[xm,ym],[xm+V1y*np.cos(alpha-np.pi/2),ym+V1y*np.sin(alpha-np.pi/2)],r'$\mathrm{v_{r1}}$','center','bottom','b',1.2)

plt.plot([xm+V2x*np.cos(alpha),xm+vm2[0]], [ym+V2x*np.sin(alpha),ym+vm2[1]],'k--',lw=0.5)
plt.plot([xm+V2y*np.cos(alpha-np.pi/2),xm+vm2[0]], [ym+V2y*np.sin(alpha-np.pi/2),ym+vm2[1]],'k--',lw=0.5)
## cotagem
dm.circlew(ax,0.15*r1, 150, 290, 0, a,r'$\omega_1$')
dm.circlew(ax,0.15*r1, -150, -290, 0, 0,r'$\omega_2$')
alpha1 = np.pi/13
dm.circlewx(ax,0.75*r2, 90+180*alpha/np.pi, 90-180*(alphaM2-alpha)/np.pi, 0, 0,r'$\alpha_2$')
dm.circlewx(ax,0.5*r2, 90+180*alpha/np.pi, 90, 0, 0,r'$\alpha$')

dm.circlewx(ax,0.75*r1, -90+180*alpha/np.pi, -90+180*(alpha-alphaM1)/np.pi, 0, a,r'$\alpha_1$')
dm.circlewx(ax,0.5*r1, -90+180*alpha/np.pi, -90, 0, a,r'$\alpha$')
dm.dimt(T1[0],T1[1],r'$\mathrm{T_1}$','left','bottom')
dm.dimt(T2[0],T2[1],r'$\mathrm{T_2}$','right','top')
dm.dimt(O1[0],O1[1],r'$\mathrm{O_1}$','left','bottom')
dm.dimt(O2[0],O2[1],r'$\mathrm{O_2}$','left','top')

plt.savefig('pdf/escorregamentoM.pdf', bbox_inches = 'tight',
    pad_inches = 0, transparent=True)


###############################################################################
## CIRCULOS OSCULADORES #######################################################
###############################################################################
## PLOT ##
plt.figure(8)
fig = plt.gcf()
ax = fig.gca(aspect='equal')
plt.axis('off')
## reta de engrenamento
T1 = dm.cxy(rb1,0,a,alpha-np.pi/2)
T2 = dm.cxy(rb2,0,0,alpha+np.pi/2)

xl=0
dm.dimln(ax,[T1[0]+xl*np.cos(alpha),T1[1]+xl*np.sin(alpha)], T2,'','k-')
dm.dimln(ax,O2,T2,'','k-')
dm.dimln(ax,O1,T1,'','k-')
dm.dimln(ax,O1,O2,'','k-')
## Speed
# dm.dimll(ax,[xm,ym],[xm+vm1[0],ym+vm1[1]],r'$\mathrm{v_{M1}}$','left','center','b',1.2)
# dm.dimll(ax,[xm,ym],[xm+vm2[0],ym+vm2[1]],r'$\mathrm{v_{M2}}$','left','center','r',1)

V1x = np.sqrt(vm1[0]**2+vm1[1]**2)*np.cos(alphaM1)
V1y = np.sqrt(vm1[0]**2+vm1[1]**2)*np.sin(alphaM1)
V2x = np.sqrt(vm2[0]**2+vm2[1]**2)*np.cos(alphaM2)
V2y = np.sqrt(vm2[0]**2+vm2[1]**2)*np.sin(alphaM2)

# dm.dimll(ax,[xm,ym],[xm+V1x*np.cos(alpha),ym+V1x*np.sin(alpha)],r'$\mathrm{v_{b1}}=\mathrm{v_{b2}}$','left','center','b',1.2)
# dm.dimll(ax,[xm,ym],[xm+V2x*np.cos(alpha),ym+V2x*np.sin(alpha)],'','center','bottom','r',1)
dm.dimll(ax,[xm,ym],[xm+V2y*np.cos(alpha-np.pi/2),ym+V2y*np.sin(alpha-np.pi/2)],r'$\mathrm{v_{r2}}$','center','bottom','r',1)
dm.dimll(ax,[xm,ym],[xm+V1y*np.cos(alpha-np.pi/2),ym+V1y*np.sin(alpha-np.pi/2)],r'$\mathrm{v_{r1}}$','center','bottom','b',1.2)

## cotagem
dm.circlew(ax,0.15*r1, 150, 290, 0, a,r'$\omega_1$')
dm.circlew(ax,0.15*r1, -150, -290, 0, 0,r'$\omega_2$')

dm.circlew(ax,0.15*r1, 150, 290, T1[0], T1[1],r'$\omega_1$')
dm.circlew(ax,0.15*r1, -150, -290, T2[0], T2[1],r'$\omega_2$')

T1T2 = a*np.sin(alpha)
T2M = np.sqrt(rM2**2-rb2**2)
T1M = T1T2 - T2M

dm.circle(T1M, 90, 290, T1[0], T1[1],'k')
dm.circle(T2M, -40, 100, T2[0], T2[1],'k')

alpha1 = np.pi/13
dm.dimt(T1[0],T1[1],r'$\mathrm{T_1}$','left','bottom')
dm.dimt(T2[0],T2[1],r'$\mathrm{T_2}$','right','top')
dm.dimt(O1[0],O1[1],r'$\mathrm{O_1}$','left','bottom')
dm.dimt(O2[0],O2[1],r'$\mathrm{O_2}$','left','top')
dm.dimt(xm,ym,r'$\mathrm{M}$','right','bottom')
# dm.dimt(0,r2,r'$\mathrm{I}$','right','bottom')
plt.savefig('pdf/osculadores.pdf', bbox_inches = 'tight',
    pad_inches = 0, transparent=True)

###############################################################################
## VELOCIDADES LINEARES #######################################################
###############################################################################
## ROTAÇÃO ##
rM2 = rl2 + m/2 
alphaM1, alphaM2, angR1, angR2 = cont(rM2)
vf1 = np.add(dm.matrix(angR1,vi1), (O1*np.ones((evsize+rasize,2))).T)
vf2 = np.add(dm.matrix(angR2,vi2), (O2*np.ones((evsize+rasize,2))).T)
## VELOCIDADE 
xm = rM2*np.sin(alphaM2-alpha)
ym = rM2*np.cos(alphaM2-alpha)
O2M = np.array([xm, ym, 0])
O1M = np.array([0, a, 0]) - O2M
w2 = np.array([0, 0, wm])
w1 = w2*z2/z1
vm1 = np.cross(w1,O1M)
vm2 = np.cross(w2,O2M) 
# ## PLOT ##
plt.figure(9)
fig = plt.gcf()
ax = fig.gca(aspect='equal')
plt.axis('off')
## reta de engrenamento
T1 = dm.cxy(rb1,0,a,alpha-np.pi/2)
T2 = dm.cxy(rb2,0,0,alpha+np.pi/2)

# raio de base
dm.circle(rb1, ang1s, ang1f, 0, a,'k--')
dm.circle(rb2, ang2s, ang2f, 0, 0,'k--')

xl=0
dm.dimln(ax,[T1[0]+xl*np.cos(alpha),T1[1]+xl*np.sin(alpha)], T2,'','k-')
dm.dimln(ax,O2,T2,'','k-')
dm.dimln(ax,O1,T1,'','k-')
dm.dimln(ax,O1,O2,'','k-')
dm.dimln(ax,O2,[xm,ym],'','k-')
dm.dimln(ax,O1,[xm,ym],'','k-')
## Speed

T1T2 = a*np.sin(alpha)
T2M = np.sqrt(rM2**2-rb2**2)
T1M = T1T2 - T2M

dm.circle(T1M, 90, 290, T1[0], T1[1],'k')
dm.circle(T2M, -40, 100, T2[0], T2[1],'k')

## cotagem
dm.circlew(ax,0.15*r1, 150, 290, 0, a,r'$\omega_1$')
dm.circlew(ax,0.15*r1, -150, -290, 0, 0,r'$\omega_2$')
alpha1 = np.pi/13
dm.circlewx(ax,0.75*r2, 90+180*alpha/np.pi, 90-180*(alphaM2-alpha)/np.pi, 0, 0,r'$\alpha_2$')
dm.circlewx(ax,0.5*r2, 90+180*alpha/np.pi, 90, 0, 0,r'$\alpha$')

dm.circlewx(ax,0.75*r1, -90+180*alpha/np.pi, -90+180*(alpha-alphaM1)/np.pi, 0, a,r'$\alpha_1$')
dm.circlewx(ax,0.5*r1, -90+180*alpha/np.pi, -90, 0, a,r'$\alpha$')
dm.dimt(T1[0],T1[1],r'$\mathrm{T_1}$','left','bottom')
dm.dimt(T2[0],T2[1],r'$\mathrm{T_2}$','right','top')
dm.dimt(O1[0],O1[1],r'$\mathrm{O_1}$','left','bottom')
dm.dimt(O2[0],O2[1],r'$\mathrm{O_2}$','left','top')

## Perfis
plt.plot(vf1[0],vf1[1],'b',lw=1.5)
plt.plot(vf2[0],vf2[1],'r',lw=1.5)
dm.dimt(xm,ym,r'$\mathrm{M}$','right','center')
plt.savefig('pdf/conceptual.pdf', bbox_inches = 'tight',
    pad_inches = 0, transparent=True)