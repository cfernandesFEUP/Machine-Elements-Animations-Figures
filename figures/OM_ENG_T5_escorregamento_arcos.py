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
alphaAA2 = np.arctan(T2A/rb2) 
deltaA1 = -0.1
deltaA2 = -0.015
def phiA(ang_inc):
    return np.linspace(0, ang_inc + dm.inv(ang_inc), evsize)

def phiB(ang0, ang_inc):
    return np.linspace(ang0 + dm.inv(ang0), ang_inc + dm.inv(ang_inc), evsize)
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

ang2s, ang2f = 60, 120
ang1s, ang1f = ang2s + 180, ang2f +180
###############################################################################
## ARCOS ######################################################################
###############################################################################
## ROTAÇÃO ##
rM2 = np.linspace(rA2,rl2,100)
vfA1 = np.zeros((vi1.shape[0],vi1.shape[1],len(rM2)))
vfA2 = np.zeros((vi2.shape[0],vi2.shape[1],len(rM2)))
alphaMA1, alphaMA2, angRA1, angRA2, xm, ym = [np.zeros((len(rM2))) for _ in range(6)]

for i in range(len(rM2)):
    alphaMA1[i], alphaMA2[i], angRA1[i], angRA2[i] = cont(rM2[i])
    vfA1[:,:,i] = np.add(dm.matrix(angRA1[i],vi1), (O1*np.ones((evsize+rasize,2))).T)
    vfA2[:,:,i] = np.add(dm.matrix(angRA2[i],vi2), (O2*np.ones((evsize+rasize,2))).T)
    xm[i] = rM2[i]*np.sin(alphaMA2[i]-alpha)
    ym[i] = rM2[i]*np.cos(alphaMA2[i]-alpha)

    
i = 0
xyiP1 = np.add(dm.matrix(angRA1[i],dm.ev(rb1,phiB(-alphaMA1[0],-alphaMA1[i]))), (O1*np.ones((len(rM2),2))).T)
xyiP2 = np.add(dm.matrix(angRA2[i],dm.ev(rb2,phiB(-alphaMA2[0],-alphaMA2[i]))), (O2*np.ones((len(rM2),2))).T)


## PLOT ##
plt.figure(1)
fig = plt.gcf()
ax = fig.gca(aspect='equal')
plt.axis('off')
ax.plot(vfA1[0,:,i], vfA1[1,:,i], 'b',lw=0.5)
ax.plot(vfA2[0,:,i], vfA2[1,:,i], 'r',lw=0.5)
ax.plot(xm[i], ym[i], 'k.',ms=5)
ax.plot(xyiP1[0], xyiP1[1], 'k',lw=1)
ax.plot(xyiP2[0], xyiP2[1], 'k',lw=1)
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
## cotagem
dm.dimt(T1[0],T1[1],r'$\mathrm{T_1}$','left','bottom')
dm.dimt(T2[0],T2[1],r'$\mathrm{T_2}$','right','top')
dm.dimt(0,r2,r'$\mathrm{I}$','right','bottom')

dm.dimtt(vfA1[0,-1,i], vfA1[1,-1,i],r'$\mathrm{p_1}$','center','top')
dm.dimtt(vfA2[0,-1,i], vfA2[1,-1,i],r'$\mathrm{p_2}$','center','bottom')
angA =  -np.pi/2-alphaA1+alpha
angB =  np.pi/2+alpha-alphaA2
dm.dimtarray(dm.cxy(ra1, O1[0], O1[1], angA),r'$\mathrm{A}$','right','top')
dm.dimtarray(dm.cxy(ra2, O2[0], O2[1], angB),r'$\mathrm{B}$','left','bottom')
plt.savefig('pdf/arcosA.pdf', bbox_inches = 'tight',
    pad_inches = 0, transparent=True)

i = -1
xyiP1 = np.add(dm.matrix(angRA1[i],dm.ev(rb1,phiB(-alphaMA1[0],-alphaMA1[i]))), (O1*np.ones((len(rM2),2))).T)
xyiP2 = np.add(dm.matrix(angRA2[i],dm.ev(rb2,phiB(-alphaMA2[0],-alphaMA2[i]))), (O2*np.ones((len(rM2),2))).T)


## PLOT ##
plt.figure(2)
fig = plt.gcf()
ax = fig.gca(aspect='equal')
plt.axis('off')
ax.plot(vfA1[0,:,i], vfA1[1,:,i], 'b',lw=0.5)
ax.plot(vfA2[0,:,i], vfA2[1,:,i], 'r',lw=0.5)
ax.plot(xm[i], ym[i], 'k.',ms=5)
ax.plot(xyiP1[0], xyiP1[1], 'k',lw=1)
ax.plot(xyiP2[0], xyiP2[1], 'k',lw=1)
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
## cotagem
dm.dimt(T1[0],T1[1],r'$\mathrm{T_1}$','left','bottom')
dm.dimt(T2[0],T2[1],r'$\mathrm{T_2}$','right','top')
dm.dimt(0,r2,r'$\mathrm{I}$','right','bottom')

dm.dimt(xyiP1[0,0],xyiP1[1,0],r'$\mathrm{a_1}$','left','bottom')
dm.dimt(xyiP2[0,0],xyiP2[1,0],r'$\mathrm{a_2}$','right','top')
dm.dimtt(vfA1[0,-1,i], vfA1[1,-1,i],r'$\mathrm{p_1}$','center','top')
dm.dimtt(vfA2[0,-1,i], vfA2[1,-1,i],r'$\mathrm{p_2}$','center','bottom')
angA =  -np.pi/2-alphaA1+alpha
angB =  np.pi/2+alpha-alphaA2
dm.dimtarray(dm.cxy(ra1, O1[0], O1[1], angA),r'$\mathrm{A}$','right','top')
dm.dimtarray(dm.cxy(ra2, O2[0], O2[1], angB),r'$\mathrm{B}$','left','bottom')
plt.savefig('pdf/arcosI.pdf', bbox_inches = 'tight',
    pad_inches = 0, transparent=True)

## ROTAÇÃO ##
rM2 = np.linspace(rl2,ra2,100)
vfA1 = np.zeros((vi1.shape[0],vi1.shape[1],len(rM2)))
vfA2 = np.zeros((vi2.shape[0],vi2.shape[1],len(rM2)))
alphaMA1, alphaMA2, angRA1, angRA2, xm, ym = [np.zeros((len(rM2))) for _ in range(6)]

for i in range(len(rM2)):
    alphaMA1[i], alphaMA2[i], angRA1[i], angRA2[i] = cont(rM2[i])
    vfA1[:,:,i] = np.add(dm.matrix(angRA1[i],vi1), (O1*np.ones((evsize+rasize,2))).T)
    vfA2[:,:,i] = np.add(dm.matrix(angRA2[i],vi2), (O2*np.ones((evsize+rasize,2))).T)
    xm[i] = rM2[i]*np.sin(alphaMA2[i]-alpha)
    ym[i] = rM2[i]*np.cos(alphaMA2[i]-alpha)

    
i = -1
xyiP1 = np.add(dm.matrix(angRA1[i],dm.ev(rb1,phiB(-alphaMA1[0],-alphaMA1[i]))), (O1*np.ones((len(rM2),2))).T)
xyiP2 = np.add(dm.matrix(angRA2[i],dm.ev(rb2,phiB(-alphaMA2[0],-alphaMA2[i]))), (O2*np.ones((len(rM2),2))).T)


## PLOT ##
plt.figure(3)
fig = plt.gcf()
ax = fig.gca(aspect='equal')
plt.axis('off')
ax.plot(vfA1[0,:,i], vfA1[1,:,i], 'b',lw=0.5)
ax.plot(vfA2[0,:,i], vfA2[1,:,i], 'r',lw=0.5)
ax.plot(xm[i], ym[i], 'k.',ms=5)
ax.plot(xyiP1[0], xyiP1[1], 'k',lw=1)
ax.plot(xyiP2[0], xyiP2[1], 'k',lw=1)
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
## cotagem
dm.dimt(T1[0],T1[1],r'$\mathrm{T_1}$','left','bottom')
dm.dimt(T2[0],T2[1],r'$\mathrm{T_2}$','right','top')
dm.dimt(0,r2,r'$\mathrm{I}$','right','bottom')
dm.dimt(xyiP1[0,0],xyiP1[1,0],r'$\mathrm{i_1}$','left','bottom')
dm.dimt(xyiP2[0,0],xyiP2[1,0],r'$\mathrm{i_2}$','right','top')

dm.dimtt(vfA1[0,-1,i], vfA1[1,-1,i],r'$\mathrm{p_1}$','center','top')
dm.dimtt(vfA2[0,-1,i], vfA2[1,-1,i],r'$\mathrm{p_2}$','center','bottom')
angA =  -np.pi/2-alphaA1+alpha
angB =  np.pi/2+alpha-alphaA2
dm.dimtarray(dm.cxy(ra1, O1[0], O1[1], angA),r'$\mathrm{A}$','right','top')
dm.dimtarray(dm.cxy(ra2, O2[0], O2[1], angB),r'$\mathrm{B}$','left','bottom')
plt.savefig('pdf/arcosB.pdf', bbox_inches = 'tight',
    pad_inches = 0, transparent=True)


## ROTAÇÃO ##
rM2 = np.linspace(rA2,.988*rl2,100)
vfA1 = np.zeros((vi1.shape[0],vi1.shape[1],len(rM2)))
vfA2 = np.zeros((vi2.shape[0],vi2.shape[1],len(rM2)))
alphaMA1, alphaMA2, angRA1, angRA2, xm, ym = [np.zeros((len(rM2))) for _ in range(6)]

for i in range(len(rM2)):
    alphaMA1[i], alphaMA2[i], angRA1[i], angRA2[i] = cont(rM2[i])
    vfA1[:,:,i] = np.add(dm.matrix(angRA1[i],vi1), (O1*np.ones((evsize+rasize,2))).T)
    vfA2[:,:,i] = np.add(dm.matrix(angRA2[i],vi2), (O2*np.ones((evsize+rasize,2))).T)
    xm[i] = rM2[i]*np.sin(alphaMA2[i]-alpha)
    ym[i] = rM2[i]*np.cos(alphaMA2[i]-alpha)

    
i = -1
xyiP1 = np.add(dm.matrix(angRA1[i],dm.ev(rb1,phiB(-alphaMA1[0],-alphaMA1[i]))), (O1*np.ones((len(rM2),2))).T)
xyiP2 = np.add(dm.matrix(angRA2[i],dm.ev(rb2,phiB(-alphaMA2[0],-alphaMA2[i]))), (O2*np.ones((len(rM2),2))).T)


## PLOT ##
plt.figure(4)
fig = plt.gcf()
ax = fig.gca(aspect='equal')
plt.axis('off')
ax.plot(vfA1[0,:,i], vfA1[1,:,i], 'b',lw=0.5)
ax.plot(vfA2[0,:,i], vfA2[1,:,i], 'r',lw=0.5)
ax.plot(xm[i], ym[i], 'k.',ms=5)
ax.plot(xyiP1[0], xyiP1[1], 'k',lw=1)
ax.plot(xyiP2[0], xyiP2[1], 'k',lw=1)
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
## cotagem
dm.dimt(T1[0],T1[1],r'$\mathrm{T_1}$','left','bottom')
dm.dimt(T2[0],T2[1],r'$\mathrm{T_2}$','right','top')
# dm.dimt(0,r2,r'$\mathrm{I}$','right','bottom')
dm.dimt(xyiP1[0,0],xyiP1[1,0],r'$\mathrm{a}$','left','bottom')
dm.dimt(xyiP1[0,-1],xyiP1[1,-1],r'$\mathrm{b}$','left','bottom')
dm.dimt(xyiP2[0,0],xyiP2[1,0],r'$\mathrm{c}$','right','top')
dm.dimt(xyiP2[0,-1],xyiP2[1,-1],r'$\mathrm{d}$','right','top')
# dm.dimtt(vfA1[0,-1,i], vfA1[1,-1,i],r'$\mathrm{p_1}$','center','top')
# dm.dimtt(vfA2[0,-1,i], vfA2[1,-1,i],r'$\mathrm{p_2}$','center','bottom')
angA =  -np.pi/2-alphaA1+alpha
angB =  np.pi/2+alpha-alphaA2
dm.dimtarray(dm.cxy(ra1, O1[0], O1[1], angA),r'$\mathrm{A}$','right','top')
dm.dimtarray(dm.cxy(ra2, O2[0], O2[1], angB),r'$\mathrm{B}$','left','bottom')
plt.savefig('pdf/arcosgs.pdf', bbox_inches = 'tight',
    pad_inches = 0, transparent=True)




## PLOT ##
plt.figure(5)
fig = plt.gcf()
ax = fig.gca(aspect='equal')
plt.axis('off')
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
## cotagem
dm.dimt(T1[0],T1[1],r'$\mathrm{T_1}$','left','bottom')
dm.dimt(T2[0],T2[1],r'$\mathrm{T_2}$','right','top')
dm.dimt(0,r2,r'$\mathrm{I}$','right','bottom')
angA =  -np.pi/2-alphaA1+alpha
angB =  np.pi/2+alpha-alphaA2
dm.dimtarray(dm.cxy(ra1, O1[0], O1[1], angA),r'$\mathrm{A}$','right','top')
dm.dimtarray(dm.cxy(ra2, O2[0], O2[1], angB),r'$\mathrm{B}$','left','bottom')
plt.savefig('pdf/gs_AB.pdf', bbox_inches = 'tight',
    pad_inches = 0, transparent=True)