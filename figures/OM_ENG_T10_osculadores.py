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
## CIRCULOS OSCULADORES #######################################################
###############################################################################
rM2 = rl2
alphaM1, alphaM2, angR1, angR2 = cont(rM2)
vf1 = np.add(dm.matrix(angR1,vi1), (O1*np.ones((evsize+rasize,2))).T)
vf2 = np.add(dm.matrix(angR2,vi2), (O2*np.ones((evsize+rasize,2))).T)
## PLOT ##
plt.figure(1)
fig = plt.gcf()
ax = fig.gca(aspect='equal')
plt.axis('off')
## reta de engrenamento
T1 = dm.cxy(rb1,0,a,alpha-np.pi/2)
T2 = dm.cxy(rb2,0,0,alpha+np.pi/2)

dm.dimln(ax,T1, T2,'','k-')
dm.dimln(ax,O2,T2,'','k-')
dm.dimln(ax,O1,T1,'','k-')
dm.dimln(ax,O1,O2,'','k-')




## cotagem
dm.circlew(ax,0.15*r1, 150, 290, 0, a,r'$\omega_1$')
dm.circlew(ax,0.15*r1, -150, -290, 0, 0,r'$\omega_2$')


dm.circlewx(ax,0.5*r2, 90+180*alpha/np.pi, 90, 0, 0,r'$\alpha^\prime$')
dm.circlewx(ax,0.5*r1, -90+180*alpha/np.pi, -90, 0, a,r'$\alpha^\prime$')

dm.dimlr(ax,T1, [T1[0]+r1*np.sin(alpha)*np.cos(np.pi/3), T1[1]+r1*np.sin(alpha)*np.sin(np.pi/3)],r'$\rho_1$')
dm.dimlr(ax,T2, [T2[0]+r2*np.sin(alpha)*np.cos(np.pi/3), T2[1]+r2*np.sin(alpha)*np.sin(np.pi/3)],r'$\rho_2$')

T1T2 = a*np.sin(alpha)
T2M = np.sqrt(rM2**2-rb2**2)
T1M = T1T2 - T2M

dm.circle(T1M, 60, 290, T1[0], T1[1],'k')
dm.circle(T2M, -40, 100, T2[0], T2[1],'k')

dm.dimt(T1[0],T1[1],r'$\mathrm{T_1}$','left','top')
dm.dimt(T2[0],T2[1],r'$\mathrm{T_2}$','right','top')
dm.dimt(O1[0],O1[1],r'$\mathrm{O_1}$','left','bottom')
dm.dimt(O2[0],O2[1],r'$\mathrm{O_2}$','left','top')
plt.plot(vf1[0],vf1[1],'b')
plt.plot(vf2[0],vf2[1],'r')
dm.dimt(0,rl2,r'$\mathrm{I}$','left','bottom')
plt.savefig('pdf/osculadoresHertz.pdf', bbox_inches = 'tight',
    pad_inches = 0, transparent=True)