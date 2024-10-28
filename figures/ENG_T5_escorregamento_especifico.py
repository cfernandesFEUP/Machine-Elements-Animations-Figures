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
u = z2/z1
alpha = np.pi/9
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
## ENGRENAMENTO COM ra1 #######################################################
###############################################################################
rMm2 = np.sqrt(T1T2**2+rb2**2)
rM2 = np.linspace(rb2,rMm2,450)
T1M, T2M, gs1, gs2 = [np.zeros(len(rM2)) for _ in range(4)]
## PLOT ## 
for i in range(len(rM2)):
    alphaM1, alphaM2, angR1, angR2 = cont(rM2[i])
    ## VELOCIDADE 
    wm = 0.5
    T1T2 = a*np.sin(alpha)
    T2M[i] = np.sqrt(rM2[i]**2-rb2**2)
    T1M[i] = T1T2 - T2M[i]
    w2 = wm
    w1 = w2*z2/z1
    gs1[i] = (T1M[i]*w1-T2M[i]*w2)/(T1M[i]*w1)
    gs2[i] = -(T1M[i]*w1-T2M[i]*w2)/(T2M[i]*w2)
# ## PLOT ##
plt.figure(5)
fig = plt.gcf()
ax = fig.gca(aspect='equal')
plt.axis('off')
## reta de engrenamento
T1 = dm.cxy(rb1,0,a,alpha-np.pi/2)
T2 = dm.cxy(rb2,0,0,alpha+np.pi/2)
gsi1 = np.array([T2M[:-17],gs1[:-17]])
gsf1 = dm.matrix(-alpha,gsi1)
gsi2 = np.array([T2M,gs2])
gsf2 = dm.matrix(-alpha,gsi2)

gs1B = np.array([T2B,(T1B*w1-T2B*w2)/(T1B*w1)])
gs2A = np.array([T2A,(T2A*w2-T1A*w1)/(T2A*w2)])
gsf1B = dm.matrix(-alpha,gs1B)
gsf2A = dm.matrix(-alpha,gs2A)

plt.plot(gsf1[0]+T2[0],gsf1[1]+T2[1],'b')
plt.plot(gsf2[0]+T2[0],gsf2[1]+T2[1],'r')
angA =  -np.pi/2-alphaA1+alpha
angB =  np.pi/2+alpha-alphaA2
dm.dimtarray(dm.cxy(ra1, O1[0], O1[1], angA),r'$\mathrm{A}$','center','bottom')
dm.dimtarray(dm.cxy(ra2, O2[0], O2[1], angB),r'$\mathrm{B}$','center','bottom')
##
dm.dimln(ax,dm.cxy(ra1, O1[0], O1[1], angA), [gsf2A[0]+T2[0],gsf2A[1]+T2[1]],'','k-')
dm.dimln(ax,dm.cxy(ra2, O2[0], O2[1], angB), [gsf1B[0]+T2[0],gsf1B[1]+T2[1]],'','k-')
dm.dimln(ax,T1, T2,'','k-')
dm.dimln(ax,O2,T2,'','k-')
dm.dimln(ax,O1,T1,'','k-')
dm.dimln(ax,O1,O2,'','k-')
## cotagem
dm.circlew(ax,0.15*r1, 150, 290, 0, a,r'$\omega_1$')
dm.circlew(ax,0.15*r1, -150, -290, 0, 0,r'$\omega_2$')

alpha1 = np.pi/13
dm.dimt(gsf2A[0]+T2[0],gsf2A[1]+T2[1],r'$\mathrm{g_{s2A}}$','center','top')
dm.dimt(gsf1B[0]+T2[0],gsf1B[1]+T2[1],r'$\mathrm{g_{s1B}}$','center','top')
dm.dimt(T1[0],T1[1],r'$\mathrm{T_1}$','left','bottom')
dm.dimt(T2[0],T2[1],r'$\mathrm{T_2}$','right','top')
dm.dimt(O1[0],O1[1],r'$\mathrm{O_1}$','left','bottom')
dm.dimt(O2[0],O2[1],r'$\mathrm{O_2}$','left','top')
dm.dimt(0,r2,r'$\mathrm{I}$','right','bottom')
plt.savefig('pdf/eespecificoA.pdf', bbox_inches = 'tight',
    pad_inches = 0, transparent=True)


# ## PLOT ##
plt.figure(6)
fig = plt.gcf()
ax = fig.gca(aspect='equal')
plt.axis('off')
## reta de engrenamento
T1 = dm.cxy(rb1,0,a,alpha-np.pi/2)
T2 = dm.cxy(rb2,0,0,alpha+np.pi/2)
# gsi1 =np.array([T2M,gs1])
# gsf1 = dm.matrix(-alpha,gsi1)
# gsi2 =np.array([T2M,gs2])
# gsf2 = dm.matrix(-alpha,gsi2)

gs1B = np.array([T2B,(T1B*w1-T2B*w2)/(T1B*w1)])
gs2A = np.array([T2A,(T2A*w2-T1A*w1)/(T2A*w2)])
gsf1B = dm.matrix(-alpha,gs1B)
gsf2A = dm.matrix(-alpha,gs2A)

from scipy import optimize
def corr(x):
    ra1 = r1 + m*(1+x)
    ra2 = r2 + m*(1-x)
    T1T2 = a*np.sin(alpha)
    T1B = np.sqrt(ra1**2-rb1**2)
    T2B = T1T2 - T1B
    T2A = np.sqrt(ra2**2-rb2**2)
    T1A = T1T2 - T2A
    gs1B = abs(1 - (1/u)*T2B/T1B)
    gs2A = abs(u*T1A/T2A - 1)
    return gs1B - gs2A
sol  =  optimize.brentq(corr,  0.0,  0.7)


## CALCULAR GEOMETRIA ##
rl1, rl2 = r1 + m*sol, r2 - m*sol
ral1, ral2 = r1 + m*(1+sol), r2 + m*(1-sol)
a = rl1+rl2
T2Bl = np.sqrt(ral2**2-rb2**2)
T1Al = np.sqrt(ral1**2-rb1**2)
T2Al = T1T2 - T1Al
T1Bl = T1T2 - T2Bl

alphaA1l = np.arctan(T1Al/rb1)
alphaA2l = np.arctan(T2Bl/rb2) 
angAl =  -np.pi/2-alphaA1l+alpha
angBl =  np.pi/2+alpha-alphaA2l

gs1Bl = np.array([T2Bl,(T1Bl*w1-T2Bl*w2)/(T1Bl*w1)])
gs2Al = np.array([T2Al,(T2Al*w2-T1Al*w1)/(T2Al*w2)])
gsf1Bl = dm.matrix(-alpha,gs1Bl)
gsf2Al = dm.matrix(-alpha,gs2Al)

plt.plot(gsf1[0]+T2[0],gsf1[1]+T2[1],'b')
plt.plot(gsf2[0]+T2[0],gsf2[1]+T2[1],'r')
angA =  -np.pi/2-alphaA1+alpha
angB =  np.pi/2+alpha-alphaA2
dm.dimtarray(dm.cxy(ra1, O1[0], O1[1], angA),'','center','bottom')
dm.dimtarray(dm.cxy(ra2, O2[0], O2[1], angB),'','center','bottom')

dm.dimtarray(dm.cxy(ral1, O1[0], O1[1], angAl),r"$\mathrm{A'}$",'center','bottom')
dm.dimtarray(dm.cxy(ral2, O2[0], O2[1], angBl),r"$\mathrm{B'}$",'center','bottom')
##
dm.dimln(ax,dm.cxy(ra1, O1[0], O1[1], angA), [gsf2A[0]+T2[0],gsf2A[1]+T2[1]],'','k-')
dm.dimln(ax,dm.cxy(ra2, O2[0], O2[1], angB), [gsf1B[0]+T2[0],gsf1B[1]+T2[1]],'','k-')

dm.dimln(ax,dm.cxy(ral1, O1[0], O1[1], angAl), [gsf2Al[0]+T2[0],gsf2Al[1]+T2[1]],'','k-')
dm.dimln(ax,dm.cxy(ral2, O2[0], O2[1], angBl), [gsf1Bl[0]+T2[0],gsf1Bl[1]+T2[1]],'','k-')

dm.dimln(ax,[gsf2Al[0]+T2[0],gsf2Al[1]+T2[1]],[gsf1Bl[0]+T2[0],gsf1Bl[1]+T2[1]],'','k--')

dm.dimln(ax,T1, T2,'','k-')
dm.dimln(ax,O2,T2,'','k-')
dm.dimln(ax,O1,T1,'','k-')
dm.dimln(ax,O1,O2,'','k-')
## cotagem

dm.dimlr(ax,O1,[O1[0]-ra1*np.sin(np.pi/9),O1[1]-ra1*np.cos(np.pi/9)],r'$\mathrm{r_{a1}}$')
dm.dimlr(ax,O1,[O1[0]-ral1*np.sin(np.pi/5),O1[1]-ral1*np.cos(np.pi/5)],r"$\mathrm{r_{a1}'}$")
dm.dimlr(ax,O2,[O2[0]+ra2*np.sin(np.pi/9),O2[1]+ra2*np.cos(np.pi/9)],r'$\mathrm{r_{a2}}$')
dm.dimlr(ax,O2,[O2[0]+ral2*np.sin(np.pi/5),O2[1]+ral2*np.cos(np.pi/5)],r"$\mathrm{r_{a2}'}$")

dm.circlet(ra1, -90, -130, 0, a,'k')
dm.circlet(ra2, 50, 90, 0, 0,'k')

dm.circlet(ral1, -90, -130, 0, a,'k--')
dm.circlet(ral2, 50, 90, 0, 0,'k--')

dm.circlew(ax,0.15*r1, 150, 290, 0, a,r'$\omega_1$')
dm.circlew(ax,0.15*r1, -150, -290, 0, 0,r'$\omega_2$')

dm.dimt(gsf2A[0]+T2[0],gsf2A[1]+T2[1],'','left','top')
dm.dimt(gsf1B[0]+T2[0],gsf1B[1]+T2[1],'','left','top')

dm.dimt(gsf2Al[0]+T2[0],gsf2Al[1]+T2[1],r"$\mathrm{g_{s2A'}}$",'center','top')
dm.dimt(gsf1Bl[0]+T2[0],gsf1Bl[1]+T2[1],r"$\mathrm{g_{s1B'}}$",'left','top')
dm.dimt(T1[0],T1[1],r'$\mathrm{T_1}$','left','bottom')
dm.dimt(T2[0],T2[1],r'$\mathrm{T_2}$','right','top')
dm.dimt(O1[0],O1[1],r'$\mathrm{O_1}$','left','bottom')
dm.dimt(O2[0],O2[1],r'$\mathrm{O_2}$','left','top')
dm.dimt(0,r2,r'$\mathrm{I}$','right','bottom')
plt.savefig('pdf/eespecificoB.pdf', bbox_inches = 'tight',
    pad_inches = 0, transparent=True)