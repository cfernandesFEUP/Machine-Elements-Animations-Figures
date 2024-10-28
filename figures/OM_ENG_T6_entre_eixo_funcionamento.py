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
z1 = 10.
z2 = 15.
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
def cont(rM2,alpha,T1T2):
    alphaM2 = np.arccos(rb2/rM2)
    angR2 = np.tan(alphaM2)-alpha
    T2M = np.sqrt(rM2**2-rb2**2)
    T1M = T1T2 - T2M
    rM1 = np.sqrt(T1M**2+rb1**2)
    alphaM1 = np.arccos(rb1/rM1)
    angR1 = np.pi + np.tan(alphaM1) - alpha
    return alphaM1, alphaM2, angR1, angR2
###############################################################################
## ENTRE-EIXO FUNCIONAMENTO ###################################################
###############################################################################
ang2s, ang2f = 60, 120
ang1s, ang1f = ang2s + 170, ang2f + 190
## ROTAÇÃO ##
rMa2 = r2
alphaM1, alphaM2, angR1, angR2 = cont(rMa2,alpha,T1T2)
vf1 = np.add(dm.matrix(angR1,vi1), (O1*np.ones((evsize+rasize,2))).T)
vf2 = np.add(dm.matrix(angR2,vi2), (O2*np.ones((evsize+rasize,2))).T)
al = a + 2
O1 = np.array([0,al])
O2 = np.array([0,0])
alpha = np.arccos(a*np.cos(alpha)/al)
rl1 = rb1/np.cos(alpha)
rl2 = rb2/np.cos(alpha)
T1T2 = al*np.sin(alpha)
## ROTAÇÃO ##
rMa2 = rl2
alphaM1, alphaM2, angR1, angR2 = cont(rMa2,alpha,T1T2)
vf1 = np.add(dm.matrix(angR1,vi1), (O1*np.ones((evsize+rasize,2))).T)
vf2 = np.add(dm.matrix(angR2,vi2), (O2*np.ones((evsize+rasize,2))).T)
## PLOT ##
plt.figure(1)
fig = plt.gcf()
ax = fig.gca(aspect='equal')
plt.axis('off')
# raio primitivo
angp1 = 2*np.pi/z1
angp2 = 2*np.pi/z2
dm.circle(rl1, ang1s, ang1f, 0, al,'k')
dm.circle(rl2, ang2s, ang2f, 0, 0,'k')
dm.dimln(ax,O1,O2,'','k-')
dm.circlewxb(ax,rl2, 90+180*angp2/np.pi, 90, 0, 0,r"$p'$")
dm.circlewx(ax,rl1, -90-180*angp1/np.pi, -90, 0, al,r"$p'$")
## cotagem
dm.dimt(O1[0],O1[1],r"$\mathrm{O_1}$",'left','bottom')
dm.dimt(O2[0],O2[1],r'$\mathrm{O_2}$','left','top')
dm.dimlr(ax,O1,dm.cxy(rl1, O1[0], O1[1], -np.pi/2+0.6*alphaA1),r"$\mathrm{r'_{1}}$")
dm.dimlr(ax,O2,dm.cxy(rl2, O2[0], O2[1], 1.8*alphaA2),r"$\mathrm{r'_{2}}$")

plt.savefig('pdf/passo_funcionamento.pdf', bbox_inches = 'tight',
    pad_inches = 0, transparent=True)

plt.figure(2)
fig = plt.gcf()
ax = fig.gca(aspect='equal')
plt.axis('off')
# raio primitivo
dm.circle(rl1, ang1s, ang1f, 0, al,'k')
dm.circle(rl2, ang2s, ang2f, 0, 0,'k')
dm.dimln(ax,O1,O2,'','k-')
angs1 = .6*angp1
angi1 = .4*angp1
angs2 = .4*angp2
angi2 = .6*angp2
dm.circlewxb(ax,rl2, 90,90+180*angi2/np.pi, 0, 0,r"$i'_2$")
dm.circlewxb(ax,rl2, 90+180*angi2/np.pi, 90+180*angp2/np.pi, 0, 0,r"$s'_2$")
dm.circlewx(ax,rl1, -90-180*angs1/np.pi, -90, 0, al,r"$s'_1$")
dm.circlewx(ax,rl1, -90-180*angs1/np.pi, -90-180*angp1/np.pi, 0, al,r"$i'_1$")

## cotagem
dm.dimt(O1[0],O1[1],r"$\mathrm{O_1}$",'left','bottom')
dm.dimt(O2[0],O2[1],r'$\mathrm{O_2}$','left','top')
dm.dimlr(ax,O1,dm.cxy(rl1, O1[0], O1[1], -np.pi/2+0.6*alphaA1),r"$\mathrm{r'_{1}}$")
dm.dimlr(ax,O2,dm.cxy(rl2, O2[0], O2[1], 1.8*alphaA2),r"$\mathrm{r'_{2}}$")

plt.savefig('pdf/s_i_funcionamento.pdf', bbox_inches = 'tight',
    pad_inches = 0, transparent=True)