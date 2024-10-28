## IMPORTAR LIBRARIAS ##
import matplotlib.pyplot as plt
import numpy as np
import functions.dimensioning as dm
from matplotlib import rcParams
# rcParams['font.family'] = 'sans-serif'
# rcParams['font.sans-serif'] = ['Fira Sans']
rcParams['font.size'] = 10
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
AB = T1T2-T2A-T1B
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
angA =  -np.pi/2-alphaA1+alpha
angB =  np.pi/2+alpha-alphaA2
###############################################################################
## VELOCIDADES ################################################################
###############################################################################
from matplotlib.backends.backend_pdf import PdfPages
# pp = PdfPages('pdf/anim/escorregamento.pdf')
pp = PdfPages('pdf/anim/vescorregamento.pdf')
rM2 = np.linspace(rA2,ra2,100)
j=0
xTM, vg = [], []
## PLOT ## 
for i in rM2:
    plt.figure()
    fig = plt.gcf()
    ax = fig.gca(aspect='equal')
    plt.axis('off')
    ## reta de engrenamento
    T1 = dm.cxy(rb1,0,a,alpha-np.pi/2)
    T2 = dm.cxy(rb2,0,0,alpha+np.pi/2)
    xl=20
    dm.dimln(ax,T1, T2,'','k-')
    dm.dimln(ax,O2,T2,'','k-')
    dm.dimln(ax,O1,T1,'','k-')
    dm.dimln(ax,O1,O2,'','k-')
    alphaM1, alphaM2, angR1, angR2 = cont(i)
    ## VELOCIDADE 
    wm = -0.75
    xm = i*np.sin(alphaM2-alpha)
    ym = i*np.cos(alphaM2-alpha)
    O2M = np.array([xm, ym, 0])
    O1M = np.array([0, a, 0]) - O2M
    w2 = np.array([0, 0, wm])
    w1 = w2*z2/z1
    vm1 = np.cross(w1,O1M)
    vm2 = np.cross(w2,O2M)
    dm.dimln(ax,O2,[xm,ym],'','k-')
    dm.dimln(ax,O1,[xm,ym],'','k-')
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
    vf1 = np.add(dm.matrix(angR1,vi1), (O1*np.ones((evsize+rasize,2))).T)
    vf2 = np.add(dm.matrix(angR2,vi2), (O2*np.ones((evsize+rasize,2))).T)
    plt.plot(vf1[0],vf1[1],'b')
    plt.plot(vf2[0],vf2[1],'r')
    ## cotagem
    dm.circlew(ax,0.15*r1, 150, 290, 0, a,r'$\omega_1$')
    dm.circlew(ax,0.15*r1, -150, -290, 0, 0,r'$\omega_2$')
    alpha1 = np.pi/13
    # dm.circlewx(ax,0.75*r2, 90+180*alpha/np.pi, 90-180*(alphaM2-alpha)/np.pi, 0, 0,r'$\alpha_2$')
    # dm.circlewx(ax,0.5*r2, 90+180*alpha/np.pi, 90, 0, 0,r'$\alpha$')
    
    # dm.circlewx(ax,0.75*r1, -90+180*alpha/np.pi, -90+180*(alpha-alphaM1)/np.pi, 0, a,r'$\alpha_1$')
    # dm.circlewx(ax,0.5*r1, -90+180*alpha/np.pi, -90, 0, a,r'$\alpha$')
    
    dm.dimtarray(dm.cxy(ra1, O1[0], O1[1], angA),r'$\mathrm{A}$','right','top')
    dm.dimtarray(dm.cxy(ra2, O2[0], O2[1], angB),r'$\mathrm{B}$','left','bottom')
    dm.dimt(T1[0],T1[1],r'$\mathrm{T_1}$','left','bottom')
    dm.dimt(T2[0],T2[1],r'$\mathrm{T_2}$','right','top')
    dm.dimt(O1[0],O1[1],r'$\mathrm{O_1}$','left','bottom')
    dm.dimt(O2[0],O2[1],r'$\mathrm{O_2}$','left','top')
    dm.dimttd(xm,ym)
    # pp.savefig(fig, bbox_inches = 'tight', transparent=True)
    plt.figure()
    fig = plt.gcf()
    # ax = fig.gca(aspect='equal')
    # plt.axis('off')
    vgx = V1y*np.cos(alpha-np.pi/2)-V2y*np.cos(alpha-np.pi/2)
    vgy = V1y*np.sin(alpha-np.pi/2)-V1y*np.sin(alpha-np.pi/2)
    vg.append(np.sqrt(vgx**2+vgy**2))
    xTM.append(np.sqrt(i**2-rb2**2)-T2A)
    plt.plot(xTM,vg,'k')
    plt.xlabel(r'$\overline{AM}$ / mm')
    plt.ylabel(r'$|v_g|$ / ms$^{-1}$')
    plt.ylim([0,vg[0]])
    plt.xlim([0,AB])
    pp.savefig(fig, bbox_inches = 'tight', transparent=True)
    j += 1
pp.close()