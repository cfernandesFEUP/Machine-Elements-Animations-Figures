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
b = 15.
alpha = np.pi/6
beta = np.pi/6
x1 = 0.
x2 = 0.
hfP = 1.
haP = 1.25
roh = 0.38
## CALCULAR GEOMETRIA ##
alphat = np.arctan(np.tan(alpha)/np.cos(beta))
r1, r2 = m*z1/(2*np.cos(beta)), m*z2/(2*np.cos(beta))
rl1, rl2 = r1 + m*x1, r2 + m*x2
ra1, ra2 = r1 + m*(1+x1), r2 + m*(1+x2)
rb1, rb2 = r1*np.cos(alphat), r2*np.cos(alphat)
a = rl1+rl2
pn = np.pi*m
pbn = pn*np.cos(alpha)
pt = pn/np.cos(beta)
pbt = pt*np.cos(alphat)
betab = np.arctan(np.tan(beta)*np.cos(alphat))
## POSIÇÃO DOS CENTROS ##
O1 = np.array([0,a])
O2 = np.array([0,0])
###############################################################################
## RAZÃO DE CONDUÇÃO  #########################################################
###############################################################################
# ang2s, ang2f = 30, 150
# ang1s, ang1f = ang2s + 180, ang2f +180
ang2s, ang2f = 50, 130
ang1s, ang1f = ang2s + 180, ang2f +180
## PLOT ##
plt.figure(1)
fig = plt.gcf()
ax = fig.gca(aspect='equal')
plt.axis('off')
# raio de base
dm.circle(rb1, ang1s, ang1f, 0, a,'k--')
dm.circle(rb2, ang2s, ang2f, 0, 0,'k--')
## raio primitivo
dm.circle(r1, ang1s, ang1f, 0, a,'k')
dm.circle(r2, ang2s, ang2f, 0, 0,'k')
## reta de engrenamento
T1 = dm.cxy(rb1,0,a,alphat-np.pi/2)
T2 = dm.cxy(rb2,0,0,alphat+np.pi/2)
dm.dimln(ax,T1,T2,'','k-')
dm.dimln(ax,O2,T2,'','k-')
dm.dimln(ax,O1,T1,'','k-')
dm.dimln(ax,O1,O2,'','k-')
## cotagem
dm.circlewx(ax,0.5*r2, 90+180*alphat/np.pi, 90, 0, 0,r'$\alpha_t$')
dm.circlewx(ax,0.5*r1, -90+180*alphat/np.pi, -90, 0, a,r'$\alpha_t$')
dm.dimt(T1[0],T1[1],r'$\mathrm{T_1}$','left','bottom')
dm.dimt(T2[0],T2[1],r'$\mathrm{T_2}$','right','top')
dm.dimt(O1[0],O1[1],r'$\mathrm{O_1}$','left','bottom')
dm.dimt(O2[0],O2[1],r'$\mathrm{O_2}$','left','top')
dm.dimt(0,r2,r'$\mathrm{I}$','left','bottom')
dm.dimlr(ax,O1,dm.cxy(r1, O1[0], O1[1], -np.pi/2-alphat),r'$\mathrm{r_{1}}$')
dm.dimlr(ax,O2,dm.cxy(r2, O2[0], O2[1], 1.5*alphat),r'$\mathrm{r_{2}}$')

dm.dimlr(ax,O1,dm.cxy(rb1, O1[0], O1[1], -np.pi/2-.4*alphat),r'$\mathrm{r_{b1}}$')
dm.dimlr(ax,O2,dm.cxy(rb2, O2[0], O2[1], 2*alphat),r'$\mathrm{r_{b2}}$')
plt.savefig('pdf/helicoidal_aparente.pdf', bbox_inches = 'tight',
    pad_inches = 0, transparent=True)
###############################################################################
## PASSOS #####################################################################
###############################################################################
AB = np.sqrt(ra1**2-rb1**2) + np.sqrt(ra2**2-rb2**2) - a*np.sin(alphat)
plt.figure(2)
fig = plt.gcf()
ax = fig.gca(aspect='equal')
plt.axis('off')
plt.plot([0, ra1],[0,0],'k-')
plt.plot([ra1, ra1],[0,b],'k-')
plt.plot([0, ra1],[b,b],'k-')
plt.plot([0, 0],[0,b],'k-')


## hélices
plt.plot([0, -pt+b*np.tan(beta)],[(-pt+b*np.tan(beta))/np.tan(beta),0],'k-')
plt.plot([0, b*np.tan(beta)],[b,0],'k-')
plt.plot([pt, pt+b*np.tan(beta)],[b,0],'k-')
plt.plot([2*pt, 2*pt+b*np.tan(beta)],[b,0],'k-')
plt.plot([3*pt, ra1],[b,-(ra1-(3*pt+b*np.tan(beta)))/np.tan(beta)],'k-')

plt.plot([ra1/2, ra1/2],[-.1*b,1.1*b],'k-.', lw=0.5)
# passos aparente
dm.diml(ax, [b*np.tan(beta), -.1*b], [pt+b*np.tan(beta), -.1*b], r'$p_t$')
plt.plot([b*np.tan(beta), b*np.tan(beta)], [0, -.15*b],'k-', lw=0.5)
plt.plot([b*np.tan(beta)+pt, b*np.tan(beta)+pt], [0, -.15*b],'k-', lw=0.5)

# passo real
dm.diml(ax, [b*np.tan(beta), 0], [(pt+b*np.tan(beta))*np.cos(beta), pt*np.sin(beta)], r'$p_n$')

# passo axial
dm.diml(ax, [pt+b*np.tan(beta), 0], [pt+b*np.tan(beta), pn/np.sin(beta)], r'$p_x$')

dm.circlewx(ax, .5*b, -90, -90 + beta*180/np.pi, 0, b, r'$\beta$')
plt.savefig('pdf/passos.pdf', bbox_inches = 'tight',
    pad_inches = 0, transparent=True)
###############################################################################
## PASSOS BASE ################################################################
###############################################################################
AB = np.sqrt(ra1**2-rb1**2) + np.sqrt(ra2**2-rb2**2) - a*np.sin(alphat)
plt.figure(3)
fig = plt.gcf()
ax = fig.gca(aspect='equal')
plt.axis('off')
plt.plot([0, rb1],[0,0],'k-')
plt.plot([rb1, rb1],[0,b],'k-')
plt.plot([0, rb1],[b,b],'k-')
plt.plot([0, 0],[0,b],'k-')


## hélices
plt.plot([0, -pbt+b*np.tan(betab)],[(-pbt+b*np.tan(betab))/np.tan(betab),0],'k-')
plt.plot([0, b*np.tan(betab)],[b,0],'k-')
plt.plot([pbt, pbt+b*np.tan(betab)],[b,0],'k-')
plt.plot([2*pbt, 2*pbt+b*np.tan(betab)],[b,0],'k-')
plt.plot([3*pbt, rb1],[b,-(rb1-(3*pbt+b*np.tan(betab)))/np.tan(betab)],'k-')

plt.plot([rb1/2, rb1/2],[-.1*b,1.1*b],'k-.', lw=0.5)
# passos aparente
dm.diml(ax, [b*np.tan(betab), -.1*b], [pbt+b*np.tan(betab), -.1*b], r'$p_{bt}$')
plt.plot([b*np.tan(betab), b*np.tan(betab)], [0, -.15*b],'k-', lw=0.5)
plt.plot([b*np.tan(betab)+pbt, b*np.tan(betab)+pbt], [0, -.15*b],'k-', lw=0.5)

# passo real
dm.diml(ax, [b*np.tan(betab), 0], [(pbt+b*np.tan(betab))*np.cos(betab), pbt*np.sin(betab)], r'$p_{bn}$')

dm.circlewx(ax, .5*b, -90, -90 + betab*180/np.pi, 0, b, r'$\beta_b$')
plt.savefig('pdf/passos_base.pdf', bbox_inches = 'tight',
    pad_inches = 0, transparent=True)
###############################################################################
## PASSOS #####################################################################
###############################################################################
plt.figure(4)
fig = plt.gcf()
ax = fig.gca(aspect='equal')
plt.axis('off')
plt.plot([0, rb1],[0,0],'k-')
plt.plot([rb1, rb1],[0,b],'k-')
plt.plot([0, rb1],[b,b],'k-')
plt.plot([0, 0],[0,b],'k-')


## hélices
plt.plot([0, -pbt+b*np.tan(betab)],[(-pbt+b*np.tan(betab))/np.tan(betab),0],'k-')
plt.plot([0, b*np.tan(betab)],[b,0],'k-')
plt.plot([pbt, pbt+b*np.tan(betab)],[b,0],'k-')
plt.plot([2*pbt, rb1],[b,-(rb1-(2*pbt+b*np.tan(betab)))/np.tan(betab)],'k-')
plt.plot([3*pbt, rb1],[b,-(rb1-(3*pbt+b*np.tan(betab)))/np.tan(betab)],'k-')

plt.plot([rb1/2, rb1/2],[-.1*b,1.1*b],'k-.', lw=0.5)
# passos aparente
dm.diml(ax, [b*np.tan(betab), -.1*b], [pbt+b*np.tan(betab), -.1*b], r'$p_{bt}$')
plt.plot([b*np.tan(betab), b*np.tan(betab)], [0, -.15*b],'k-', lw=0.5)
plt.plot([b*np.tan(betab)+pbt, b*np.tan(betab)+pbt], [0, -.15*b],'k-', lw=0.5)

# passo real
dm.diml(ax, [b*np.tan(betab), 1.1*b], [0, 1.1*b], r'$b\cdot\tan\beta_b$')
plt.plot([b*np.tan(betab), b*np.tan(betab)], [0, 1.15*b],'k-', lw=0.5)
plt.plot([0, 0], [b, 1.15*b],'k-', lw=0.5)

dm.circlewx(ax, .5*b, -90, -90 + betab*180/np.pi, 0, b, r'$\beta_b$')
plt.savefig('pdf/razao_conducao_helicoidal.pdf', bbox_inches = 'tight',
    pad_inches = 0, transparent=True)