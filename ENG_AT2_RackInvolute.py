## IMPORTAR LIBRARIAS ##
import numpy as np
import figures.functions.dimensioning as dm
import matplotlib.animation as anm
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
p = np.pi*m
###############################################################################
## CREMALHEIRA 2 ##############################################################
###############################################################################
vi, vp = dm.rack(m,z2,rb2,x2,hfP,haP,roh,alpha)
xra = hfP*m*np.tan(alpha)
yra = hfP*m
yta = haP*m
Tx = r1+rb1*np.sin(alpha)
Ty = r1-rb1*np.cos(alpha)
begin = 250
end = -250
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

phiEc = np.linspace(-alphaA,-np.pi/2.6,n)
xi1c =  dm.ev(rb1,phiEc)[0]
yi1c =  dm.ev(rb1,phiEc)[1]
vi1c =  np.array([xi1c,yi1c])

ang = np.linspace(alpha1 + dm.inv(alpha1),0,20)

from matplotlib.backends.backend_pdf import PdfPages
pp = PdfPages('figures/pdf/anim/cremalheira.pdf')

def init():
        line.set_data([], [])
        linep.set_data([], []) 
        lineq.set_data([], [])
        pointM.set_data([], [])
        return line, linep, lineq, pointM,
    
def animate(i):
    delta = ang[i]*r1
    xp = vi[0,begin:end]+delta
    yp = -vi[1,begin:end]
    vf1 = dm.matrix(phi1-ang[i],vi1)
    vfc1 = dm.matrix(phi1-ang[i],vi1c)
    vM =  dm.ev(rb1,-(alpha1+dm.inv(alpha1))+ang[i])
    vfM = dm.matrix(phi1-ang[i],vM)
    line.set_data(xp, yp)
    linep.set_data(vf1[0]+r1,vf1[1]+r1)
    lineq.set_data(vfc1[0]+r1,vfc1[1]+r1)
    pointM.set_data(vfM[0]+r1, vfM[1]+r1)
    pp.savefig(fig, bbox_inches = 'tight', transparent=True)
    return line, linep, lineq, pointM,

## PLOT ##
plt.figure(2)
fig = plt.gcf()
ax = fig.gca()
ax.set_aspect("equal")
plt.axis('off')
plt.xlim([-.5*r1,2*r1])
plt.ylim=([-.75*r1,r1])
line, = ax.plot([], [], 'k')
linep, = ax.plot([], [], 'r')
lineq, = ax.plot([], [], 'k',lw=0.5)
pointM, = ax.plot([], [], 'k.',ms=6)
plt.plot([-0.5*r1,1.5*r1], [0,0],'k',lw=0.5)
## secantes
dm.secant(Tx,r1,Ty,0,0,1.5*r1)
## raio de base
dm.circle(rb1, -25, -155, r1, r1,'k--')
dm.circle(r1, -25, -155, r1, r1,'k')
dm.circle(ra1, -25, -155, r1, r1,'k:')
## cotagem
dm.dimln(ax,O1, [Tx,Ty],'','k-')
dm.dimt(O1[0],O1[1],r'$\mathrm{O}$','left','bottom')
dm.dimt(Tx,Ty,r'$\mathrm{T}$','left','bottom')
dm.dimt(r1,0,r'$\mathrm{I}$','left','bottom')

ani = anm.FuncAnimation(fig, animate, init_func=init,
                                   frames=len(ang), blit=True, repeat=False, interval=100)
plt.show()
pp.close()