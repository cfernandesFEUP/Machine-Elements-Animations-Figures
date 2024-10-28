import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as anm
m = 2.
z = 20.
r = m*z/2
rb = r*np.cos(np.pi/9)
x = 0.1817
alpha = np.pi/9
n= 50

phi = np.linspace(0,.9*np.pi/2,n)
phiC = np.linspace(-np.pi/2,np.pi/2,100)
phase = 0#-np.pi/4
xp, yp = [], []
bap, bar, mcal, angle = [], [], [], []
## involuta
xINV = rb*(np.sin(phi+phase)-phi*np.cos(phi+phase))
yINV = rb*(np.cos(phi+phase)+phi*np.sin(phi+phase))
# dydx = np.diff(yINV)/np.diff(xINV)

from matplotlib.backends.backend_pdf import PdfPages
pp = PdfPages('figures/pdf/anim/involute.pdf')

def init():
    line.set_data([], [])    
    linet.set_data([], []) 
    liner.set_data([], []) 
    linep.set_data([], [])
    linemt.set_data([], [])
    linemr.set_data([], [])
    linemm.set_data([], [])
    return line, linet, liner, linep, linemt, linemr, linemm,
    
def animate(i):
    xp.append(xINV[i])
    yp.append(yINV[i])
    x1, y1 = xp[-1], yp[-1]
    alphaM = phi[i] + phase
    ## TANGENTE AO CIRCULO DE BASE
    Tx = rb*np.sin(alphaM)
    Ty = rb*np.cos(alphaM)
    ## TANGENTE À EVOLVENTE
    leng = 6
    xi = x1 - leng*np.sin(alphaM)
    yi = y1 - leng*np.cos(alphaM)
    xe = x1 + leng*np.sin(alphaM)
    ye = y1 + leng*np.cos(alphaM)
    br = (y1-Ty)/(x1-Tx)
    be = (x1-xe)/(y1-ye)
    mcal.append(br/be)
    bar.append(br)
    bap.append(1/be)
    angle.append(phi[i])
    line.set_data(xp, yp)
    linet.set_data([Tx,x1],[Ty,y1])
    liner.set_data([0,Tx], [0,Ty])
    linep.set_data([xi,xe], [yi,ye])
    linemt.set_data(angle, bap)
    linemr.set_data(angle, bar)
    linemm.set_data(angle, mcal)
    pp.savefig(fig, bbox_inches = 'tight', transparent=True)
    return line, linet, liner, linep, linemt, linemr, linemm,

## circulo de base
xA, yA = rb*np.sin(phiC), rb*np.cos(phiC)

# th0 = np.arccos(xINV[0]/rb)*180/np.pi
# th1 = np.arccos(x1/np.sqrt(x1**2+y1**2))*180/np.pi
# th2 = np.arccos(Tx/rb)*180/np.pi

fig = plt.figure(2)
ax = fig.add_subplot(211, aspect='equal', autoscale_on=True, xlim= (-rb, rb),ylim=(-0.1*rb,2*rb))
plt.axis('off')
plt.plot(xA,yA,'k', lw=0.5)
line, = ax.plot([], [], 'r', lw=2.5)
linet, = ax.plot([], [], 'g-', lw = 0.75)
liner, = ax.plot([], [], 'k-', lw = 0.5)
linep, = ax.plot([], [], 'b-', lw = 0.75)

ax = fig.add_subplot(212, autoscale_on=True)
plt.yscale('symlog')
plt.ylim([-10, 30])
plt.xlim([0, 1.4])
plt.xlabel('Ângulo / rad')
# ax.set_yscale("log")
linemt, = ax.plot([], [], 'b-', label=r'Declive de $t_e$')
linemr, = ax.plot([], [], 'g-', label=r'Declive $t_r$')
linemm, = ax.plot([], [], 'k-', label=r' $\frac{m t_r}{m t_e}$', lw=2)
ax.legend()
# ANIMATION ##
ani = anm.FuncAnimation(fig, animate, init_func=init,
                                   frames=(len(phi)-1), blit=False, repeat=False, interval=100)
plt.show()

pp.close()