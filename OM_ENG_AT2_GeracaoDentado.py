import matplotlib.pyplot as plt
import matplotlib.animation as anm
import numpy as np
import sys

## INPUT ##
alpha = 20
m = float(input("Módulo:"))
z = float(input("Número de dentes:"))
xdesv = float(input("Correcção de Dentado:"))
haP = 1.25
hfP = 1.
roh = 0.38
## DISCRETIZAÇÃO ##
nl = 30
n = 1000

phi = np.linspace(-6*np.pi/z,6*np.pi/z,nl)
phi2 = np.linspace(-3*np.pi/z,3*np.pi/z,nl)
alphaR = alpha*np.pi/180
p = np.pi*m
if p<abs(4*m*np.tan(alphaR)):
    sys.exit('Escolher ângulo de pressão mais baixo')
c = haP - hfP
hx = m*xdesv
r = m*z/2
rl = m*z/2+xdesv*m
ra = rl + m
rb = r*np.cos(alphaR)
xra, yra = ra*np.sin(phi2), ra*np.cos(phi2)
xri, yri = 0.4*ra*np.sin(phi2), 0.4*ra*np.cos(phi2)
xA, yA = r*np.sin(phi), r*np.cos(phi)
xAl, yAl = rl*np.sin(phi), rl*np.cos(phi)
xB,yB = rb*np.sin(phi), rb*np.cos(phi)
ylU = np.linspace(0,hfP*m,n)
xlU = ylU*np.tan(alphaR)-p/4
ylB = np.linspace(0,-hfP*m,n)
xlB = -p/4+ylB*np.tan(alphaR)
xl = np.concatenate((xlB[::-1],xlU), axis=0)
yl = np.concatenate((ylB[::-1],ylU), axis=0)
xt, yt = np.linspace(xl[-1],-xl[-1],n//2), np.linspace(yl[-1], yl[-1],n//2)
yc = yl[0] -c*m + roh*m
teta = np.arcsin((yc-yl[0])/(roh*m))
xc = xl[0] - roh*m*np.cos(teta)
if xc<-p/2:
    xfl = np.linspace(-p/2,xl[0],n//2)
else:
    xfl = np.linspace(xc,xl[0],n//2)     
yfl = yc-((roh*m)**2-(xfl-xc)**2)**.5
xr, yr = -xl[::-1], yl[::-1]
xfr, yfr = -xfl[::-1], yfl[::-1]
x1 = np.concatenate((xfl, xl, xt, xr, xfr), axis=0)
y1 = np.concatenate((yfl, yl, yt, yr, yfr), axis=0)
x_1 = x1-2*p
y_1 = y1
x0 = x1-p
y0 = y1
x2 = x1+p
y2 = y1
x3 = x1+2*p
y3 = y1
x_10 = np.linspace(x_1[-1],x0[0],n//3)
y_10 = np.linspace(y_1[-1],y0[0],n//3)
x01 = np.linspace(x0[-1],x1[0],n//3)
y01 = np.linspace(y0[-1],y1[0],n//3)
x12 = np.linspace(x1[-1],x2[0],n//3)
y12 = np.linspace(y1[-1],y2[0],n//3)
x23 = np.linspace(x2[-1],x3[0],n//3)
y23 = np.linspace(y2[-1],y3[0],n//3)
xRack = np.concatenate((x_1, x_10, x0, x01, x1, x12, x2, x23, x3), axis=0)
yRack = np.concatenate((y_1, y_10, y0, y01, y1, y12, y2, y23, y3), axis=0)
zRack = np.zeros(len(xRack))
tRack = np.ones(len(xRack))
vi = np.array([xRack, yRack, zRack, tRack])
xprim = np.linspace(min(xRack),max(xRack),nl)
yprim = np.zeros(nl)
zprim = np.zeros(nl)
tprim = np.ones(nl)
vp = np.array([xprim, yprim, zprim, tprim])
vf = np.zeros((len(vi),len(xRack),len(phi)))
vfp = np.zeros((len(vp),len(xprim),len(phi)))
vl = np.array([xRack, m*np.ones(len(xRack)), zRack, tRack])
vup = np.zeros((len(vl),len(xRack),len(phi)))

fig = plt.figure()
# fig.canvas.set_window_title('Geração de Dentado')
ax = fig.add_subplot(111, aspect='equal', autoscale_on=False, xlim=(-ra,ra), ylim=(0,1.5*ra))
plt.title('Geração de dentado por uma cremalheira')
ax.fill(np.append(xra,xri[::-1]), np.append(yra,yri),'b')

# xp, yp = [], []
def init():
    line.set_data([], [])    
    linep.set_data([], [])
    return line, linep,

def animate(i):
    xp = vf[0,:,i]
    yp = vf[1,:,i]
    line.set_data(xp, yp)
    linep.set_data(vfp[0,:,i], vfp[1,:,i])
    ax.fill(np.append(xp, vup[0,::-1,i]), np.append(yp, vup[1,::-1,i]), 'w')
    return line, linep,

## DEFINIÇÃO DA MATRIZ DE TRANSFORMAÇÃO ##
def matrix(i,vi):
    M12 = np.array([[ np.cos(i),  np.sin(i), 0, (r+hx)*np.sin(i)-r*i*np.cos(i)], 
                    [-np.sin(i),  np.cos(i), 0, (r+hx)*np.cos(i)+r*i*np.sin(i)],
                    [         0,          0, 1,                              0],
                    [         0,          0, 0,                              1]])
    return M12.dot(vi)

## APLICAÇÃO DA MATRIZ ##
for i in range(len(phi)):
    vf[:,:,i] = matrix(phi[i],vi)
    vfp[:,:,i] = matrix(phi[i],vp)
    vup[:,:,i] = matrix(phi[i],vl)
    
## PLOT ##
plt.plot(xAl,yAl,'k:', linewidth = 1, label = 'Primitivo  (x=0)')
plt.plot(xA,yA,'r-', linewidth = 1,label = 'Primitvo')
plt.plot(xB,yB,'b-.', linewidth = 1, label = 'Círculo de Base')

ax.axis('off')
plt.legend(loc=3)

plt.suptitle('Roda dentada2: '+'z='+str(z)+', m='+str(m)+r', x='+str(xdesv)+'\n'+
          'Cremalheira: '+ r'$\alpha$='+str(alpha)+r'$^\circ$,'+' ha='+str(haP)+', hf='+str(hfP)+
          r', $\rho$='+ str(roh)+'\n',y=0.1)
line, = ax.plot([], [], 'k')
linep, = ax.plot([], [], 'k:')
# ANIMATION ##
figManager = plt.get_current_fig_manager()
figManager.window.showMaximized()

ani = anm.FuncAnimation(fig, animate, init_func=init,
                               frames=len(phi), blit=False, interval=500)

plt.show()