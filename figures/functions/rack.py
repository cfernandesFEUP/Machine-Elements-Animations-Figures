import matplotlib.pyplot as plt
import matplotlib.animation as anm
import numpy as np
import sys
plt.rcParams["animation.html"] = "jshtml"

def animation(alpha, haP, hfP, roh, m, z, xdesv, nl, n):
    phi = np.linspace(-6*np.pi/z,6*np.pi/z,nl)
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
    xra, yra = ra*np.sin(phi), ra*np.cos(phi)
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
    x0 = x1-p
    y0 = y1
    x2 = x1+p
    y2 = y1
    x01 = np.linspace(x0[-1],x1[0],n//3)
    y01 = np.linspace(y0[-1],y1[0],n//3)
    x12 = np.linspace(x1[-1],x2[0],n//3)
    y12 = np.linspace(y1[-1],y2[0],n//3)
    xRack = np.concatenate((x0, x01, x1, x12, x2), axis=0)
    yRack = np.concatenate((y0, y01, y1, y12, y2), axis=0)
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
    ## PLOT ##
    fig = plt.figure(2)
    ax = fig.add_subplot(111, aspect='equal', autoscale_on=False, xlim=(-ra,ra), ylim=(0,1.5*ra))
    plt.plot(xAl,yAl,'r:', linewidth = 1, label = 'Tangente à reta de referência')
    plt.plot(xA,yA,'k-', linewidth = 1,label = 'Primitivo de Corte')
    plt.plot(xB,yB,'b-.', linewidth = 1, label = 'Círculo de Base')
    plt.plot(xra,yra,'g:', linewidth = 1, label = 'Círculo de Cabeça')
    ax.axis('off')
    plt.legend(loc=3)
    plt.title('Geração de dentado por uma cremalheira')
    plt.suptitle('Engrenagem: '+'z='+str(z)+', m='+str(m)+r', x='+str(xdesv)+'\n'+
              'Cremalheira: '+ r'$\alpha$='+str(alpha)+r'$^\circ$,'+' ha='+str(haP)+', hf='+str(hfP)+
              r', $\rho$='+ str(roh)+'\n',y=0.45)
    line, = ax.plot([], [], 'k',lw=0.5)
    linep, = ax.plot([], [], 'r-')
    # ANIMATION ##
    ani = anm.FuncAnimation(fig, animate, init_func=init,
                                   frames=len(phi), blit=False, interval=100)
    
    plt.show()
    return ani