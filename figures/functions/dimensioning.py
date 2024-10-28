## DIMENSIONING ##
import matplotlib.pyplot as plt
import numpy as np

def diml(ax, xyfrom, xyto, text):
    ax.annotate('',xyfrom,xyto,arrowprops=dict(arrowstyle='<|-|>,head_width=0.1,head_length=0.2',\
                            color='k',shrinkA=0,shrinkB=0,lw=0.5))
    ax.text((xyto[0]+xyfrom[0])/2, (xyto[1]+xyfrom[1])/2, text, fontsize=8,horizontalalignment='center',verticalalignment='bottom')

def dimlv(ax, xyfrom, xyto, text):
    ax.annotate('',xyfrom,xyto,arrowprops=dict(arrowstyle='<|-|>,head_width=0.1,head_length=0.1',\
                            color='k',shrinkA=0,shrinkB=0,lw=0.5))
    ax.text((xyto[0]+xyfrom[0])/2, (xyto[1]+xyfrom[1])/2, text, fontsize=8,horizontalalignment='right',verticalalignment='center')

def dimlnr(ax, xyfrom, xyto):
    ax.annotate('',xyfrom,xyto,arrowprops=dict(arrowstyle='<|-,head_width=0.1,head_length=0.2',\
                            color='k',shrinkA=0,shrinkB=0,lw=0.5))

def dimlr(ax, xyfrom, xyto, text):
    text_angle = 180*np.arctan((xyto[1]-xyfrom[1])/(xyto[0]-xyfrom[0]))/np.pi
    ax.annotate('',xyfrom,xyto,arrowprops=dict(arrowstyle="<|-,head_width=0.1,head_length=0.2",\
                            color='k',shrinkA=0,shrinkB=0,lw=0.5))
    ax.text((xyto[0]+xyfrom[0])/2, (xyto[1]+xyfrom[1])/2, text,\
            rotation=text_angle,horizontalalignment='center')

def dimll(ax, xyfrom, xyto, text, hor, ver,col,lw):
    text_angle = 180*np.arctan((xyto[1]-xyfrom[1])/(xyto[0]-xyfrom[0]))/np.pi
    ax.annotate('',xyfrom,xyto,arrowprops=dict(arrowstyle="<|-,head_width=0.1,head_length=0.2",\
                            color=col,shrinkA=0,shrinkB=0,lw=lw))
    ax.text(xyto[0], xyto[1], text,\
            rotation=text_angle,horizontalalignment=hor,verticalalignment=ver,fontsize=12)

def dimln(ax, xyfrom, xyto, text, form):
    ax.plot([xyfrom[0],xyto[0]],[xyfrom[1],xyto[1]],form,lw=0.5)
    ax.text((xyto[0]+xyfrom[0])/2, (xyto[1]+xyfrom[1])/2, text)

def dimt(xt,yt,text,hp,vp):
    plt.plot(xt, yt, 'k.', ms=3)
    plt.text(xt, yt,text, horizontalalignment=hp,
      verticalalignment=vp,fontsize=8)

def dimtarray(xyt,text,hp,vp):
    plt.plot(xyt[0], xyt[1], 'k.', ms=3)
    plt.text(xyt[0], xyt[1], text, horizontalalignment=hp,
      verticalalignment=vp)
    
def dimtd(xt,yt,text,hp,vp):
    plt.plot(xt, yt, 'k.', ms=3)
    plt.text(xt, yt,text, horizontalalignment=hp,
      verticalalignment=vp,fontsize=10)

def dimttd(xt,yt):
    plt.plot(xt, yt, 'k.', ms=6)

def dimtt(xt,yt,text,hp,vp):
    plt.text(xt, yt,text, horizontalalignment=hp,
      verticalalignment=vp)
    
def circle(r, ang0, ang, x0, y0, col):
    angle = np.linspace(ang0, ang,100)
    angle = angle*np.pi/180
    x = x0 + r*np.cos(angle)
    y = y0 + r*np.sin(angle)
    return plt.plot(x,y,col,lw=0.75)


def circlet(r, ang0, ang, x0, y0, col):
    angle = np.linspace(ang0, ang,100)
    angle = angle*np.pi/180
    x = x0 + r*np.cos(angle)
    y = y0 + r*np.sin(angle)
    return plt.plot(x,y,col,lw=0.5)

def circleb(r, ang0, ang, x0, y0, col):
    angle = np.linspace(ang0, ang,100)
    angle = angle*np.pi/180
    x = x0 + r*np.cos(angle)
    y = y0 + r*np.sin(angle)
    return plt.plot(x,y,col)

def circlec(ax, r, ang0, ang, x0, y0, text):
    angle = np.linspace(ang0, ang,100)
    angle = angle*np.pi/180
    x = x0 + r*np.cos(angle)
    y = y0 + r*np.sin(angle)
    ax.text((x[0]+x[-1])/2, (y[0]+y[-1])/2, text, fontsize=8,horizontalalignment='center',verticalalignment='bottom')
    return plt.plot(x,y,'k:',lw=0.5)


def circlew(ax, r, ang0, ang, x0, y0, text):
    angle = np.linspace(ang0, ang,100)
    angle = angle*np.pi/180
    x = x0 + r*np.cos(angle)
    y = y0 + r*np.sin(angle)
    ax.text(x[-1], y[-1], text, fontsize=8,horizontalalignment='left',verticalalignment='center')
    ax.annotate('',[x[-1],y[-1]],[x[-2],y[-2]],arrowprops=dict(arrowstyle='-|>,head_width=0.1,head_length=0.2',\
                            color='k',shrinkA=0,shrinkB=0,lw=0.1))
    return plt.plot(x,y,'k-',lw=0.5)

def circlewx(ax, r, ang0, ang, x0, y0, text):
    angle = np.linspace(ang0, ang,100)
    angle = angle*np.pi/180
    x = x0 + r*np.cos(angle)
    y = y0 + r*np.sin(angle)
    delta = np.pi*(ang+ang0)/(2*180)
    xc = x0 + r*np.cos(delta)
    yc = y0 + r*np.sin(delta)
    ax.text(xc, yc, text, fontsize=8,horizontalalignment='center',verticalalignment='bottom')
    ax.annotate('',[x[-1],y[-1]],[x[-2],y[-2]],arrowprops=dict(arrowstyle='-|>,head_width=0.1,head_length=0.2',\
                            color='k',shrinkA=0,shrinkB=0,lw=0))
    ax.annotate('',[x[0],y[0]],[x[1],y[1]],arrowprops=dict(arrowstyle='-|>,head_width=0.1,head_length=0.2',\
                            color='k',shrinkA=0,shrinkB=0,lw=0))
    return plt.plot(x,y,'k-',lw=0.5)


def circlewxb(ax, r, ang0, ang, x0, y0, text):
    angle = np.linspace(ang0, ang,100)
    angle = angle*np.pi/180
    x = x0 + r*np.cos(angle)
    y = y0 + r*np.sin(angle)
    delta = np.pi*(ang+ang0)/(2*180)
    xc = x0 + r*np.cos(delta)
    yc = y0 + r*np.sin(delta)
    ax.text(xc, yc, text, fontsize=8,horizontalalignment='center',verticalalignment='top')
    ax.annotate('',[x[-1],y[-1]],[x[-2],y[-2]],arrowprops=dict(arrowstyle='-|>,head_width=0.1,head_length=0.2',\
                            color='k',shrinkA=0,shrinkB=0,lw=0))
    ax.annotate('',[x[0],y[0]],[x[1],y[1]],arrowprops=dict(arrowstyle='-|>,head_width=0.1,head_length=0.2',\
                            color='k',shrinkA=0,shrinkB=0,lw=0))
    return plt.plot(x,y,'k-',lw=0.5)


def circlewxh(ax, r, ang0, ang, x0, y0, text):
    angle = np.linspace(ang0, ang,100)
    angle = angle*np.pi/180
    x = x0 + r*np.cos(angle)
    y = y0 + r*np.sin(angle)
    delta = np.pi*(ang+ang0)/(2*180)
    xc = x0 + r*np.cos(delta)
    yc = y0 + r*np.sin(delta)
    ax.text(xc, yc, text, fontsize=8,horizontalalignment='left',verticalalignment='center')
    ax.annotate('',[x[-1],y[-1]],[x[-2],y[-2]],arrowprops=dict(arrowstyle='-|>,head_width=0.1,head_length=0.2',\
                            color='k',shrinkA=0,shrinkB=0,lw=0))
    ax.annotate('',[x[0],y[0]],[x[1],y[1]],arrowprops=dict(arrowstyle='-|>,head_width=0.1,head_length=0.2',\
                            color='k',shrinkA=0,shrinkB=0,lw=0))
    return plt.plot(x,y,'k-',lw=0.5)

def cxy(radius, x0, y0, phi):
    return [x0 + radius*np.cos(phi), y0 + radius*np.sin(phi)]

def dxy(radius, x0, y0, phi):
    return [x0 + radius*np.sin(phi), y0 + radius*np.cos(phi)]


def dimang(ax, x1, y1, x2, y2, connectionstyle):
    ax.annotate("",
                xy=(x1, y1), xycoords='data',
                xytext=(x2, y2), textcoords='data',
                arrowprops=dict(arrowstyle="<|-|>", color="black",
                                shrinkA=0, shrinkB=0,
                                patchA=None, patchB=None,
                                connectionstyle=connectionstyle,
                                ),
                )

def dimangr(ax, x1, y1, x2, y2, connectionstyle):
    ax.annotate("",
                xy=(x1, y1), xycoords='data',
                xytext=(x2, y2), textcoords='data',
                arrowprops=dict(arrowstyle="<|-", color="black",
                                shrinkA=0, shrinkB=0,
                                patchA=None, patchB=None,
                                connectionstyle=connectionstyle,
                                ),
                )

def ev(rb,phi):
    x = rb*(np.sin(phi)-phi*np.cos(phi))
    y = rb*(np.cos(phi)+phi*np.sin(phi))
    return np.array([x, y])

def inv(alpha):
    return np.tan(alpha)-alpha

def matrix(i,vi):
    M12 = np.array([[ np.cos(i),  np.sin(i)], 
                    [-np.sin(i),  np.cos(i)]])
    return M12.dot(vi)

def secant(x1,x2,y1,y2,xlim1,xlim2):
    coef = np.polyfit([x1,x2],[y1,y2],1)
    poly = np.poly1d(coef)
    xline = np.linspace(xlim1,xlim2)
    yline = poly(xline)
    plt.plot(xline, yline,'k:',lw=0.5)

def rack(m,z,rb,xdesv,hfP,haP,roh,alpha):
    nl = 30
    n = 100
    p = np.pi*m
    c = haP - hfP
    ylU = np.linspace(0,hfP*m,n)
    xlU = ylU*np.tan(alpha)-p/4
    ylB = np.linspace(0,-hfP*m,n)
    xlB = -p/4+ylB*np.tan(alpha)
    xl = np.concatenate((xlB[::-1],xlU), axis=0)
    yl = np.concatenate((ylB[::-1],ylU), axis=0)
    xt, yt = np.linspace(xl[-1],-xl[-1],n//2), np.linspace(yl[-1], yl[-1],n//2)
    xc, yc = np.sqrt((roh*m)**2-((roh-c)*m)**2), yl[0]-c*m + roh*m
    yfl = np.linspace(yl[0]-(haP-1)*m,yl[0],n//2)
    xfl = ((roh*m)**2-(yfl-yc)**2)**.5+(xl[0]-xc)
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
    return vi, vp