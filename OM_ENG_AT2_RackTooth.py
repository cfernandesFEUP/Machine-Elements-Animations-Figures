## IMPORTAR LIBRARIAS ##
import numpy as np
import matplotlib.pyplot as plt
import sys
## INPUT ##
alpha = float(input(u"Ângulo de Pressão \u00b0:"))
m = float(input("Módulo:"))
z = float(input("Número de dentes:"))
xdesv = float(input("Correcção de Dentado:"))
haP = float(input("Altura de Cabeça da Cremalheira = Altura de Pé da Roda Dentada:"))
hfP = float(input("Altura de Pé da Cremalheira = Altura de Cabeça da Roda Dentada:"))
roh = float(input("Raio da Ferramenta:"))
## DISCRETIZAÇÃO ##
nl = 30
n = 1000
phi = np.linspace(-5*np.pi/z,5*np.pi/z,nl)
alphaR = alpha*np.pi/180
p = np.pi*m
if p<abs(4*m*np.tan(alphaR)):
    sys.exit('Escolher ângulo de pressão mais baixo')
c = haP - hfP
hx = m*xdesv
r = m*z/2
rl = m*z/2+xdesv*m
rb = r*np.cos(alphaR)
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
plt.figure(1)
plt.plot(xRack,yRack)
plt.gca().set_aspect('equal', adjustable='box')
plt.axis('off')

vi = np.array([xRack, yRack, zRack, tRack])
vf = np.zeros((len(vi),len(xRack),len(phi)))

def matrix(i,vi):
    M12 = np.array([[ np.cos(i),  np.sin(i), 0, (r+hx)*np.sin(i)-r*i*np.cos(i)], 
                    [-np.sin(i),  np.cos(i), 0, (r+hx)*np.cos(i)+r*i*np.sin(i)],
                    [         0,          0, 1,                              0],
                    [         0,          0, 0,                              1]])
    return M12.dot(vi)

for i in range(len(phi)):
    vf[:,:,i] = matrix(phi[i],vi)
plt.figure(2)
plt.gca().set_aspect('equal', adjustable='box')
plt.axis('off')
plt.autoscale(enable=True)
plt.plot(vf[0,:,:],vf[1,:,:],'k', lw=0.75)
plt.plot(xAl,yAl,'r:', linewidth = 0.5)
plt.plot(xA,yA,'k-', linewidth = 0.5)
plt.plot(xB,yB,'b-.', linewidth = 0.5)
plt.title('Geração de dentado por uma cremalheira\n'+'z='+str(z)+\
            ', m='+str(m)+', x='+str(xdesv)+'\n'+\
            r'$\alpha$='+str(alpha)+r'$^\circ$,'+' haP='+str(haP)+', hfP='+str(hfP)+r', $\rho$='+ str(roh))
plt.savefig('figures/pdf/dentado_cremalheira.pdf', bbox_inches = 'tight',
    pad_inches = 0, transparent=True)
plt.show()
