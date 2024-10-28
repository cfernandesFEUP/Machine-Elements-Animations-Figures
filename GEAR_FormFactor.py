import numpy as np
from scipy import optimize

alpha = 20.0
beta = 0.0
m = 4.5
z = np.array([16., 24.])
x = np.array([0.1817, 0.1715])
b = 14.0
dsh = 30.0

betab = np.arcsin(np.sin(beta)*np.cos(alpha)) 


e_v = 0
alpha = np.radians(alpha)                               # alpha in rad
beta = np.radians(beta)                                 # beta in rad
mt = m/np.cos(beta)                                     # transverse module
d = z*mt                                                # reference diameter
alpha_t = np.arctan(np.tan(alpha)/np.cos(beta))         # transverse pressure angle
betab = np.arcsin(np.sin(beta)*np.cos(alpha))           # base helix angle
pb = np.pi*m*np.cos(alpha)                              # base pitch
pt = np.pi*mt                                           # transverse pitch
pbt = np.pi*m*np.cos(alpha)/np.cos(betab)               # transverse base pitch
db = d*np.cos(alpha_t)                                  # base diameter
u = z[1]/z[0]                                           # transmission ratio
inv_alpha_tw = np.tan(alpha_t) - alpha_t + (2*np.tan(alpha)*np.sum(x)/np.sum(z))
def inv(xx):
    return inv_alpha_tw - np.tan(xx) + xx 
sol  =  optimize.brentq(inv,  0.0,  0.7)
alpha_tw = sol
al = np.sum(db)/(2*np.cos(alpha_tw)) + e_v              # working axis distance    
dl1 = 2*al/(u + 1)                                      # working pitch diameter
dl2 = 2*u*al/(u + 1)
dl = np.array([dl1, dl2])
alpha_tw = np.arccos(np.sum(db)/(2*al))
haP = 1
hfP = 1.25
k = np.sum(z)/2*((((np.tan(alpha_tw) - alpha_tw) - \
                (np.tan(alpha_t) - alpha_t))/np.tan(alpha)) - \
1/np.cos(beta)*(np.cos(alpha_t)/np.cos(alpha_tw) - 1))
da = d + 2*m*(haP + x - k)                              # tip diameter
df = d + 2*m*(x - hfP)                                  # root diameter
rb = db/2
rl = dl/2
rf = df/2
ra = da/2
r = d/2 
alpha_a = np.arccos(db/da);                             # transverse profile angle at tooth tip
epslon_a = z*(np.tan(alpha_a) - np.tan(alpha_tw))/(2*np.pi) # addendum contact ratio
epslon_alpha = np.sum(epslon_a)
epslon_beta = b*np.tan(betab)/pbt
epslon_gama = epslon_beta + epslon_alpha
galphai = rb*(np.tan(alpha_a) - np.tan(alpha_tw))       # length of path of addendum contact individual
galpha = np.sum(galphai)                                # length of path of contact

if beta == 0:
    rfer = 0.375
else:
    rfer = 0.3

rfP = rfer * m

hfP = 1.25 * m

spr = 0

zn1 = z[0]/(np.cos(betab)**2*np.cos(beta))

zn2 = z[1]/(np.cos(betab)**2*np.cos(beta))
 
zn = np.array([zn1,zn2])

dn = m*zn

dbn = dn * np.cos(alpha)

G = rfP/m - hfP/m + x

Es = np.pi/4*m - hfP*np.tan(alpha) + spr/np.cos(alpha)\
    - (1 - np.sin(alpha))*rfP/np.cos(alpha)

H = 2./zn*(np.pi/2 - Es/m) - np.pi/3

def eq1(xx):
    return xx - 2 * G[0] / zn[0] * np.tan(xx) + H[0]

s1 = optimize.brentq(eq1, 0, np.pi / 2.1)

sol1 = s1

vu1 = sol1

def eq2(xx):
    return xx - 2 * G[1] / zn[1] * np.tan(xx) + H[1]

s2 = optimize.brentq(eq2, 0, np.pi / 2.1)

sol2 = s2

vu2 = sol2

vu = np.array([vu1, vu2])

sFn = m*(zn*np.sin(np.pi/3 - vu) + np.sqrt(3) * (G/np.cos(vu) - rfP/m))

dan = dn + 2*ra - 2*r

epslon_alphan = epslon_alpha/(np.cos(betab)**2)

FCT1 = ((dan/2)**2 - (dbn/2)**2)**0.5

FCT2 = np.pi*2*r*np.cos(beta)*np.cos(alpha)/z*(epslon_alphan - 1)

FCT3 = (dbn/2)**2

den = 2 * ((FCT1 - FCT2) ** 2 + FCT3) ** 0.5

alpha_en = np.zeros(2)

gamma_e = np.zeros(2)

hFe = np.zeros(2)

alpha_en[0] = np.arccos(dbn[0] / den[0])

alpha_en[1] = np.arccos(dbn[1] / den[1])

gamma_e[0] = (np.pi/2 + 2*x[0]*np.tan(alpha))/zn[0] + np.tan(alpha) - alpha - np.tan(
    alpha_en[0]) + alpha_en[0]

gamma_e[1] = (np.pi/2 + 2*x[1]*np.tan(alpha))/zn[1] + np.tan(alpha) - alpha - np.tan(
    alpha_en[1]) + alpha_en[1]

alphaFen = alpha_en - gamma_e

hFe[0] = 0.5*m*((np.cos(gamma_e[0]) - np.sin(gamma_e[0])*np.tan(alphaFen[0]))*den[0]/m\
                - zn[0]*np.cos(np.pi/3 - vu[0]) - G[0]/np.cos(vu[0]) + rfP/m)

hFe[1] = 0.5*m*((np.cos(gamma_e[1]) - np.sin(gamma_e[1])*np.tan(alphaFen[1]))*den[1]/m\
                - zn[1]*np.cos(np.pi/3 - vu[1]) - G[1]/np.cos(vu[1]) + rfP/m)

YF = 6*(hFe/ m)*np.cos(alphaFen)/((sFn/m)**2*np.cos(alpha))

