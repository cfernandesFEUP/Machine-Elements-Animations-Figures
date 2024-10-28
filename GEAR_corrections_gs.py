from math import atan, tan, cos, sin, pi, sqrt
from scipy import optimize

z1 = 25
z2 = 40
m = 2
alpha = pi/9
beta = 0

alphat = atan(tan(alpha)/cos(beta)) 

r1 = m*z1/(2*cos(beta))
r2 = m*z2/(2*cos(beta))

rb1 = r1*cos(alphat)
rb2 = r2*cos(alphat)

u = z2/z1

def COR(x, r1, r2, alphat, u):
    x1 = x
    x2 = -x
    a = r1 + r2
    ra1 = r1 + m*(1+x1)
    ra2 = r2 + m*(1+x2)
    T1T2 = a*sin(alphat)
    T1B = sqrt(ra1**2-rb1**2)
    T2B = T1T2 - T1B
    T2A = sqrt(ra2**2-rb2**2)
    T1A = T1T2 - T2A
    gs1B = abs(1 - (1/u)*T2B/T1B)
    gs2A = abs(u*T1A/T2A - 1)
    return gs1B - gs2A

sol = optimize.brentq(COR, 0, 1., args=(r1, r2, alphat, u))

def WK(z, alpha, m, x):
    s = pi*m/2 + 2*x*m*tan(alpha)
    k = round(z*alpha/pi + 0.5,0)
    print(k,s)
    Wk = m*cos(alpha)*((k-1)*pi + s/m + z*(tan(alpha)-alpha))
    return Wk

Wk1 = WK(z1, alpha, m, sol)
Wk2 = WK(z2, alpha, m, -sol)

print('x1:', sol)
print('x2:', -sol)
print('Wk1:', Wk1)
print('Wk2:', Wk2)