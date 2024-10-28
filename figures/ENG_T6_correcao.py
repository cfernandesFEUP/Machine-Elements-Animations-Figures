from scipy import optimize
import numpy as np
## z1 +z2 > 60
alpha = np.pi/9
ratio = np.array([1.1,1.2,1.3,1.4,1.5,2,2.5,3,4,5,10])
z = np.linspace(10,45,36)
m = 2 

def geo(z, ratio, x):
    z1 = int(60/(1+ratio))
    z2 = int(z1*ratio)
    r1, r2 = m*z1/2, m*z2/2 
    rb1, rb2 = r1*np.cos(alpha), r2*np.cos(alpha)
    a = r1+r2
    ra1, ra2 = r1+x*m, r1-x*m
    T2B = np.sqrt(ra2**2-rb2**2)
    T1B = a*np.sin(alpha) - np.sqrt(ra2**2-rb2**2)
    T1A = np.sqrt(ra1**2-rb1**2)
    T2A = a*np.sin(alpha) - np.sqrt(ra1**2-rb1**2)
    return z1, z2, T2B, T1B, T2A, T1A

def gs(x):
    z1, z2, T2B, T1B, T2A, T1A = geo(i, j, x)
    gs1B = abs(1-T2B*z1)/(T1B*z2)
    gs2A = abs(T1A*z2/(T2A*z1)-1)
    return gs2A-gs1B

sol = []
for i in z:
        for j in ratio:
            sol.append(optimize.brentq(gs,  0.0,  0.7))

# def f(x, y):
#     return x * x - 3 + y

# def main():
#     x0 = .1
#     y = 1
#     res = optimize.newton(f, x0, args=(y,))