from math import atan, tan, cos, acos, pi
from scipy import optimize

z = float(input('Número de dentes - z: '))

m = float(input('Módulo - m / mm: '))

alpha = float(input('Ângulo de pressão - \u03B1 / \xb0: '))*pi/180

DM = float(input('Dimensão dos calibres - DM / mm: '))

Md = float(input('Dimensão sobre os calibres - Md / mm: '))

beta = 0

alphat = atan(tan(alpha)/cos(beta))

r = m*z/(2*cos(beta))

rb = r*cos(alphat)

ra = r + m

rd = r -1.25*m

ld = (Md-DM)/2

alphaK = acos(rb/ld)

def INV(angle):
    return tan(angle) - angle

def FUNCTION(x, DM, z, m, alpha, alphat, alphaK):
    RES = INV(alphat) + DM/(z*m*cos(alphat)) - (pi - 4*x*tan(alpha))/(2*z) - INV(alphaK)
    return RES

sol = optimize.newton(FUNCTION, 0.1, args=(DM, z, m, alpha, alphat, alphaK))

rdl = rd + sol*m
ral = ra + sol*m

print('r: ', r)
print('rb: ', rb)
print('ra: ', ra)
print('rd: ', rd)
print('x:', sol)
print('rd*:', rdl)
print('ra*:', ral)