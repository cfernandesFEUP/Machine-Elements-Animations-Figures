import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize

## MATERIAL
sigmaCed = 324
SF = 2
sigmaE = sigmaCed/SF

## DIMENSIONS
ro = 150
ri = 50
h = ro-ri
R = 50

## AREA
A = np.pi*R**2

## CENTROIDAL AXIS

r = ri + h/2

## NEUTRAL AXIS

rn = (ro**(1/2) + ri**(1/2))**2/4

## DISTANCE BETWEEN BOTH AXIS

e = r - rn


## LIMITS Y
h2 = (h/2)
h1 = (h/2)

y = np.linspace(-h2,h1,100)


def solution(F,r,e,A,y,sigmaE):
    ## MOMENT
    M = F*r

    ## BENDING STRESSES

    sigmaB = M*(y-e)/(e*A*(r-y))

    ## NORMAL STRESSES

    sigmaN = F/A

    ## TOTAL STRESS

    sigma = sigmaB + sigmaN
    
    return sigma[-1] - sigmaE


F = optimize.brentq(solution, 40000., 100000., args=(r, e, A, y, sigmaE))


sigma0 = solution(F,r,e,A,np.array([0]),sigmaE)

## MOMENT
M = F*r

## BENDING STRESSES

sigmaB = M*(y-e)/(e*A*(r-y))

## NORMAL STRESSES

sigmaN = F/A

## TOTAL STRESS

sigma = sigmaB + sigmaN


plt.figure(1)
plt.plot(y,sigma,'k-',label='Analytical')
plt.ylabel(r'$\sigma$ / MPa')
plt.xlabel('r / inch')
plt.grid()
plt.legend()