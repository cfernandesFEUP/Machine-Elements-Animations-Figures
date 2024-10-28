import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import scienceplots

plt.style.use(['science','ieee'])

## DIMENSIONS

inch = 25.4 #mm
ro = 6*inch
ri = 2*inch
h = ro-ri
t = 0.75*inch

## AREA

A = h*t

## LOAD

F = 5000*4.44822162

## CENTROIDAL AXIS

r = ri + h/2

## NEUTRAL AXIS

rn = h/np.log(ro/ri)

## DISTANCE BETWEEN BOTH AXIS

e = r - rn

## MOMENT
M = F*r

## LIMITS Y

y = np.linspace(h/2,-h/2,201)

## BENDING STRESSES

sigmaB = M*(y-e)/(e*A*(r-y))

## NORMAL STRESSES

sigmaN = F/A


## TOTAL STRESS

sigma = sigmaB + sigmaN

## 1st ORDER FEM
dfL = pd.read_csv('dados/cranehookLIN.csv')
xL = dfL[dfL.columns[16]]
sigmaFEML = dfL[dfL.columns[2]]

## 2nd ORDER FEM
df2 = pd.read_csv('dados/cranehook2ND.csv')
x2 = df2[df2.columns[16]]
sigmaFEM2 = df2[df2.columns[2]]

## 2nd ORDER FEM
dfR = pd.read_csv('dados/cranehook2OR.csv')
xR = dfR[dfR.columns[16]]
sigmaFEMR = dfR[dfR.columns[2]]

plt.figure(1)
plt.plot((r-y)/inch,sigma,'k-',label='Analytical')
plt.plot(xL/inch,sigmaFEML,'k:',label='FEM 1st order')
plt.plot(xL/inch,sigmaFEM2,'k-.',label='FEM 2nd order')
plt.plot(xL/inch,sigmaFEMR,'k--',label='FEM 2nd order reduced integration')
plt.ylabel(r'$\sigma$ / MPa')
plt.xlabel('r / inch')
#plt.grid()
plt.legend()
plt.savefig('pdf/cranehookFEM.pdf', bbox_inches = 'tight',
    pad_inches = 0.1, transparent=True)
