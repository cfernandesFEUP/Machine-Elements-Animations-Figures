'''MIT License

Copyright (c) 2021 Carlos M.C.G. Fernandes

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE. '''

import sys
sys.dont_write_bytecode = True
import numpy as np
from functions import MAAG, involute
## GEAR SELECTION ##################################################################
alpha = 20.0
beta = 0.0
m = 2.
z = np.array([30., 30.])
x = np.array([[-0.4, -0.4],[0.,0.],[0.4,0.4]])
b = 10.0
dsh = 20.

phi = 7.5

color = ['r', 'k', 'b']
corr = ['x<0','x=0', 'x>0']
## DRAW CIRCLES #####################################################################
def circlew(r, ang0, ang, x0, y0):
    angle = np.linspace(ang0, ang,100)
    angle = angle*np.pi/180
    x = x0 + r*np.cos(angle)
    y = y0 + r*np.sin(angle)
    return x, y

## MATRIX TRANSFORMATION ###########################################################
def matrix(i,vi):
    M12 = np.array([[ np.cos(i),  np.sin(i)], 
                    [-np.sin(i),  np.cos(i)]])
    return M12.dot(vi)


import matplotlib.pyplot as plt 
plt.figure(1002)
for i in range(len(x)):
    ## MAAG CALCULATION ##
    mt, pt, pb, pbt, betab, al, r, rl, ra, rb, rf, alpha_t, alpha_tw, epslon_alpha,\
    epslon_a, epslon_beta, epslon_gama, galpha, galphai, Req, u, T1T2, T1A, T2A, \
    AB, AC, AD, AE, rA1, rA2, rB1, rB2, rD1, rD2 = MAAG.calc(alpha, beta, m, z, x[i], b)
    
    ## TOOTH THICKNESS #################################################################
    s0 = np.pi/2
    s = s0 + 2*x[i]*m*np.tan(alpha*np.pi/180) 
    ds = (s-s0)/2
    ## GEAR GEOMETRY AND MESH ##########################################################
    n = 1000
    xGEO, yGEO, xRF, yRF, xI, yI = involute.geo(z,alpha_t,beta,rb,ra,rf,x[i],m,mt,n,al,0)
    ##
    rot = ds/r
    xyF = matrix(rot[0],np.array([xGEO,yGEO]))
## PLOT ############################################################################
    plt.plot(xyF[0], xyF[1], color[i],label=corr[i])
x, y  = circlew(r[0],90+phi,90-phi,0,0)
plt.plot(x,y,'k-.',label='círculo primitivo',lw=0.5)
x, y  = circlew(rb[0],90+phi,90-phi,0,0)
plt.plot(x,y,'k--',label='círculo de base',lw=0.5)
plt.legend(loc=0, prop={'size': 8})
plt.axis('off')
plt.axis('equal')
plt.savefig('pdf/corr_profiles.pdf', bbox_inches = 'tight',
pad_inches = 0, transparent=True)