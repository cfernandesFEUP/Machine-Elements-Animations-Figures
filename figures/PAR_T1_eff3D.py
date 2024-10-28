from sympy import *
import numpy as np
Mt, H, dm, F, tanp, gama, beta, betan, eta, N, P, Wf, Wu = symbols('M_t H d_m F tan_phi gamma \
beta beta_n eta N P W_f W_u')

eq1 = Eq(Mt, H*dm/2)
display(eq1)

eq2 = Eq(H*cos(gama), F*sin(gama)+tanp*(F*cos(gama)+H*sin(gama)))
display(eq2)
res1 = simplify(solve(eq2,H)[0])
display(res1)


eq1.subs(H,res1)


eq3 = Eq(tan(betan),tan(beta)*cos(gama))
display(eq3)

eq4 = Eq(N,F*cos(gama)/cos(betan)+H*sin(gama)/cos(betan))
display(eq4)
eq5 = Eq(H*cos(gama), F*sin(gama)+tanp*(eq4.rhs))
display(eq5)
res2 = simplify(solve(eq5,H)[0])
display(res2)

eqmt = eq1.subs(H,res2)
display(eqmt)

eq6 = Eq(tan(gama),P/(2*pi*dm/2))
display(eq6)
res3 = solve(eq6,P)[0]

eq7 = Eq(Wf,Mt*2*pi)
display(eq7)

eq8 = Eq(Wu,F*P)
display(eq8)

eqe = Eq(eta, Wu/Wf)
display(eqe)

eqe2 = eqe.subs(Wu,eq8.rhs).subs(Wf,eq7.rhs)
display(eqe2)

expr = simplify(eqe2.subs(Mt,eqmt.rhs).subs(P,res3))
display(expr)

Pa = 1.5
D = 10
di = 8.16 
rm = (D+di)/2
gam = np.arctan(Pa/(2*np.pi*rm))
val = expr.subs(gama,gam)
fun = lambdify([betan, tanp], val.rhs)

size = 50
ba = np.linspace(0,60,size)
bn = np.arctan(np.tan(ba*np.pi/180)*np.cos(gam))
tp = np.linspace(0.0,0.2,size) 
bplot = np.outer(bn,np.ones(size)).T
tplot = np.outer(tp,np.ones(size))
zplot = fun(bplot,tplot)

from mpl_toolkits import mplot3d
# %matplotlib inline
import matplotlib.pyplot as plt
# fig = plt.figure(figsize=(20,10))
plt.figure(1)
ax = plt.axes(projection='3d')
ax.plot_wireframe(tplot, bplot*180/np.pi, zplot, color='black')
#ax.plot3D(tplot, bplot, zplot, 'gray')
#ax = plt.axes(projection='3d')
ax.set_ylabel(r'$\beta_n$ / $^{\circ}$')
ax.set_xlabel(r'$tan \phi$')
ax.set_zlabel(r'$\eta$ / $\%$')
ax.set_facecolor('w')
plt.tick_params(axis='both', which='major')
# ax.view_init(30, 30)
plt.savefig('pdf/PAR_eff3D.pdf', bbox_inches = 'tight',
    pad_inches = 0, transparent=True)
