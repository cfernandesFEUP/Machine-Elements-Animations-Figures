import matplotlib.pyplot as plt
import numpy as np
# import matplotlib as mpl
# mpl.rcParams["font.size"] = 11

v = 0.3

K = np.linspace(1,4,5000)

smax = (K**2-1)/(K**2+1)

emax = (K**2-1)/((1-2*v)+(1+v)*K**2)

emax0 = (K**2-1)/((1-v)+(1+v)*K**2)

tmax = (K**2-1)/(2*K**2)

EDmax = (K**2-1)/(np.sqrt(3)*K**2)

EDmax0 = (K**2-1)/np.sqrt(1+3*K**4)

ISO = 2*(K-1)/(K+1)

FP = np.log(K)

Hamburg = K-1

plt.figure(1)
plt.plot(K,Hamburg,'k-',label=r'$\sigma_{max}$')
plt.plot(K,ISO,'k-.',label=r'$\tau_{max}$')
plt.plot(K,smax,'b-',label=r'$\sigma_{max}$')
plt.plot(K,tmax,'b:',label=r'$\tau_{max}$')
plt.plot(K,emax,'b--',label=r'$\varepsilon_{max}$')
plt.plot(K,EDmax,'b-.',label=r'$ED_{max}$')
plt.ylim((0,1.6))
plt.xlabel('K')
plt.ylabel(r'$\frac{p_i}{\sigma_{y}}$')
plt.legend()
plt.savefig('pdf/criteria_cylinders.pdf', bbox_inches = 'tight',
    pad_inches = 0, transparent=True)

plt.figure(2)
idx = np.argwhere(np.diff(np.sign(Hamburg - EDmax))).flatten()
plt.plot(K,Hamburg,'k-',label=r'thin-wall $\sigma_{max}$ (Hamburg formula)')
plt.plot(K,EDmax,'b-',label=r'thick-wall $ED_{max}$')
plt.plot(K[idx[1]], Hamburg[idx[1]], 'ro')
word = 'K=' + str("%.3f" % K[idx[1]])
tdi = (K[idx[1]]-1)/2
word1 = r'$\frac{t}{d_i}$=' + str("%.3f" % tdi)
plt.text(K[idx[1]], Hamburg[idx[1]]-.015, word)
plt.text(K[idx[1]], Hamburg[idx[1]]-.030, word1)
plt.ylim((0,0.2))
plt.xlim((1,1.2))
plt.xlabel('K')
plt.ylabel(r'$\frac{p_i}{\sigma_{y}}$')
plt.legend()
plt.savefig('pdf/hamburg_thinwall.pdf', bbox_inches = 'tight',
    pad_inches = 0, transparent=True)


plt.figure(3)
idx = np.argwhere(np.diff(np.sign(Hamburg - emax))).flatten()
plt.plot(K,Hamburg,'k-',label=r'thin-wall $\sigma_{max}$ (Hamburg formula)')
plt.plot(K,emax,'b-',label=r'thick-wall $\varepsilon_{max}$')
plt.plot(K[idx[1]], Hamburg[idx[1]], 'ro')
word = 'K=' + str("%.3f" % K[idx[1]])
tdi = (K[idx[1]]-1)/2
word1 = r'$\frac{t}{d_i}$=' + str("%.3f" % tdi)
plt.text(K[idx[1]], Hamburg[idx[1]]-.015, word)
plt.text(K[idx[1]], Hamburg[idx[1]]-.030, word1)
plt.ylim((0,0.2))
plt.xlim((1,1.2))
plt.xlabel('K')
plt.ylabel(r'$\frac{p_i}{\sigma_{y}}$')
plt.legend()
plt.savefig('pdf/hamburg_thinwall_epsilon.pdf', bbox_inches = 'tight',
    pad_inches = 0, transparent=True)

plt.figure(4)
idx = np.argwhere(np.diff(np.sign(ISO - EDmax))).flatten()
plt.plot(K,ISO,'k-',label=r'thin-wall $\tau_{max}$ (ISO formula)')
plt.plot(K,EDmax,'b-',label=r'thick-wall $ED_{max}$')
plt.plot(K[idx[1]], ISO[idx[1]], 'ro')
word = 'K=' + str("%.3f" % K[idx[1]])
tdi = (K[idx[1]]-1)/2
word1 = r'$\frac{t}{d_i}$=' + str("%.3f" % tdi)
plt.text(K[idx[1]], ISO[idx[1]]-.015, word)
plt.text(K[idx[1]], ISO[idx[1]]-.030, word1)
plt.ylim((0,0.2))
plt.xlim((1,1.2))
plt.xlabel('K')
plt.ylabel(r'$\frac{p_i}{\sigma_{y}}$')
plt.legend()
plt.savefig('pdf/ctresca_thinwall.pdf', bbox_inches = 'tight',
    pad_inches = 0, transparent=True)


plt.figure(5)
plt.plot(K,FP,'k-',label=r'$\tau_{max}$ (Fully plastic)')
plt.plot(K,EDmax,'k-.',label=r'$ED_{max}$')
plt.plot(K,tmax,'k:',label=r'$\tau_{max}$')
plt.ylim((0,1.6))
# plt.xlim((1,1.2))
plt.xlabel('K')
plt.ylabel(r'$\frac{p_i}{\sigma_{y}}$')
plt.legend()
plt.savefig('pdf/FP_thickwall.pdf', bbox_inches = 'tight',
    pad_inches = 0, transparent=True)