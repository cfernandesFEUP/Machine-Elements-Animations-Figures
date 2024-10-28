# import tikzplotlib as tik
import numpy as np
import matplotlib.pyplot as plt
n = 150
sigmaR = np.linspace(400,1800,n)

def surface(a,b,c):
    return a * c**b

Ca_polido = np.linspace(1,1,n)
Ca_retificado = surface(1.58,-0.085,sigmaR)
Ca_maquinado = surface(4.51,-0.265,sigmaR)
Ca_laminado = surface(57.7,-0.718,sigmaR)
Ca_forjado = surface(272.,-0.955,sigmaR)

plt.figure(1)
plt.plot(sigmaR,Ca_polido,'k',linewidth=5)
plt.plot(sigmaR,Ca_retificado,'k-')
plt.plot(sigmaR,Ca_maquinado,'k-.')
plt.plot(sigmaR,Ca_laminado,'k--')
plt.plot(sigmaR,Ca_forjado,'k:')
plt.xlabel(r'$\mathtt{\sigma_R}$ / MPa')
plt.ylabel(r'$\mathtt{C_3}$')
plt.xlim([400,1800])
plt.ylim([0,1])
plt.grid(True)
plt.legend(['Polido','Retificado','Maquinado','Laminado','Forjado'],loc=3, prop={'size': 11})
plt.savefig('pdf/FAD_T1_acabamento.pdf', bbox_inches = 'tight',
    pad_inches = 0, transparent=True)
# tik.save('acabamento.tex')