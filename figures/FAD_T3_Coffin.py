import matplotlib as mpl
import numpy as np
import matplotlib.pyplot as plt

N = np.linspace(1, 1e4, 100)

a = 0.5
b = 1

eP = b/N**a

fig = plt.figure(figsize=(5, 2))
ax = fig.add_subplot(1, 1, 1)
plt.loglog(N, eP, 'k-')
plt.xlim(1, 1e4)
plt.ylim(0.001, 1.)
plt.ylabel(r'$\Delta \epsilon_p$')
plt.xlabel(r'Ciclos')
plt.grid(b=True, which='major', color='gray', linestyle='-', alpha=0.35)
plt.grid(b=True, which='minor', color='gray', linestyle='-', alpha=0.35)
# ax.yaxis.set_minor_formatter(mpl.ticker.StrMethodFormatter('{x:.3f}'))
ax.yaxis.set_major_formatter(mpl.ticker.StrMethodFormatter('{x:.3f}'))
# plt.legend()
plt.savefig('pdf/FAD_T3_Coffin.pdf', bbox_inches='tight',
            pad_inches=0, transparent=True)
