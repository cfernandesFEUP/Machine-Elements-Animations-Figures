# import tikzplotlib as tik
import numpy as np
import matplotlib.pyplot as plt
n = 150
ymax = 1120
xmax = 2100
sigmaR = np.linspace(0,xmax,n)
sigmaf0 = np.linspace(0,ymax,n)

import csv

xf, xa, xl = [], [], []
yf, ya, yl = [], [], []

with open('dados/ferros.csv','r') as csvfile:
    plots = csv.reader(csvfile, delimiter = ';')
    for row in plots:
        xf.append(float(row[0].replace(',', '.')))
        yf.append(float(row[1].replace(',', '.')))

with open('dados/acos.csv','r') as csvfile:
    plots = csv.reader(csvfile, delimiter=';')
    for row in plots:
        xa.append(float(row[0].replace(',', '.')))
        ya.append(float(row[1].replace(',', '.')))
        
with open('dados/acos_liga.csv','r') as csvfile:
    plots = csv.reader(csvfile, delimiter=';')
    for row in plots:
        xl.append(float(row[0].replace(',', '.')))
        yl.append(float(row[1].replace(',', '.')))


plt.figure(figsize=(8,6))
plt.plot(xa,ya, 'o', markeredgecolor='k', markerfacecolor='w', label='Aços-Carbono')
plt.plot(xl,yl, 'ko',  label='Aços-Liga')
plt.plot(xf,yf, 'k+', label='Ferros forjados')
plt.plot(sigmaR,0.5*sigmaR, 'k--')
plt.plot(sigmaR,0.4*sigmaR, 'k-')
plt.plot(sigmaR,0.6*sigmaR, 'k-')
plt.plot([1470, 2100],[735, 735], 'k--')
plt.xticks(np.arange(0, xmax, 140))
plt.yticks(np.arange(0, ymax, 140))
#plt.title('Interesting Graph\nCheck it out')
plt.xlabel(r'$\mathtt{\sigma_R}$ / MPa')
plt.ylabel(r'$\mathtt{\sigma_{f0}}$ / MPa')
plt.xlim([0,xmax])
plt.ylim([0,ymax])
plt.grid(True)
plt.legend(loc=2, prop={'size': 12})
plt.savefig('pdf/FAD_T1_limite_fadiga.pdf', bbox_inches = 'tight',
    pad_inches = 0, transparent=True)
# plt.show()
# tik.save('limite_fadiga.tex')