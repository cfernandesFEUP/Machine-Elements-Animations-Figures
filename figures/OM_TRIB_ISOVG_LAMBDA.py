import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
## Velocidade tangencial ######################################################

vt = np.linspace(0.1,250,1000)

Lambda = (2.68863/vt+0.47767)**(-1)


## LEI ASTM ###################################################################
def ASTM(T,m,n):
    return -0.7+10**(10**(m - n*np.log10(T+273.15)))

xdata = [40, 100]
p01 = [9., 3.5]

ISOVG32 = [27.18, 4.29]
p32, pc32 = curve_fit(ASTM, xdata, ISOVG32, p01)
ISOVG46 = [39.36, 5.44]
p46, pc46 = curve_fit(ASTM, xdata, ISOVG46, p01)
ISOVG68 = [58.65, 7.06]
p68, pc68 = curve_fit(ASTM, xdata, ISOVG68, p01)
ISOVG100 = [86.91, 9.25]
p100, pc100 = curve_fit(ASTM, xdata, ISOVG100, p01)
ISOVG150 = [131.43, 12.28]
p150, pc150 = curve_fit(ASTM, xdata, ISOVG150, p01)
ISOVG220 = [194.24, 15.98]
p220, pc220 = curve_fit(ASTM, xdata, ISOVG220, p01)
ISOVG320 = [284.63, 20.61]
p320, pc320 = curve_fit(ASTM, xdata, ISOVG320, p01)
ISOVG460 = [412.08, 26.34]
p460, pc460 = curve_fit(ASTM, xdata, ISOVG460, p01)
ISOVG680 = [613.83, 34.24]
p680, pc680 = curve_fit(ASTM, xdata, ISOVG680, p01)
ISOVG1000 = [909.48, 42.56]
p1000, pc1000 = curve_fit(ASTM, xdata, ISOVG1000, p01)
ISOVG1500 = [1374.93, 52.59]
p1500, pc1500 = curve_fit(ASTM, xdata, ISOVG1500, p01)

Tlub = np.linspace(0,140,1000)
roh = 900/1000
k = 1.0471
s = 0.1348
## Figura #####################################################################
props = dict(boxstyle='round4', facecolor='white', alpha=1)
plt.figure(1)
# plt.title('Curva de 5% de probabilidade de avaria')
plt.loglog(vt, Lambda,'k')
plt.xlabel('Velocidade tangencial / m/s')
plt.ylabel(r'$\Lambda_{5\%}$')
plt.text(15,0.5,r'$\Lambda_{5\%}=\left(\frac{2.68863}{V} + 0.47767\right)^{-1}$',horizontalalignment='left',
verticalalignment='center', bbox=props)
plt.xlim([0.1,1000])
plt.ylim([0.01,10])
plt.grid(b=True, which='minor', color='gray', linestyle='-', lw=0.25)
plt.grid(b=True, which='major', color='gray', linestyle='-', lw=0.5)
plt.savefig('pdf/lambda5.pdf', bbox_inches = 'tight',
    pad_inches = 0.1, transparent=True)
plt.show()


plt.figure(2)
# plt.semilogy(Tlub,ASTM(Tlub,p1500[0],p1500[1])*roh,'k',label='11')
# plt.text(max(Tlub)+1, min(ASTM(Tlub,p1500[0],p1500[1])*roh), '11', fontsize=5, horizontalalignment='center',verticalalignment='center')
# plt.semilogy(Tlub,ASTM(Tlub,p1000[0],p1000[1])*roh,'k',label='10')
# plt.text(max(Tlub)+1, min(ASTM(Tlub,p1000[0],p1000[1])*roh), '10', fontsize=5, horizontalalignment='center',verticalalignment='center')
plt.semilogy(Tlub,ASTM(Tlub,p680[0],p680[1])*roh,'k',label='9')
plt.text(max(Tlub)+1, min(ASTM(Tlub,p680[0],p680[1])*roh), '9', fontsize=5, horizontalalignment='center',verticalalignment='center')
plt.semilogy(Tlub,ASTM(Tlub,p460[0],p460[1])*roh,'k',label='8')
plt.text(max(Tlub)+1, min(ASTM(Tlub,p460[0],p460[1])*roh), '8', fontsize=5, horizontalalignment='center',verticalalignment='center')
plt.semilogy(Tlub,ASTM(Tlub,p320[0],p320[1])*roh,'k',label='7')
plt.text(max(Tlub)+1, min(ASTM(Tlub,p320[0],p320[1])*roh), '7', fontsize=5, horizontalalignment='center',verticalalignment='center')
plt.semilogy(Tlub,ASTM(Tlub,p220[0],p220[1])*roh,'k',label='6')
plt.text(max(Tlub)+1, min(ASTM(Tlub,p220[0],p220[1])*roh), '6', fontsize=5, horizontalalignment='center',verticalalignment='center')
plt.semilogy(Tlub,ASTM(Tlub,p150[0],p150[1])*roh,'k',label='5')
plt.text(max(Tlub)+1, min(ASTM(Tlub,p150[0],p150[1])*roh), '5', fontsize=5, horizontalalignment='center',verticalalignment='center')
plt.semilogy(Tlub,ASTM(Tlub,p100[0],p100[1])*roh,'k',label='4')
plt.text(max(Tlub)+1, min(ASTM(Tlub,p100[0],p100[1])*roh), '4', fontsize=5, horizontalalignment='center',verticalalignment='center')
plt.semilogy(Tlub,ASTM(Tlub,p68[0],p68[1])*roh,'k',label='3')
plt.text(max(Tlub)+1, min(ASTM(Tlub,p68[0],p68[1])*roh), '3', fontsize=5, horizontalalignment='center',verticalalignment='center')
plt.semilogy(Tlub,ASTM(Tlub,p46[0],p46[1])*roh,'k',label='2')
plt.text(max(Tlub)+1, min(ASTM(Tlub,p46[0],p46[1])*roh), '2', fontsize=5, horizontalalignment='center',verticalalignment='center')
plt.semilogy(Tlub,ASTM(Tlub,p32[0],p32[1])*roh,'k',label='1')
plt.text(max(Tlub)+1, min(ASTM(Tlub,p32[0],p32[1])*roh), '1', fontsize=5, horizontalalignment='center',verticalalignment='center')
plt.xlabel('Temperatura / \u00b0C')
plt.ylabel(r'$\eta$ / mPas')
plt.xlim([0,140])
plt.xticks(range(0,140,10))
plt.ylim([1,100000])
plt.grid(b=True, which='minor', color='gray', linestyle='-', lw=0.25)
plt.grid(b=True, which='major', color='gray', linestyle='-', lw=0.5)
plt.legend((['9 - ISO VG 680',\
             '8 - ISO VG 460','7 - ISO VG 320',\
             '6 - ISO VG 220', '5 - ISO VG 150','4 - ISO VG 100', \
                 '3 - ISO VG 68', '2 - ISO VG 46','1 - ISO VG 32']),handletextpad=0.0, handlelength=0,prop={'family':'cursive','weight':'roman','size':'x-small'})
plt.savefig('pdf/ISOVG.pdf', bbox_inches = 'tight',
    pad_inches = 0.1, transparent=True)
plt.show()

plt.figure(3)
# plt.semilogy(Tlub,k*(ASTM(Tlub,p1500[0],p1500[1])*roh)**(1+s),'k',label='11')
# plt.semilogy(Tlub,k*(ASTM(Tlub,p1000[0],p1000[1])*roh)**(1+s),'k',label='10')
plt.semilogy(Tlub,k*(ASTM(Tlub,p680[0],p680[1])*roh)**(1+s),'k',label='9')
plt.semilogy(Tlub,k*(ASTM(Tlub,p460[0],p460[1])*roh)**(1+s),'k',label='8')
plt.semilogy(Tlub,k*(ASTM(Tlub,p320[0],p320[1])*roh)**(1+s),'k',label='7')
plt.semilogy(Tlub,k*(ASTM(Tlub,p220[0],p220[1])*roh)**(1+s),'k',label='6')
plt.semilogy(Tlub,k*(ASTM(Tlub,p150[0],p150[1])*roh)**(1+s),'k',label='5')
plt.semilogy(Tlub,k*(ASTM(Tlub,p100[0],p100[1])*roh)**(1+s),'k',label='4')
plt.semilogy(Tlub,k*(ASTM(Tlub,p68[0],p68[1])*roh)**(1+s),'k',label='3')
plt.semilogy(Tlub,k*(ASTM(Tlub,p46[0],p46[1])*roh)**(1+s),'k',label='2')
plt.semilogy(Tlub,k*(ASTM(Tlub,p32[0],p32[1])*roh)**(1+s),'k',label='1')
# plt.text(max(Tlub)+1, min(k*(ASTM(Tlub,p1500[0],p1500[1])*roh)**(1+s)), '11', fontsize=5, horizontalalignment='center',verticalalignment='center')
# plt.text(max(Tlub)+1, min(k*(ASTM(Tlub,p1000[0],p1000[1])*roh)**(1+s)), '10', fontsize=5, horizontalalignment='center',verticalalignment='center')
plt.text(max(Tlub)+1, min(k*(ASTM(Tlub,p680[0],p680[1])*roh)**(1+s)), '9', fontsize=5, horizontalalignment='center',verticalalignment='center')
plt.text(max(Tlub)+1, min(k*(ASTM(Tlub,p460[0],p460[1])*roh)**(1+s)), '8', fontsize=5, horizontalalignment='center',verticalalignment='center')
plt.text(max(Tlub)+1, min(k*(ASTM(Tlub,p320[0],p320[1])*roh)**(1+s)), '7', fontsize=5, horizontalalignment='center',verticalalignment='center')
plt.text(max(Tlub)+1, min(k*(ASTM(Tlub,p220[0],p220[1])*roh)**(1+s)), '6', fontsize=5, horizontalalignment='center',verticalalignment='center')
plt.text(max(Tlub)+1, min(k*(ASTM(Tlub,p150[0],p150[1])*roh)**(1+s)), '5', fontsize=5, horizontalalignment='center',verticalalignment='center')
plt.text(max(Tlub)+1, min(k*(ASTM(Tlub,p100[0],p100[1])*roh)**(1+s)), '4', fontsize=5, horizontalalignment='center',verticalalignment='center')
plt.text(max(Tlub)+1, min(k*(ASTM(Tlub,p68[0],p68[1])*roh)**(1+s)), '3', fontsize=5, horizontalalignment='center',verticalalignment='center')
plt.text(max(Tlub)+1, min(k*(ASTM(Tlub,p46[0],p46[1])*roh)**(1+s)), '2', fontsize=5, horizontalalignment='center',verticalalignment='center')
plt.text(max(Tlub)+1, min(k*(ASTM(Tlub,p32[0],p32[1])*roh)**(1+s)), '1', fontsize=5, horizontalalignment='center',verticalalignment='center')

plt.xticks(range(0,140,10))
plt.xlabel('Temperatura / \u00b0C')
plt.ylabel('LP / s')
plt.xlim([0,140])
plt.ylim([1,10000])
plt.legend((['9 - ISO VG 680',\
             '8 - ISO VG 460','7 - ISO VG 320',\
             '6 - ISO VG 220', '5 - ISO VG 150','4 - ISO VG 100', \
                 '3 - ISO VG 68', '2 - ISO VG 46','1 - ISO VG 32']),handletextpad=0.0, handlelength=0,prop={'family':'cursive','weight':'roman','size':'x-small'})
plt.grid(b=True, which='minor', color='gray', linestyle='-', lw=0.25)
plt.grid(b=True, which='major', color='gray', linestyle='-', lw=0.5)
plt.savefig('pdf/ISOVGLP.pdf', bbox_inches = 'tight',
    pad_inches = 0.1, transparent=True)
plt.show()

plt.figure(4)
# plt.semilogy(Tlub,ASTM(Tlub,p1500[0],p1500[1])*roh,'k',label='11')
# plt.text(max(Tlub)+1, min(ASTM(Tlub,p1500[0],p1500[1])*roh), '11', fontsize=5, horizontalalignment='center',verticalalignment='center')
# plt.semilogy(Tlub,ASTM(Tlub,p1000[0],p1000[1])*roh,'k',label='10')
# plt.text(max(Tlub)+1, min(ASTM(Tlub,p1000[0],p1000[1])*roh), '10', fontsize=5, horizontalalignment='center',verticalalignment='center')
plt.semilogy(Tlub,ASTM(Tlub,p680[0],p680[1])*roh,'k',label='9')
plt.text(max(Tlub)+1, min(ASTM(Tlub,p680[0],p680[1])*roh), '9', fontsize=5, horizontalalignment='center',verticalalignment='center')
plt.semilogy(Tlub,ASTM(Tlub,p460[0],p460[1])*roh,'k',label='8')
plt.text(max(Tlub)+1, min(ASTM(Tlub,p460[0],p460[1])*roh), '8', fontsize=5, horizontalalignment='center',verticalalignment='center')
plt.semilogy(Tlub,ASTM(Tlub,p320[0],p320[1])*roh,'k',label='7')
plt.text(max(Tlub)+1, min(ASTM(Tlub,p320[0],p320[1])*roh), '7', fontsize=5, horizontalalignment='center',verticalalignment='center')
plt.semilogy(Tlub,ASTM(Tlub,p220[0],p220[1])*roh,'k',label='6')
plt.text(max(Tlub)+1, min(ASTM(Tlub,p220[0],p220[1])*roh), '6', fontsize=5, horizontalalignment='center',verticalalignment='center')
plt.semilogy(Tlub,ASTM(Tlub,p150[0],p150[1])*roh,'k',label='5')
plt.text(max(Tlub)+1, min(ASTM(Tlub,p150[0],p150[1])*roh), '5', fontsize=5, horizontalalignment='center',verticalalignment='center')
plt.semilogy(Tlub,ASTM(Tlub,p100[0],p100[1])*roh,'k',label='4')
plt.text(max(Tlub)+1, min(ASTM(Tlub,p100[0],p100[1])*roh), '4', fontsize=5, horizontalalignment='center',verticalalignment='center')
plt.semilogy(Tlub,ASTM(Tlub,p68[0],p68[1])*roh,'k',label='3')
plt.text(max(Tlub)+1, min(ASTM(Tlub,p68[0],p68[1])*roh), '3', fontsize=5, horizontalalignment='center',verticalalignment='center')
plt.semilogy(Tlub,ASTM(Tlub,p46[0],p46[1])*roh,'k',label='2')
plt.text(max(Tlub)+1, min(ASTM(Tlub,p46[0],p46[1])*roh), '2', fontsize=5, horizontalalignment='center',verticalalignment='center')
plt.semilogy(Tlub,ASTM(Tlub,p32[0],p32[1])*roh,'k',label='1')
plt.text(max(Tlub)+1, min(ASTM(Tlub,p32[0],p32[1])*roh), '1', fontsize=5, horizontalalignment='center',verticalalignment='center')
plt.xlabel('Temperature / \u00b0C')
plt.ylabel(r'$\eta$ / mPas')
plt.xlim([0,140])
plt.xticks(range(0,140,10))
plt.ylim([1,100000])
plt.grid(b=True, which='minor', color='gray', linestyle='-', lw=0.25)
plt.grid(b=True, which='major', color='gray', linestyle='-', lw=0.5)
plt.legend((['9 - ISO VG 680',\
             '8 - ISO VG 460','7 - ISO VG 320',\
             '6 - ISO VG 220', '5 - ISO VG 150','4 - ISO VG 100', \
                 '3 - ISO VG 68', '2 - ISO VG 46','1 - ISO VG 32']),handletextpad=0.0, handlelength=0,prop={'family':'cursive','weight':'roman','size':'x-small'})
plt.savefig('pdf/ISOVG_en.pdf', bbox_inches = 'tight',
    pad_inches = 0.1, transparent=True)
plt.show()

plt.figure(5)
# plt.semilogy(Tlub,k*(ASTM(Tlub,p1500[0],p1500[1])*roh)**(1+s),'k',label='11')
# plt.semilogy(Tlub,k*(ASTM(Tlub,p1000[0],p1000[1])*roh)**(1+s),'k',label='10')
plt.semilogy(Tlub,k*(ASTM(Tlub,p680[0],p680[1])*roh)**(1+s),'k',label='9')
plt.semilogy(Tlub,k*(ASTM(Tlub,p460[0],p460[1])*roh)**(1+s),'k',label='8')
plt.semilogy(Tlub,k*(ASTM(Tlub,p320[0],p320[1])*roh)**(1+s),'k',label='7')
plt.semilogy(Tlub,k*(ASTM(Tlub,p220[0],p220[1])*roh)**(1+s),'k',label='6')
plt.semilogy(Tlub,k*(ASTM(Tlub,p150[0],p150[1])*roh)**(1+s),'k',label='5')
plt.semilogy(Tlub,k*(ASTM(Tlub,p100[0],p100[1])*roh)**(1+s),'k',label='4')
plt.semilogy(Tlub,k*(ASTM(Tlub,p68[0],p68[1])*roh)**(1+s),'k',label='3')
plt.semilogy(Tlub,k*(ASTM(Tlub,p46[0],p46[1])*roh)**(1+s),'k',label='2')
plt.semilogy(Tlub,k*(ASTM(Tlub,p32[0],p32[1])*roh)**(1+s),'k',label='1')
# plt.text(max(Tlub)+1, min(k*(ASTM(Tlub,p1500[0],p1500[1])*roh)**(1+s)), '11', fontsize=5, horizontalalignment='center',verticalalignment='center')
# plt.text(max(Tlub)+1, min(k*(ASTM(Tlub,p1000[0],p1000[1])*roh)**(1+s)), '10', fontsize=5, horizontalalignment='center',verticalalignment='center')
plt.text(max(Tlub)+1, min(k*(ASTM(Tlub,p680[0],p680[1])*roh)**(1+s)), '9', fontsize=5, horizontalalignment='center',verticalalignment='center')
plt.text(max(Tlub)+1, min(k*(ASTM(Tlub,p460[0],p460[1])*roh)**(1+s)), '8', fontsize=5, horizontalalignment='center',verticalalignment='center')
plt.text(max(Tlub)+1, min(k*(ASTM(Tlub,p320[0],p320[1])*roh)**(1+s)), '7', fontsize=5, horizontalalignment='center',verticalalignment='center')
plt.text(max(Tlub)+1, min(k*(ASTM(Tlub,p220[0],p220[1])*roh)**(1+s)), '6', fontsize=5, horizontalalignment='center',verticalalignment='center')
plt.text(max(Tlub)+1, min(k*(ASTM(Tlub,p150[0],p150[1])*roh)**(1+s)), '5', fontsize=5, horizontalalignment='center',verticalalignment='center')
plt.text(max(Tlub)+1, min(k*(ASTM(Tlub,p100[0],p100[1])*roh)**(1+s)), '4', fontsize=5, horizontalalignment='center',verticalalignment='center')
plt.text(max(Tlub)+1, min(k*(ASTM(Tlub,p68[0],p68[1])*roh)**(1+s)), '3', fontsize=5, horizontalalignment='center',verticalalignment='center')
plt.text(max(Tlub)+1, min(k*(ASTM(Tlub,p46[0],p46[1])*roh)**(1+s)), '2', fontsize=5, horizontalalignment='center',verticalalignment='center')
plt.text(max(Tlub)+1, min(k*(ASTM(Tlub,p32[0],p32[1])*roh)**(1+s)), '1', fontsize=5, horizontalalignment='center',verticalalignment='center')

plt.xticks(range(0,140,10))
plt.xlabel('Temperature / \u00b0C')
plt.ylabel('LP / s')
plt.xlim([0,140])
plt.ylim([1,10000])
plt.legend((['9 - ISO VG 680',\
             '8 - ISO VG 460','7 - ISO VG 320',\
             '6 - ISO VG 220', '5 - ISO VG 150','4 - ISO VG 100', \
                 '3 - ISO VG 68', '2 - ISO VG 46','1 - ISO VG 32']),handletextpad=0.0, handlelength=0,prop={'family':'cursive','weight':'roman','size':'x-small'})
plt.grid(b=True, which='minor', color='gray', linestyle='-', lw=0.25)
plt.grid(b=True, which='major', color='gray', linestyle='-', lw=0.5)
plt.savefig('pdf/ISOVGLP_en.pdf', bbox_inches = 'tight',
    pad_inches = 0.1, transparent=True)
plt.show()

plt.figure(6)
# plt.title('Curva de 5% de probabilidade de avaria')
plt.loglog(vt, Lambda,'k')
plt.xlabel('Pitch line velocity / m/s')
plt.ylabel(r'$\Lambda_{5\%}$')
plt.text(15,0.5,r'$\Lambda_{5\%}=\left(\frac{2.68863}{V} + 0.47767\right)^{-1}$',horizontalalignment='left',
verticalalignment='center', bbox=props)
plt.xlim([0.1,1000])
plt.ylim([0.01,10])
plt.grid(b=True, which='minor', color='gray', linestyle='-', lw=0.25)
plt.grid(b=True, which='major', color='gray', linestyle='-', lw=0.5)
plt.savefig('pdf/lambda5_en.pdf', bbox_inches = 'tight',
    pad_inches = 0.1, transparent=True)
plt.show()
