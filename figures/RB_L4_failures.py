import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rcParams
# rcParams['font.family'] = 'sans-serif'
# rcParams['font.sans-serif'] = ['Fira Sans']
rcParams['font.size'] = 10
# rcParams['legend.fontsize'] = 'medium'
params = {'mathtext.default': 'regular' }  

labels = 'Insufficient lubricant', 'solid contamination', 'liquid contamination',\
    'mounting faults', 'consequential damage', 'unsuitable bearing', 'material and production faults',\
        'aged lubricant', 'unsuitable lubricant'
sizes = [15, 20, 5, 5, 5, 10, 1, 20, 20]
explode = (0, 0, 0, 0.15, 0.3, 0.4, 0.5, 0, 0) 

cmap = plt.cm.tab10
cs = cmap(np.linspace(0., 1., len(sizes)))

fig1, ax1 = plt.subplots()
ax1.pie(sizes, explode=explode, labels=labels, colors=cs, autopct='%1.f%%',
        shadow=True, startangle=70)
ax1.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
fig1.savefig('pdf/bearingfailures.pdf', bbox_inches = 'tight',
    pad_inches = 0.1, transparent=True)
plt.show()

