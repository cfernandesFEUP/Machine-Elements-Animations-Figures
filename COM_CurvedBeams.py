import numpy as np
import matplotlib.pyplot as plt

E = 210000 # MPa
dphi = 0.1  # degree
phi = 30    # degree
rn = 25     # mm
y = np.linspace(-10,10,10)


sigma = E*y*dphi/((rn+y)*phi)

plt.figure(1)
plt.plot(y,sigma)
