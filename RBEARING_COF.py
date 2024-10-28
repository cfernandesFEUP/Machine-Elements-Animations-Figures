# Import necessary libraries
import numpy as np
import matplotlib.pyplot as plt

# Define variables
# Define an array of 200 evenly spaced values between 0 and 1000
n = np.linspace(0, 1000, 200)

# Define a fixed value for niu
niu = 70

# Define a fixed value for dm (as the average of 30 and 60)
dm = (30+60)/2

# Compute the ratio of the lubrication provided by the boundary layer to the total lubrication
phi_bl = 1/(np.exp(2.6e-8*(n*niu)**1.4*dm))

# Define the coefficients of friction for the boundary layer and EHL
mu_bl = 0.12
mu_EHL = 0.05

# Compute the coefficient of friction
COF = mu_bl * phi_bl + (1 - phi_bl)*mu_EHL

# Create a new figure and Plot
plt.figure()
plt.plot(n, COF, 'k')
plt.xlabel('$n$ / rpm')
plt.ylabel('COF')
plt.grid()
plt.savefig('figures/pdf/COF_SKF.pdf')
