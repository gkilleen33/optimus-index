#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Power isoquants for null of 2 independent N(0,1) variables using Wald, KLK, and Optimus tests 
Assumes a 2-sided test with 95% confidence and 80% power 
"""

import os 
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import norm 
from scipy.spatial import ConvexHull
import math
from copy import deepcopy

#%% 
# Set the working directory to the parent folder of this script 
abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)

#%%
# The Wald isoquant is given by (hat beta1^2 + hat beta2^2) = 9.65  
wald = plt.Circle((0, 0), 3.106, edgecolor='black', fill=False)

# Add 50% power 
wald_50 = plt.Circle((0, 0), 2.22935, edgecolor='black', fill=False)

# Add 20% power 
wald_20 = plt.Circle((0, 0), 1.319, edgecolor='black', fill=False)

#%%
# The KLK isoquant is given by |(hat beta1 + hat beta2)/2| = 3.96  
klk = plt.Rectangle((3, -6.96), 5.60029, 15, angle=45, edgecolor='blue', fill=False)

# Add 50% power 
klk_50 = plt.Rectangle((3, -5.77), 3.917, 15, angle=45, edgecolor='blue', fill=False)

# Add 20% power 
klk_20 = plt.Rectangle((4, -5.58), 2.234, 15, angle=45, edgecolor='blue', fill=False)


#%% 
# Now calculate the optimus isoquant, which needs to be done numerically 

# Create an empty numpy array with the beta hat values 
beta1 = np.arange(-2.82, 2.82, step=0.01)
beta1 = beta1.repeat(564)

beta2 = np.arange(-2.82, 2.82, step=0.01) 
beta2 = np.tile(beta2, 564)

betas = np.empty((564*564, 2))
betas[:,0] = beta1 
betas[:,1] = beta2 

def optimus_weights(betas):
    beta1 = betas[0]
    beta2 = betas[1]
    # If signs don't match we only want to weight one variable to maximize the t-stat
    if beta1 < 0 and beta2 >= 0: 
        if abs(beta1) > abs(beta2):
            w = 1 
        else:
            w = 0
    elif beta1 >= 0 and beta2 < 0:
        if abs(beta1) > abs(beta2):
            w = 1 
        else:
            w = 0
    else:
        w = beta1/(beta1 + beta2)  # Closed form solution to weights
    return w

weights = np.apply_along_axis(optimus_weights, 1, betas)
weights = np.reshape(weights, (len(weights), 1))

# Add the weights to the betas matrix
betas_w_weights = np.hstack([betas, weights])

# Now determine whether the power is (approximately) 0.8 for each set of betas 
def power(row):
    beta1 = row[0]
    beta2 = row[1]
    w = row[2]
    mean = w*beta1 + (1-w)*beta2 
    variance = w**2 + (1-w)**2  # Since beta1 and beta2 are independent with variance 1
    sd = math.sqrt(variance)  
    crit_val = norm.isf(0.025, scale=sd)
    power = 1 - norm.cdf((crit_val - abs(mean))/sd)
    return power
    
power_vector = np.apply_along_axis(power, 1, betas_w_weights)

indicators = np.where(abs(power_vector - 0.8) < .0025, 1, 0)

# Extract the accepted beta values 
isoquat = betas[indicators == 1]
hull = ConvexHull(isoquat)
isoquat_region = np.empty((len(hull.vertices), 2))
i = 0
for simplex in hull.vertices:
    isoquat_region[i, 0] = isoquat[simplex, 0] 
    isoquat_region[i, 1] = isoquat[simplex, 1]
    i += 1

optimus = plt.Polygon(isoquat_region, edgecolor='red', fill=False)


# 50% power 
indicators = np.where(abs(power_vector - 0.5) < .0025, 1, 0)

# Extract the accepted beta values 
isoquat = betas[indicators == 1]
hull = ConvexHull(isoquat)
isoquat_region_50 = np.empty((len(hull.vertices), 2))
i = 0
for simplex in hull.vertices:
    isoquat_region_50[i, 0] = isoquat[simplex, 0] 
    isoquat_region_50[i, 1] = isoquat[simplex, 1]
    i += 1

# Plot the optimus accept region 
optimus_50 = plt.Polygon(isoquat_region_50, edgecolor='red', fill=False)


# 20% power 
indicators = np.where(abs(power_vector - 0.2) < .0025, 1, 0)

# Extract the accepted beta values 
isoquat = betas[indicators == 1]
hull = ConvexHull(isoquat)
isoquat_region_20 = np.empty((len(hull.vertices), 2))
i = 0
for simplex in hull.vertices:
    isoquat_region_20[i, 0] = isoquat[simplex, 0] 
    isoquat_region_20[i, 1] = isoquat[simplex, 1]
    i += 1

# Plot the optimus accept region 
optimus_20 = plt.Polygon(isoquat_region_20, edgecolor='red', fill=False)

#%% 
# Finally generate the actual plots
optimus_1 = deepcopy(optimus)
wald_1 = deepcopy(wald)
klk_1 = deepcopy(klk)
optimus_2 = deepcopy(optimus)
wald_2 = deepcopy(wald)
klk_2 = deepcopy(klk)

plt.rcParams['text.usetex'] = True  # Renders LaTex 

fig, ax = plt.subplots() 

ax.add_patch(optimus_1)
ax.add_patch(wald_1)

plt.axis([-4, 4, -4, 4])
ax.set_aspect('equal', adjustable='box')
plt.xlabel(r'$\beta_1$')
plt.ylabel(r'$\beta_2$', rotation=0)
ax.spines['left'].set_position('zero')
ax.spines['right'].set_color('none')
ax.spines['bottom'].set_position('zero')
ax.spines['top'].set_color('none')
ax.xaxis.set_label_coords(0.98, 0.57)
ax.yaxis.set_label_coords(0.55, 0.97)
ax.legend(['Optimus 80\% power', 'Wald 80\% power'], frameon=False, loc='upper right', prop={'size': 6})

plt.savefig('power-isoquant-optimus-wald.png', dpi=300)

fig, ax = plt.subplots() 
ax.add_patch(optimus_2)
ax.add_patch(klk_1)

plt.axis([-4, 4, -4, 4])
ax.set_aspect('equal', adjustable='box')
plt.xlabel(r'$\beta_1$')
plt.ylabel(r'$\beta_2$', rotation=0)
ax.spines['left'].set_position('zero')
ax.spines['right'].set_color('none')
ax.spines['bottom'].set_position('zero')
ax.spines['top'].set_color('none')
ax.xaxis.set_label_coords(0.98, 0.57)
ax.yaxis.set_label_coords(0.55, 0.97)
ax.legend(['Optimus 80\% power', 'KLK 80\% power'], frameon=False, loc='upper right', prop={'size': 6})

plt.savefig('power-isoquant-optimus-klk.png', dpi=300)


fig, ax = plt.subplots() 
ax.add_patch(wald_2)
ax.add_patch(klk_2)

plt.axis([-4, 4, -4, 4])
ax.set_aspect('equal', adjustable='box')
plt.xlabel(r'$\beta_1$')
plt.ylabel(r'$\beta_2$', rotation=0)
ax.spines['left'].set_position('zero')
ax.spines['right'].set_color('none')
ax.spines['bottom'].set_position('zero')
ax.spines['top'].set_color('none')
ax.xaxis.set_label_coords(0.98, 0.57)
ax.yaxis.set_label_coords(0.55, 0.97)
ax.legend(['Wald 80\% power', 'KLK 80\% power'], frameon=False, loc='upper right', prop={'size': 6})

plt.savefig('power-isoquant-wald-klk.png', dpi=300)


#%% 
# 50% power 
optimus_50_1 = deepcopy(optimus_50)
wald_50_1 = deepcopy(wald_50)
klk_50_1 = deepcopy(klk_50)
optimus_50_2 = deepcopy(optimus_50)
wald_50_2 = deepcopy(wald_50)
klk_50_2 = deepcopy(klk_50)

plt.rcParams['text.usetex'] = True  # Renders LaTex 

fig, ax = plt.subplots() 

ax.add_patch(optimus_50_1)
ax.add_patch(wald_50_1)

plt.axis([-4, 4, -4, 4])
ax.set_aspect('equal', adjustable='box')
plt.xlabel(r'$\beta_1$')
plt.ylabel(r'$\beta_2$', rotation=0)
ax.spines['left'].set_position('zero')
ax.spines['right'].set_color('none')
ax.spines['bottom'].set_position('zero')
ax.spines['top'].set_color('none')
ax.xaxis.set_label_coords(0.98, 0.57)
ax.yaxis.set_label_coords(0.55, 0.97)
ax.legend(['Optimus 50\% power', 'Wald 50\% power'], frameon=False, loc='upper right', prop={'size': 6})

plt.savefig('power-isoquant-optimus-wald-50.png', dpi=300)

fig, ax = plt.subplots() 
ax.add_patch(optimus_50_2)
ax.add_patch(klk_50_1)

plt.axis([-4, 4, -4, 4])
ax.set_aspect('equal', adjustable='box')
plt.xlabel(r'$\beta_1$')
plt.ylabel(r'$\beta_2$', rotation=0)
ax.spines['left'].set_position('zero')
ax.spines['right'].set_color('none')
ax.spines['bottom'].set_position('zero')
ax.spines['top'].set_color('none')
ax.xaxis.set_label_coords(0.98, 0.57)
ax.yaxis.set_label_coords(0.55, 0.97)
ax.legend(['Optimus 50\% power', 'KLK 50\% power'], frameon=False, loc='upper right', prop={'size': 6})

plt.savefig('power-isoquant-optimus-klk-50.png', dpi=300)


fig, ax = plt.subplots() 
ax.add_patch(wald_50_2)
ax.add_patch(klk_50_2)

plt.axis([-4, 4, -4, 4])
ax.set_aspect('equal', adjustable='box')
plt.xlabel(r'$\beta_1$')
plt.ylabel(r'$\beta_2$', rotation=0)
ax.spines['left'].set_position('zero')
ax.spines['right'].set_color('none')
ax.spines['bottom'].set_position('zero')
ax.spines['top'].set_color('none')
ax.xaxis.set_label_coords(0.98, 0.57)
ax.yaxis.set_label_coords(0.55, 0.97)
ax.legend(['Wald 50\% power', 'KLK 50\% power'], frameon=False, loc='upper right', prop={'size': 6})

plt.savefig('power-isoquant-wald-klk-50.png', dpi=300)


#%% 
# 20% power 
optimus_20_1 = deepcopy(optimus_20)
wald_20_1 = deepcopy(wald_20)
klk_20_1 = deepcopy(klk_20)
optimus_20_2 = deepcopy(optimus_20)
wald_20_2 = deepcopy(wald_20)
klk_20_2 = deepcopy(klk_20)

plt.rcParams['text.usetex'] = True  # Renders LaTex 

fig, ax = plt.subplots() 

ax.add_patch(optimus_20_1)
ax.add_patch(wald_20_1)

plt.axis([-4, 4, -4, 4])
ax.set_aspect('equal', adjustable='box')
plt.xlabel(r'$\beta_1$')
plt.ylabel(r'$\beta_2$', rotation=0)
ax.spines['left'].set_position('zero')
ax.spines['right'].set_color('none')
ax.spines['bottom'].set_position('zero')
ax.spines['top'].set_color('none')
ax.xaxis.set_label_coords(0.98, 0.57)
ax.yaxis.set_label_coords(0.55, 0.97)
ax.legend(['Optimus 20\% power', 'Wald 20\% power'], frameon=False, loc='upper right', prop={'size': 6})

plt.savefig('power-isoquant-optimus-wald-20.png', dpi=300)

fig, ax = plt.subplots() 
ax.add_patch(optimus_20_2)
ax.add_patch(klk_20_1)

plt.axis([-4, 4, -4, 4])
ax.set_aspect('equal', adjustable='box')
plt.xlabel(r'$\beta_1$')
plt.ylabel(r'$\beta_2$', rotation=0)
ax.spines['left'].set_position('zero')
ax.spines['right'].set_color('none')
ax.spines['bottom'].set_position('zero')
ax.spines['top'].set_color('none')
ax.xaxis.set_label_coords(0.98, 0.57)
ax.yaxis.set_label_coords(0.55, 0.97)
ax.legend(['Optimus 20\% power', 'KLK 20\% power'], frameon=False, loc='upper right', prop={'size': 6})

plt.savefig('power-isoquant-optimus-klk-20.png', dpi=300)


fig, ax = plt.subplots() 
ax.add_patch(wald_20_2)
ax.add_patch(klk_20_2)

plt.axis([-4, 4, -4, 4])
ax.set_aspect('equal', adjustable='box')
plt.xlabel(r'$\beta_1$')
plt.ylabel(r'$\beta_2$', rotation=0)
ax.spines['left'].set_position('zero')
ax.spines['right'].set_color('none')
ax.spines['bottom'].set_position('zero')
ax.spines['top'].set_color('none')
ax.xaxis.set_label_coords(0.98, 0.57)
ax.yaxis.set_label_coords(0.55, 0.97)
ax.legend(['Wald 20\% power', 'KLK 20\% power'], frameon=False, loc='upper right', prop={'size': 6})

plt.savefig('power-isoquant-wald-klk-20.png', dpi=300)
