#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plots rejection regions for null of 2 independent N(0,1) variables using Wald, KLK, and Optimus tests 
Assumes a 1-sided test with 95% confidence 
"""

import os 
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import minimize 
from scipy.stats import norm 
from scipy.spatial import ConvexHull
import math

#%% 
# Set the working directory to the parent folder of this script 
abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)

#%%
# The Wald rejection region is given by (hat beta1^2 + hat beta2^2) > 5.99  
wald = plt.Circle((0, 0), 2.45, edgecolor='black', fill=False, hatch='////')

#%%
# The KLK rejection region is given by |(hat beta1 + hat beta2)/2| > 1.41  
klk = plt.Rectangle((3, -5.82), 3.98, 15, angle=45, edgecolor='black', fill=False, hatch='....')

#%% 
# Now calculate the optimus rejection region, which needs to be done numerically 

# Create an empty numpy array with the beta hat values 
beta1 = np.arange(-2.82, 2.82, step=0.01)
beta1 = beta1.repeat(564)

beta2 = np.arange(-2.82, 2.82, step=0.01) 
beta2 = np.tile(beta2, 564)

betas = np.empty((564*564, 2))
betas[:,0] = beta1 
betas[:,1] = beta2 

# Define a function to calculate optimus weights as a function of beta1, beta2. Returns weight on beta1 
# First define an objective function to optimize over when both covariates have the same sign 
# Sign we'll use scipy.optimize.minimize, return negative of the objective for positive t-stats
def optimus_same_sign_objective(w, beta1, beta2):
    numerator = w*beta1 + (1-w)*beta2
    denominator = math.sqrt(w**2 + (1-w)**2)  # This is the standard deviation since we have independent betas with variance 1
    return -abs(numerator/denominator)

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
        res = minimize(optimus_same_sign_objective, 0.5, args=(beta1, beta2), bounds=[(0,1)])
        w = res.x 
    return w

weights = np.apply_along_axis(optimus_weights, 1, betas)

# Add the weights to the betas matrix
betas_w_weights = np.hstack([betas, weights])

# Now determine whether the test would accept or reject for each set of weights, betas 
def accept_or_reject(row):
    beta1 = row[0]
    beta2 = row[1]
    w = row[2]
    t = (w*beta1 +(1-w)*beta2)/(math.sqrt(w**2 + (1-w)**2)) 
    p = norm.sf(abs(t))*2
    if p < 0.05: 
        return 1 
    else: 
        return 0 
    
rejections = np.apply_along_axis(accept_or_reject, 1, betas_w_weights)

# Extract the accepted beta values 
accepted_betas = betas[rejections == 0]
hull = ConvexHull(accepted_betas)
accept_region = np.empty((66, 2))
i = 0
for simplex in hull.vertices:
    accept_region[i, 0] = accepted_betas[simplex, 0] 
    accept_region[i, 1] = accepted_betas[simplex, 1]
    i += 1

# Plot the optimus accept region 
optimus = plt.Polygon(accept_region, edgecolor='black', fill=False, hatch='xxx')

#%% 
# Finally generate the actual plot 

plt.rcParams['text.usetex'] = True  # Renders LaTex 

fig, ax = plt.subplots() 

ax.add_patch(optimus)
ax.add_patch(wald)
ax.add_patch(klk)

plt.axis([-3.5, 3.5, -3.5, 3.5])
ax.set_aspect('equal', adjustable='box')
plt.xlabel(r'$\hat{\beta}_1$')
plt.ylabel(r'$\hat{\beta}_2$', rotation=0)
ax.spines['left'].set_position('zero')
ax.spines['right'].set_color('none')
ax.spines['bottom'].set_position('zero')
ax.spines['top'].set_color('none')
ax.xaxis.set_label_coords(0.98, 0.57)
ax.yaxis.set_label_coords(0.55, 0.97)
ax.legend(['Optimus accepts', 'Wald accepts', 'KLK accepts'], frameon=False, loc='upper right', prop={'size': 6})

plt.savefig('rejection-regions.png', dpi=300)

