# -*- coding: utf-8 -*-
"""
Created on Wed Nov 30 12:23:45 2016

@author: veronica
"""

import data
import numpy as np
import NMF
import matplotlib.pyplot as pp

sigma = 0.5
k = 7
V = data.get_synthetic_data(gaussian_noise = sigma)
sparseness_coefficients = 0.1
sparseness_atoms = 0

W, H = NMF.nmf_sparsness_constraint_hoyer(V, k, sparseness_atoms, sparseness_coefficients)
reconstruction_error =  0.5 * np.sum((V - np.dot(W, H))**2)

print("Reconstruction error", reconstruction_error)
for i in range(0, k):
    pp.figure()
    pp.plot(W[:, i])
    pp.show()