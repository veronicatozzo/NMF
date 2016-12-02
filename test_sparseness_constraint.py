# -*- coding: utf-8 -*-
"""
Created on Wed Nov 30 12:23:45 2016

@author: veronica
"""

import data
import numpy as np
import NMF
import matplotlib.pyplot as pp


pp.close("all")

sigma = 0.0000001
k = 7
V = data.get_synthetic_data(gaussian_noise = sigma)
sparseness_coefficients = 0.7
sparseness_atoms = 0.4

W, H = NMF.nmf_sparsness_constraint_hoyer(V, k, sparseness_atoms, sparseness_coefficients)
reconstruction_error =  0.5 * np.sum((V - np.dot(W, H))**2)/np.sum(V**2)

print("Reconstruction error", reconstruction_error)
for i in range(0, k):
    fig = pp.figure()
    fig.canvas.set_window_title('Hoyer gradient descent pattern'+ str(i))   
    pp.plot(W[:, i])
    pp.show()
    
#W, H = NMF.non_negative_matrix_factorization(V, k)#.nmf_sparsness_constraint_hoyer(V, k, sparseness_atoms, sparseness_coefficients)
#reconstruction_error =  0.5 * np.sum((V - np.dot(W, H))**2)/np.sum(V**2)
#
#print("Reconstruction error", reconstruction_error)
#for i in range(0, k):
#    fig = pp.figure()
#    fig.canvas.set_window_title('Lee and seung pattern'+ str(i))     
#    pp.plot(W[:, i])
#    pp.show()
#    
#W, H = NMF.non_negative_sparse_matrix_factorization(V,k, sparseness_coefficients)#.non_negative_matrix_factorization(V, k)#.nmf_sparsness_constraint_hoyer(V, k, sparseness_atoms, sparseness_coefficients)
#reconstruction_error =  0.5 * np.sum((V - np.dot(W, H))**2)/np.sum(V**2)
#
#print("Reconstruction error", reconstruction_error)
#for i in range(0, k):
#    fig = pp.figure()
#    fig.canvas.set_window_title('hoyer multiplicative pattern'+ str(i))   
#    pp.plot(W[:, i])
#    pp.show()