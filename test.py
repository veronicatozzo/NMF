# -*- coding: utf-8 -*-
"""
Created on Fri Nov 11 16:51:13 2016

@author: veronica
"""
import numpy as np
import matplotlib.pyplot as pp
import NMF
import data

pp.close("all")




#obtain reconstruction error for different type of noise
sigmas = np.arange(0.01, 1, 0.01)
reconstruction_errors_sigmas = np.zeros((len(sigmas),3))

simple_NMF = 0
sparse_NMF = 1
NMF_sparseness_constraint = 2

k = 7 #right number of coefficients
sparseness_coefficients = 0.1
sparseness_atoms = 0
for i in range(0, len(sigmas)):
    print("sigma", sigmas[i])
    V = data.get_synthetic_data(gaussian_noise = sigmas[i])
    W, H = NMF.non_negative_matrix_factorization(V, k)
    reconstruction_errors_sigmas[i, simple_NMF] = 0.5 * np.sum((V - np.dot(W, H))**2)
    print("Finished simple NMF")
    W, H = NMF.non_negative_sparse_matrix_factorization(V, k, sparseness_coefficients)
    reconstruction_errors_sigmas[i, sparse_NMF] =  0.5 * np.sum((V - np.dot(W, H))**2)
    print("Finished sparse NMF")    
    #W, H = NMF.nmf_sparsness_constraint_hoyer(V, k, sparseness_atoms, sparseness_coefficients)
    #reconstruction_errors[i, NMF_sparseness_constraint] =  0.5 * np.sum((V - np.dot(W, H))**2)

print("Obtained all data")
fig = pp.figure()
pp.plot(sigmas, reconstruction_errors_sigmas[:, simple_NMF], label = "NMF Lee and Seung 2001")
pp.plot(sigmas, reconstruction_errors_sigmas[:, sparse_NMF], label = "NMF Hoyer 2002")
fig.suptitle('Reconstruction error in function of different gaussian noise - synthetic data', fontsize=18)
pp.xlabel('Sigma', fontsize=16)
pp.ylabel('Reconstruction error - Frobenious norm', fontsize=16)
pp.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
       ncol=2, mode="expand", borderaxespad=0.)
pp.show()
   


number_of_atoms = np.array(range(1,20)) #right number of coefficients
sigma = 0.5

reconstruction_errors_atoms = np.zeros((len(number_of_atoms),3))
for i in range(0, len(number_of_atoms)):
    print("number of atoms", number_of_atoms[i])
    V = data.get_synthetic_data(gaussian_noise = sigma)
    W, H = NMF.non_negative_matrix_factorization(V, number_of_atoms[i])
    reconstruction_errors_sigmas[i, simple_NMF] = 0.5 * np.sum((V - np.dot(W, H))**2)
    print("Finished simple NMF")
    W, H = NMF.non_negative_sparse_matrix_factorization(V, number_of_atoms[i], sparseness_coefficients)
    reconstruction_errors_sigmas[i, sparse_NMF] =  0.5 * np.sum((V - np.dot(W, H))**2)
    print("Finished sparse NMF")    
    #W, H = NMF.nmf_sparsness_constraint_hoyer(V, k, sparseness_atoms, sparseness_coefficients)
    #reconstruction_errors[i, NMF_sparseness_constraint] =  0.5 * np.sum((V - np.dot(W, H))**2)

print("Obtained all data")
fig = pp.figure()
pp.plot(sigmas, reconstruction_errors_sigmas[:, simple_NMF], label = "NMF Lee and Seung 2001")
pp.plot(sigmas, reconstruction_errors_sigmas[:, sparse_NMF], label = "NMF Hoyer 2002")
fig.suptitle('Reconstruction error in function of different number of atoms - synthetic data', fontsize=18)
pp.xlabel('Sigma', fontsize=16)
pp.ylabel('Reconstruction error - Frobenious norm', fontsize=16)
pp.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
       ncol=2, mode="expand", borderaxespad=0.)
pp.show()

#    
#k = 7
#retrieved_atoms, retrieved_coefficients = NMF.nmf_sparsness_constraint_hoyer(V, k, atoms_sparseness=0.1, coefficients_sparseness=0.5)#, _lambda)
##check results
#
#for i in range(0, k):
#    pp.figure()
#    pp.plot(retrieved_atoms[:, i])
#    pp.show()
#    
#    
##check that both results are non negative
#assert (np.min(retrieved_atoms) >= 0  and np.min(retrieved_coefficients) >= 0)
#    
#    
#error   = 
#print("Reconstruction error:", error)