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

simple_NMF = 0
sparse_NMF = 1
NMF_sparseness_constraint = 2

sparseness_coefficients = 0.1
sparseness_atoms = 0

number_of_atoms = np.array(range(10,10))#right number of coefficients
sigma = 0.5

reconstruction_errors = np.zeros((len(number_of_atoms),3))
for i in range(0, len(number_of_atoms)):
    print("number of atoms", number_of_atoms[i])
    V = data.get_breast_cancer_data()
    W_1, H_1 = NMF.non_negative_matrix_factorization(V, number_of_atoms[i])
    reconstruction_errors[i, simple_NMF] = 0.5 * np.sum((V - np.dot(W_1, H_1))**2)
    print("Finished simple NMF")
    W_2, H_2 = NMF.non_negative_sparse_matrix_factorization(V, number_of_atoms[i], sparseness_coefficients)
    reconstruction_errors[i, sparse_NMF] =  0.5 * np.sum((V - np.dot(W_2, H_2))**2)
    print("Finished sparse NMF")    
    #W, H = NMF.nmf_sparsness_constraint_hoyer(V, k, sparseness_atoms, sparseness_coefficients)
    #reconstruction_errors[i, NMF_sparseness_constraint] =  0.5 * np.sum((V - np.dot(W, H))**2)

print("Obtained all data")
fig = pp.figure()
pp.plot(number_of_atoms, reconstruction_errors[:, simple_NMF], label = "NMF Lee and Seung 2001")
pp.plot(number_of_atoms, reconstruction_errors[:, sparse_NMF], label = "NMF Hoyer 2002")
fig.suptitle('Reconstruction error in function of different number of atoms - breast cancer data', fontsize=18)
pp.xlabel('Number of atoms', fontsize=16)
pp.ylabel('Reconstruction error - Frobenious norm', fontsize=16)
pp.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
       ncol=2, mode="expand", borderaxespad=0.)
pp.show()