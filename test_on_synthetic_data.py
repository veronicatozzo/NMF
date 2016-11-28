# -*- coding: utf-8 -*-
"""
Created on Fri Nov 11 16:51:13 2016

@author: veronica
"""
import numpy as np
from scipy import signal
import matplotlib.pyplot as pp
import random
import NMF


#generates true patterns with 96 features
number_of_features = 96
number_of_samples = 80
number_of_atoms = 7


atoms = np.empty([number_of_features, number_of_atoms])
atoms[:, 0] = np.transpose(np.concatenate((np.ones([30,1]), np.zeros([66, 1]))))
atoms[:, 1] = np.transpose(np.concatenate((np.zeros([60,1]), np.ones([36, 1]))))
atoms[:, 2] = np.transpose(np.concatenate((np.zeros([24,1]), np.ones([30, 1]), np.zeros([42,1]))))
atoms[:, 3] = signal.gaussian(96, 5)
atoms[:, 4] = np.transpose(np.concatenate((np.zeros([17,1]), np.ones([15, 1]), np.zeros([30,1]), np.ones([24, 1]), np.zeros([10,1]))))
atoms[:, 5] = np.roll(signal.gaussian(96, 5), 30)
atoms[:, 6] = signal.gaussian(96, 8)
atoms[0:50, 6] = 0

sums = np.sum(atoms, axis = 0)
atoms = atoms/sums

#for i in range(0, 7):
# pp.figure()
# pp.plot(atoms[:, i])
# pp.show()
 
 

#create sparse coefficients 
coefficients = np.zeros([number_of_atoms, number_of_samples])
for i in range(0, number_of_samples):
    number_of_nonzero_elements = random.randint(2, 4)
    indices = random.sample(range(0, 7), number_of_nonzero_elements)   
    coeffs = random.sample(range(0, 100), number_of_nonzero_elements)
    coefficients[indices, i] = coeffs
    
#print(coefficients[:, 0:10])

#create matrix
V = np.dot(atoms, coefficients)


#add noise
V_tilde = V + np.random.normal(0, 0.5, (number_of_features, number_of_samples))
V_tilde[np.where(V_tilde < 0)] = 0 

#run decomposition
#res = NMF.non_negative_matrix_factorization(V_tilde, 10)
_lambda = 0.1
k = 11
retrieved_atoms, retrieved_coefficients = NMF.nmf_sparsness_constraint_hoyer(V_tilde, k, atoms_sparseness=0.1, coefficients_sparseness=0.1)#, _lambda)
#check results

for i in range(0, k):
    pp.figure()
    pp.plot(retrieved_atoms[:, i])
    pp.show()
    
    
#check that both results are non negative
assert (np.min(retrieved_atoms) >= 0  and np.min(retrieved_coefficients) >= 0)
    
    
V_tilde = np.dot(retrieved_atoms, retrieved_coefficients)
error   = np.sum((V - V_tilde)**2)
print("Reconstruction error:", error)