# -*- coding: utf-8 -*-
"""
Created on Tue Nov 29 12:41:21 2016

@author: veronica
"""


import numpy as np
from scipy import signal
import random
from scipy.io import loadmat



def get_synthetic_data(gaussian_noise = 1):
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
    V_tilde = V + np.random.normal(0, gaussian_noise, (number_of_features, number_of_samples))
    V_tilde[np.where(V_tilde < 0)] = 0 
    
    return V_tilde
    
    
def get_breast_cancer_data():
    filename = "/home/veronica/Desktop/Progetto Uli/alexandrov data/input/21_WTSI_BRCA_whole_genome_substitutions.mat"
    data = loadmat(filename,  appendmat=False)
    V = data["originalGenomes"]
    return V
    
    