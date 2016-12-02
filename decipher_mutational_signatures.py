# -*- coding: utf-8 -*-
"""
Created on Thu Nov 24 11:27:23 2016

@author: veronica
"""

import bootstrap
import NMF
import numpy as np
import clustering_and_compactness


def decipher_mutational_signatures(mutational_catalogue, N):
    
    
    epsilon = 0.0001
    difference = 1
    rows = np.shape(mutational_catalogue)[0]
    cols = np.shape(mutational_catalogue)[1]
    
    D = np.random.random(rows, N)
    C = np.random.random(N, cols)
    all_D = []
    all_C = []
    
    mean_D = D
    #until convergence
    while(difference > epsilon):
        old_mean_D = mean_D
        
        #random bootstrap with replacement
        sampled_catalogue = bootstrap.boostrap(mutational_catalogue)   
        
    
        #do NMF on sampled genomes
        D, C = NMF.non_negative_matrix_factorization(sampled_catalogue, N)
        all_D.append(D)
        all_C.append(C)
        
        #cluster this and previous solutions to obtain averaging matrix, silouhette and reconstruction error
        all_D, all_C, mean_D, silhouettes = clustering_and_compactness.compute_clusters_and_silhouettes(all_D, all_C)
        
        difference = np.sum((old_mean_D - mean_D)**2)
    
    coefficients = all_C#compute_coefficients_with_sparsity()
    #return averaged D and C with silouhette and reconstruction error
    return mean_D, coefficients
    
    
    