# -*- coding: utf-8 -*-
"""
Created on Thu Nov 24 11:27:23 2016

@author: veronica
"""

import bootstrap
import NMF
import numpy as np
import clustering_and_compactness
import matplotlib.pyplot as pp

def decipher_mutational_signatures(mutational_catalogue, N, plot_probabilities = False, sparse = False, _lambda = 0.1):
    
    
    epsilon = 0.01
    difference = 10
    rows = np.shape(mutational_catalogue)[0]
    cols = np.shape(mutational_catalogue)[1]
    
    D = np.random.random((rows, N))
    C = np.random.random((N, cols))
    all_D = []
    all_C = []
    
    sampling_probabilities = []
    sampling_probabilities.append(get_features_probability(mutational_catalogue))    
    
    mean_D = D
    #until convergence
    while(difference > epsilon):
    #for i in range(0,1):
        old_mean_D = mean_D
        
        #random bootstrap with replacement
        sampled_catalogue = bootstrap.boostrap(mutational_catalogue)   
        sampling_probabilities.append(get_features_probability(sampled_catalogue))
    
        #do NMF on sampled genomes
        if(sparse):
            D, C = NMF.non_negative_sparse_matrix_factorization(sampled_catalogue, N, _lambda = _lambda )
        else:
            D, C = NMF.non_negative_matrix_factorization(sampled_catalogue, N)
            
        print("Factorization computed..")
        all_D.append(D)
        all_C.append(C)
        
        #cluster this and previous solutions to obtain averaging matrix, silouhette and reconstruction error
        all_D, all_C, mean_D, mean_C,  stability = clustering_and_compactness.compute_clusters_and_silhouettes(all_D, all_C)
        print("Finished clustering and silhouette")
        
        for i in range(0, mean_D.shape[1]):
            mean_D[:,i] = mean_D[:,i]/np.max(mean_D[:,i])
            
        difference = np.sum((old_mean_D - mean_D)**2)
        print("Difference", difference)
        
    
    print("Finished the iteration, found the dictionary and the coefficients")
    
    if(plot_probabilities):
        fig = pp.figure()
        pp.plot(range(0,rows), sampling_probabilities[0], label = "Initial probabilities")      
        mean = np.zeros_like(sampling_probabilities[0])
        variance_min = np.zeros_like(sampling_probabilities[0])
        variance_max = np.zeros_like(sampling_probabilities[0])
        for i in range(0, rows):
            sum_row = 0
            minimum = float("inf")
            maximum = -float("inf")
            for vector in sampling_probabilities[1:]:
                sum_row = sum_row + vector[i]
                if(vector[i] < minimum):
                    minimum = vector[i]
                if(vector[i] > maximum):
                    maximum = vector[i]
            mean[i] = sum_row/len(sampling_probabilities)
            variance_min[i] = minimum
            variance_max[i] = maximum
        pp.plot(range(0,rows), mean, label = "Mean")
        pp.plot(range(0,rows), variance_min, label = "Minimum values")
        pp.plot(range(0,rows), variance_max, label = "Maximum values")
        #for vector in sampling_probabilities[1:]:
        #    pp.plot(range(0,rows), vector)      
            
        fig.suptitle('Changes in features probabilities due to random sampling', fontsize=18)
        pp.xlabel('Features', fontsize=16)
        pp.ylabel('Probability', fontsize=16)
        pp.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
                  ncol=2, mode="expand", borderaxespad=0.)
        pp.show()  
        
    #return averaged D and C with silouhette and reconstruction error
    return mean_D, mean_C, stability
    
    


def get_features_probability(M):
    return np.sum(M, axis = 1)/ np.sum(M)