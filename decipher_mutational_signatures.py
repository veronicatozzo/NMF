# -*- coding: utf-8 -*-
"""
Created on Thu Nov 24 11:27:23 2016

@author: veronica
"""

import bootstrap
import NMF


def decipher_mutational_signatures(mutational_catalogue, N):
    
    
    epsilon = 0.0001
    difference = 1
    rows = np.shape(mutational_catalogue)[0]
    cols = np.shape(mutational_catalogue)[1]
    
    D = np.random.random(rows, N)
    C = np.random.random(N, cols)
    
    #until convergence
    while(difference > epsilon)
        old_D = D
        
        #random bootstrap with replacement
        sampled_catalogue = bootstrap.boostrap(mutational_catalogue)   
        
    
        #do NMF on sampled genomes
        res = NMF.non_negative_matrix_factorization(sampled_catalogue, N)
        atoms = res["atoms"]
        coefficients = res["coefficients"]
    
        #cluster this and previous solutions to obtain averaging matrix, silouhette and reconstruction error
    
        
        np.linalg.norm(W - W_old, ord = 'fro')
    
    #return averaged D and C with silouhette and reconstruction error
    
    
    
    