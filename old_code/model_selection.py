# -*- coding: utf-8 -*-
"""
Created on Tue Dec  6 14:52:51 2016

@author: veronica
"""

import numpy as np
from decipher_mutational_signatures import *

def find_model(mutational_catalogue, min_N = 2, max_N = 20):
    possible_Ns = np.array(range(min_N, max_N+1))
    stabilities = np.zeros((len(possible_Ns)))
    reconstruction_errors = np.zeros((len(possible_Ns)))
    Ds = []
    Cs = []
    
    for n in possible_Ns:
        print("_________________________________________________________")
        print("Trying with N=", n, "\n\n")        
        D, C, stability = decipher_mutational_signatures(mutational_catalogue, n, plot_probabilities = False)
        Ds.append(D)
        Cs.append(C)        
        stabilities[n-min_N] = stability
        reconstruction_errors[n - min_N] = (np.sum((mutational_catalogue - np.dot(D, C))**2))/np.sum(mutational_catalogue**2) 
    
    return Ds, Cs, stabilities, reconstruction_errors