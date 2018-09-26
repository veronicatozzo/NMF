# -*- coding: utf-8 -*-
"""
Created on Wed Nov 30 17:06:00 2016

@author: veronica
"""

import numpy as np
import math
import NMF

dims = np.array([2,3,5,10,50,100,500,1000,5000,10000])
#dims = np.array([96, 80])

_ds   = np.array([0.1, 0.3, 0.5, 0.7, 0.9])
_is   = np.array([0.1, 0.3, 0.5, 0.7, 0.9])

ntests = 50

for dsiter in range(0, len(_ds)):
    for isiter in range(0, len(_is)):
        desidered_sparseness = _ds[dsiter]
        initial_sparseness   = _is[isiter]
        print("Desidered sparseness: ", desidered_sparseness)
        print("Initial sparseness: ", initial_sparseness)
        
        for dimiter in range(0, len(dims)):
            N = dims[dimiter]
            print("Dimension: ", N)
            
            for testcase in range(0, ntests):
                sqrt_N = math.sqrt(N) 
                L1_d = sqrt_N - (sqrt_N -1)*desidered_sparseness
                L1_i = sqrt_N - (sqrt_N -1)*initial_sparseness

                x = np.random.normal(size=(N))
                x = x/ np.linalg.norm(x)                
                x = NMF.project_vector(x, L1_d, 1)
                
                s = np.random.normal(size=(N))
                s = s/ np.linalg.norm(s)                
                s = NMF.project_vector(s, L1_i, 1)
                
                v = NMF.project_vector(s, L1_d, 1)
                
               # print("Norm L1", np.sum(abs(v)))
               # print("Norm L2", np.sum(v**2))
                if( (abs(np.sum(abs(v))) - L1_d) > 1e-08 or (abs(np.sum(v**2))-1) > 1e-08):
                    raise Exception("Error, the norms do not correspond")
                    
                if(np.min(v) < 0):
                    raise Exception("Error, the non-negativity is not respected")
                    
                if(np.linalg.norm(x-s) < (np.linalg.norm(v-s) - 1e-10)):
                    raise Exception("Error, non clostest point")
                
                
                
        