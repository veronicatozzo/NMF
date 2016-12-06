# -*- coding: utf-8 -*-
"""
Created on Fri Nov 11 14:43:30 2016

@author: veronica
"""
import numpy as np
import math
import cmath 

def non_negative_matrix_factorization(V, k, debug = False):
    """
     Function that, given a matrix V (d x n) of non-negative factors and a number
     k, decomposes the matrix V in k atoms and a matrix of coefficients. 
     
     Parameters: 
     
     V: matrix (d x n) of non-negative factors
     
     k: chosen number of atoms
     
     
     Returns: 
     A matrix W (d x k) of k atoms and a matrix H (k x n) of coefficients.  
    """
    if np.min(V) < 0:
        raise ValueError("The argument matrix is not positive")
    
    d = np.shape(V)[0]
    n = np.shape(V)[1]
    epsilon = 0.00001
    
    # devo controllare se k sia minore di d e n or not???

    #create and initialize two matrices with random numbers between 0 and 1
    W = np.random.random((d, k))
    H = np.random.random((k, n))
    
    for iteration in range(1, 100000):    
        W_old = W
        H_old = H
        
        #iterate
        H = H * np.dot(W.T, V)/np.dot(np.dot(W.T,W),H)
        W = W * np.dot(V, H.T)/np.dot(np.dot(W, H), H.T)
        
        #difference
        diff_W = np.linalg.norm(W - W_old, ord = 'fro')
        diff_H = np.linalg.norm(H - H_old, ord = 'fro')
        
        #print("Iteration:", iteration)
        #check convergence
        if diff_W <  epsilon and diff_H < epsilon:
            #print("W difference", diff_W)
            #print("H difference", diff_H)
            break
    
    return  W, H
    
    
def non_negative_sparse_matrix_factorization(V, k, _lambda, debug = False):
    """
     Function that, given a matrix V (d x n) of non-negative factors and a number
     k, decomposes the matrix V in k atoms with unit norm and a sparse matrix of coefficients. 
     
     Parameters: 
     
     V: matrix (d x n) of non-negative factors
     
     k: chosen number of atoms
     
     
     Returns: 
     A matrix W (d x k) of k atoms with each column having unit norm
     and a sparse matrix H (k x n) of coefficients.  
    """
     
    if np.min(V) < 0:
        raise ValueError("The argument matrix is not positive")
    
    d = np.shape(V)[0]
    n = np.shape(V)[1]
    
    epsilon = 0.00001
    mu      = 0.000001 #più grande di 0.00001 non converge   # maybe we should use adaptive steps 
    # devo controllare se k sia minore di d e n or not???

    #create and initialize two matrices with random numbers between 0 and 1
    W = np.random.random((d, k))
    H = np.random.random((k, n))
    
    for iteration in range(1, 100000):    
        
        W_old = W
        H_old = H
        
        #iterate
        W = W - mu* np.dot(np.dot(W, H) - V, H.T)
        W[np.where(W < 0)] = 0
        sums = np.sum(W, axis = 0)
        W = W/sums
        
        H = H * np.dot(W.T, V) / (np.dot(np.dot(W.T, W), H) + _lambda)
        
        #difference
        diff_W = np.sum((W - W_old)**2)
        diff_H = np.sum((H - H_old)**2)
        
        if debug:
            print("Iteration:", iteration)
            print("difference W:", diff_W)
            print("difference_H:", diff_H)
        
        #check convergence
        if diff_W <  epsilon and diff_H < epsilon:
            print("W difference", diff_W)
            print("H difference", diff_H)
            break
    
    return  W, H
    
    
    
    
def nmf_sparsness_constraint_hoyer(V, k, atoms_sparseness, coefficients_sparseness, debug = False):
    if np.min(V) < 0:
        raise ValueError("The argument matrix is not positive")

    if atoms_sparseness < 0 or atoms_sparseness >= 1:
        raise ValueError("The sparseness constant for the atoms must be in the interval (0,1)")
    
    
    if coefficients_sparseness < 0 or coefficients_sparseness >= 1:
        raise ValueError("The sparseness constant for the coefficients must be in the interval (0,1)")
        
    d = np.shape(V)[0]
    n = np.shape(V)[1]
    
    #rescale matrix to avoid overflow and underflow 
    #V = V/np.max(V)

    #create and initialize two matrices with random numbers between 0 and 1
    W = np.random.random((d, k))
    H = np.random.random((k, n))
    
    #project matrices to have specific norm and sparseness
    L1_W = math.sqrt(d) - (math.sqrt(d) - 1)*atoms_sparseness
    W = project_columns(W, L1_W, 1)
    print("finished atoms")
    L1_H = math.sqrt(n) - (math.sqrt(n) - 1)*coefficients_sparseness        
    H = project_columns(H, L1_H, None)
    print("finished coeffcients")
    #compute initial reconstruction error
    reconstruction_error = 0.5 * np.sum((V - np.dot(W, H))**2)
    
    stepsize_W    = 1 
    stepsize_H    = 1
    
    iteration = 0
    while True:
        
        #W_old = W
        #H_old = H
        
        #update H, matrix of coefficients
        dH        = np.dot(W.T, (np.dot(W, H)-V))
        old_error = reconstruction_error        

        while True:
            Hnew = H - stepsize_H * dH
            #Hnew = H * np.dot(W.T, V)/np.dot(np.dot(W.T,W),H)
            Hnew = project_columns(H, L1_H, None)
            
            new_error = 0.5 * np.sum((V - np.dot(W, Hnew))**2)
            if new_error < old_error:
                break
            
            stepsize_H = stepsize_H/2
            if stepsize_H < 1e-200:
                if debug:
                    print("Error", reconstruction_error)
                    np.savetxt("atoms_matrix.csv", W, delimiter=",")
                    np.savetxt("coefficients_matrix.csv", H, delimiter=",")
                return W, H#algorithm has converged
        
        stepsize_H = stepsize_H*1.2 #increase slightly the stepsize
        H = Hnew
        print("finished H")
        
        
        #update W, matrix of coefficients
        dW        = np.dot((np.dot(W, H) -  V), H.T)
        old_error = reconstruction_error        

        while True:
            Wnew = W - stepsize_W * dW
            #Wnew = W = W * np.dot(V, H.T)/np.dot(np.dot(W, H), H.T)
            Wnew = project_columns(W, L1_W, 1)
            
            new_error = 0.5 * np.sum((V - np.dot(Wnew, H))**2)
            if new_error < old_error:
                break
            
            stepsize_W = stepsize_W/2
            if stepsize_W < 1e-200:  
                if debug:
                    print("Error", reconstruction_error)
                    np.savetxt("atoms_matrix.csv", W, delimiter=",")
                    np.savetxt("coefficients_matrix.csv", H, delimiter=",")
                return W, H #algorithm has converged
        
        stepsize_W = stepsize_W*1.2 #increase slightly the stepsize
        W = Wnew
        
        
        reconstruction_error = 0.5 * np.sum((V - np.dot(W, H))**2)
        print("Iteration:", iteration)
        print("Reconstruction error", reconstruction_error)
    
    
    


def project_columns(matrix, sparseness, norm):
    
    new_matrix = np.zeros_like(matrix)    
    
    for i in range(0, matrix.shape[1]):
        new_matrix[:,i] = project_vector(matrix[:,i], sparseness, norm)
    
    return new_matrix
    
   
   
def project_rows(matrix, sparseness, norm):
    
    new_matrix = np.zeros_like(matrix)    
    
    for i in range(0, matrix.shape[0]):
        new_matrix[i,:] = project_vector(matrix[i,:], sparseness, norm)
    
    return new_matrix
    
    
    
def project_vector(x, L1, L2):
    
    #if there isn't any assigned norm, use the norm of the vector    
    if(L2 == None):
        L2 = np.sum(x**2)
        
    #print("L2 norm", L2)
    
    N = len(x)
    
    v = x + (L1 - np.sum(x))/N
    zerocoeffs = np.array([])
    
    #for i in range(1, 10000):
   # print("doing it..")
    while(True):  
       midpoint = np.ones((N)) * L1/(N - len(zerocoeffs))
       if(len(zerocoeffs) != 0):
           midpoint[zerocoeffs] = 0

       aux  = v - midpoint   
       a = np.sum(aux**2)
       b = 2* np.dot(aux.T, v)
       c = np.sum(v**2) - L2
       alpha = (-b + cmath.sqrt(complex(b**2 - 4*a*c,0)).real)/(2*a)
       v = alpha * aux + v

#       a = np.sum((v-midpoint)**2)
#       b = 2* np.dot(midpoint,v-midpoint)
#       c = np.sum(v**2) - L2
#       alpha = (-b + cmath.sqrt(complex(b**2 - 4*a*c,0)).real)/(2*a)
#       v = v + alpha * (v - midpoint)
       
       
       if(np.min(v) >= 0):
           #print("final L2 norm", np.sum(v**2))
    #       print("Done")
           return v
    
       zerocoeffs = np.unique(np.append(zerocoeffs, np.asarray(np.where(v < 0)))).astype(int)
       #print(len(zerocoeffs))
       if(len(zerocoeffs) != 0):
           v[zerocoeffs] = 0
       
       constant = (L1 - np.sum(v))/(N - len(zerocoeffs))
       v = v + constant
       if(len(zerocoeffs) != 0):
           v[zerocoeffs] = 0
       
    #return abs(v)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
  