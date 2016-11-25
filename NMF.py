# -*- coding: utf-8 -*-
"""
Created on Fri Nov 11 14:43:30 2016

@author: veronica
"""
import numpy as np
import math

def non_negative_matrix_factorization(V, k):
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
        H = H * np.dot(np.transpose(W), V)/np.dot(np.dot(np.transpose(W),W),H)
        W = W * np.dot(V, np.transpose(H))/np.dot(np.dot(W, H), np.transpose(H))
        
        #difference
        diff_W = np.linalg.norm(W - W_old, ord = 'fro')
        diff_H = np.linalg.norm(H - H_old, ord = 'fro')
        
        #print("Iteration:", iteration)
        #check convergence
        if diff_W <  epsilon and diff_H < epsilon:
            print("W difference", diff_W)
            print("H difference", diff_H)
            break
    
    return  W, H
    
    
def non_negative_sparse_matrix_factorization(V, k, _lambda):
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
        W = W - mu* np.dot(np.dot(W, H) - V, np.transpose(H))
        W[np.where(W < 0)] = 0
        sums = np.sum(W, axis = 0)
        W = W/sums
        
        H = H * np.dot(np.transpose(W), V) / (np.dot(np.dot(np.transpose(W), W), H) + _lambda)
        
        #difference
        diff_W = np.sum((W - W_old)**2)
        diff_H = np.sum((H - H_old)**2)
        
        print("Iteration:", iteration)
        print("difference W:", diff_W)
        print("difference_H:", diff_H)
        #check convergence
        if diff_W <  epsilon and diff_H < epsilon:
            print("W difference", diff_W)
            print("H difference", diff_H)
            break
    
    return {"atoms": W, "coefficients":H}
    
    
    
    
def nmf_sparsness_constraint_hoyer(V, k, atoms_sparseness=None, coefficients_sparseness=None, coefficients_norm=None):
    if np.min(V) < 0:
        raise ValueError("The argument matrix is not positive")
    
    d = np.shape(V)[0]
    n = np.shape(V)[1]

    #create and initialize two matrices with random numbers between 0 and 1
    W = np.random.random((d, k))
    H = np.random.random((k, n))
    atoms_sparsity = atoms_sparseness != None
    coefficients_sparsity = coefficients_sparseness != None    
    
    if(atoms_sparsity):
        W = project_columns(W, atoms_sparseness)
    
    if(coefficients_sparsity):
        H = project_rows(H, coefficients_sparseness, coefficients_norm)
    
    
    epsilon = 0.00001
    mu_W    = 0.000001 # maybe we should use adaptive steps
    mu_H    = 0.000001
    
    
    for iteration in range(1, 100000):    
        
        W_old = W
        H_old = H
        
        #iterate
        if(atoms_sparsity):
            W = W - mu_W * np.dot(np.dot(W, H) - V, np.transpose(H))
            W = project_columns(W, atoms_sparseness)
        else:
            W = W * np.dot(V, np.transpose(H)) / np.dot(np.dot(W, H), np.transpose(H))
        
        
        if(coefficients_sparsity):
            H = H - mu_H * np.dot(np.transpose(W),np.dot(W, H) -V)
            H = project_rows(H, coefficients_sparseness, coefficients_norm)
        else:
            H = H * np.dot(np.transpose(W), V) / np.dot(np.transpose(W), np.dot(W, H))
            
            
        #difference
        diff_W = np.sum((W - W_old)**2)
        diff_H = np.sum((H - H_old)**2)
        
        print("Iteration:", iteration)
        print("difference W:", diff_W)
        print("difference_H:", diff_H)
        #check convergence
        if diff_W <  epsilon and diff_H < epsilon:
            print("W difference", diff_W)
            print("H difference", diff_H)
            break
    
    
    


def project_columns(matrix, sparseness):
   return None
   
   
def project_rows(matrix, sparseness, norm):
    return None
    
    
    
def project_vector(x, L1, L2):
    s = x - (L1 - np.sum(x))/len(x)
    Z = []

    while(True):
       midpoint = np.zeros_like(x)
       for i in range(0, len(x)):
           if i not in Z:
               midpoint[i] = L1/(len(x)- len(Z))
       
       aux   = s - midpoint
       aux1  = np.sum(aux^2)
       aux2  = 2*np.dot(np.transpose(aux), s)
       aux3  = sum(s^2) - L2
       alpha = (-aux2+math.sqrt(aux2^2 - 4* aux1 * aux3))/2*aux1
       s = alpha * aux + midpoint
       
       
       if(np.min(s) > 0):
           return s
    
       less_than_zero_indices = np.where(s < 0)
       Z.append(less_than_zero_indices)
       s[less_than_zero_indices] = 0
       
       constant = (np.sum(s) - L1)( len(x) - len(Z))
       mask = np.ones(len(s), np.bool)
       mask[Z] = 0
       s[mask] = s[mask] - constant
       
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
  