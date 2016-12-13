# -*- coding: utf-8 -*-
"""
Created on Thu Nov 24 15:28:21 2016

@author: veronica
"""

import numpy as np
import math

def compute_clusters_and_silhouettes(all_dictionaries, all_coefficients):
    """
    Given a list of dictionaries and a list of coefficients that fulfill the requirements specified 
    in the function clustering(all_dictionaries, all_coefficients) it returns the ordered dictionaries 
    and coefficients, the computed dictionary centroids of the clusters and the relative silhouettes
    """
    all_dictionaries, all_coefficients, centroids, clusters = clustering(all_dictionaries, all_coefficients)    
    #silhouettes = np.zeros(len(clusters))    
    #if(len(all_dictionaries) != 1):    
    silhouettes = silhouette(clusters)
    
    sums = np.zeros_like(all_coefficients[0])
    for i in range(0, len(all_coefficients)):
        sums = sums + all_coefficients[i]
    mean_coeffs = sums/len(all_coefficients)
    
    
    return all_dictionaries, all_coefficients, centroids, mean_coeffs, np.mean(silhouettes)



def clustering(all_dictionaries, all_coefficients):
    """
    Function that computes the clusters of all the dictionaries and all the coefficients. 
    It assumes that the dictionaries and the coefficients, except the last one, were already
    clustered and ordered in such way that all the dictionaries contained in the list 
    all_dictionaries have an internal order for example let D_1 = all_dictionaries[0]
    and D_2 = all_dictionaries[1] then if we access the atom in the first position D_1[:,0]
    and D_2[:,0] they both belong to the same cluster. 
    This is true for all the dictionaries and all the atoms. 
    
    It returns the list of all_dictionaries and all_coefficients with the last one ordered
    w.r.t. the previous clusters and the new centroids (of the dictionaries).
    """
    #take the last element (not clustered yet)
    iterations = len(all_dictionaries)
    last   = iterations - 1
    D_last = all_dictionaries[last]
    C_last = all_coefficients[last]
    
    if(iterations == 1):
        #we only have one dictionary and the centroids are the atoms themeselves
        return all_dictionaries, all_coefficients, D_last, D_last
    
    
    #take dimensions
    number_of_features = D_last.shape[0]
    number_of_atoms    = D_last.shape[1]  
    
    #recompute centroids and clusters
    clusters = []
    sums = np.zeros((number_of_features, number_of_atoms))
    for k in range(0, number_of_atoms):
        cluster = []        
        for d in range(0, last):
            sums[:, k] =  sums[:, k] + all_dictionaries[d][:,k]
            cluster.append(all_dictionaries[d][:,k])
        clusters.append(cluster)
    
    centroids = sums/last
    #find the nearest element in D_last to each centroids
    #do not consider cases where an atoms it's the most similar to two or more centroids
    ordered_D  = np.zeros(number_of_atoms)
    used_atoms = []
    for c in range(0,number_of_atoms):
        distances    = compute_distances(centroids[:,c], D_last, used_atoms)
        index      = np.argmax(distances)
        ordered_D[c] = index
        centroids[:,c] = (sums[:,c] + D_last[:, index])/iterations
        clusters[c].append(D_last[:, index])
        used_atoms.append(index)
        
        
    #order D_last and H_last w.r.t. ordered_D
    new_D = np.zeros_like(D_last)
    new_C = np.zeros_like(C_last)
    for i in range(0, number_of_atoms):
        new_D[:,i] = D_last[:,ordered_D[i]]
        new_C[i,:] = C_last[ordered_D[i], :]
    
    all_dictionaries[last] = new_D
    all_coefficients[last] = new_C
    
    return all_dictionaries, all_coefficients, centroids, clusters
    


def silhouette(clusters):
    """ 
    Computes the average silhouette for each cluster where a cluster
    is identified by the same row over all the dictionaries in the list
    all_dictionaries    
    """
    number_of_clusters = len(clusters)
    #final silhouette for each cluster
    silhouettes = np.zeros(number_of_clusters)
    
    #for each cluster 
    for c in range(0, number_of_clusters):
        #take current cluster
        A = clusters[c]
        others = [x for i,x in enumerate(clusters) if i!=c]
        silhouettes[c] = single_cluster_silhouette(A, others)
    return silhouettes

            
            
 
def single_cluster_silhouette(A, others):
    single_point_silhouette = np.zeros(len(A))
    for i in range(0, len(A)):
        a_i = average_distance(A[i], A)
        d_i = np.zeros(len(others))        
        for c in range(0, len(others)):
            d_i[c] = average_distance(A[i], others[c])
        b_i = np.max(d_i)
        
        if(a_i > b_i):
            single_point_silhouette[i] = 1 - b_i/a_i
        elif(a_i == b_i):
            single_point_silhouette[i] = 0
        else:
            single_point_silhouette[i] = a_i/b_i - 1
    return np.mean(single_point_silhouette)
        
    
def compute_distances(centroid, D, used):
    """
    It computes the distances from the point centroid to all the columns in the
    matrix D. If the index of a columns is in the list used that distance with 
    centroid is set to infinity.    
    """
    
    distances = np.zeros(D.shape[1])
    for atom in range(0, D.shape[1]):
        if(atom in used):
            distances[atom] = -float("inf")
        else:            
            centroid = centroid / np.max(centroid)
            vector   = D[:,atom]/ np.max(D[:,atom])
            distances[atom] = (np.dot(centroid, vector.T))/(math.sqrt(np.sum(centroid**2))*math.sqrt(np.sum(vector**2)))
            #distances[atom ] = np.sum((centroid -  vector)**2)
    return distances
        


def average_distance(point, cluster):
    sum_of_distances = 0
    
    for new_point in cluster:
        sum_of_distances = sum_of_distances + (np.dot(point, new_point))/(math.sqrt(np.sum(point**2))*math.sqrt(np.sum(new_point**2)))
    
    return sum_of_distances/len(cluster)    
            
        
        
        
        
        