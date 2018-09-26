# -*- coding: utf-8 -*-
"""
Created on Mon Nov 21 16:40:20 2016

@author: veronica
"""
import numpy as np

def remove_weak_mutatios(mutational_catalogue, percentage):
    
    total_sum = np.sum(mutational_catalogue)
    
    rows_sum  = np.sum(mutational_catalogue, 1)
    sorted_indices = np.argsort(rows_sum)
    
    partial_sum = 0
    for i in range(0,len(rows_sum)):
        partial_sum = partial_sum + rows_sum[sorted_indices[i]]
        if(partial_sum/total_sum >= percentage):
            break
    
    rows_to_remove = np.sort(sorted_indices[0:i-1])
    new_mutational_catalogue = np.delete(mutational_catalogue,rows_to_remove, axis=0)

    return {"mutational_catalogue": new_mutational_catalogue, "removed_rows": rows_to_remove}
    


def add_removed_rows(dictionary, removed_rows):
    rows_dictionary = dictionary.shape[0]
    rows = rows_dictionary + len(removed_rows)
    cols = dictionary.shape[1]
    
    D = np.zeros((rows, cols))
    
    count = 0    
    for i in range(0, rows_dictionary):
        if(count < len(removed_rows) and i==removed_rows[count]):
            count = count+1
        D[i+count,:] = dictionary[i,:]
    return D
    
    
def ordering_for_types(D, types):
    possible_types = ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"]
    mutations = dict.fromkeys(possible_types)
    
    for mut in mutations:
        mutations[mut] = np.where(types == mut)
        

    new_D = np.zeros_like(D)
    
    count = 0
    for mut in possible_types:
        print(mut)
        new_D[count:count+len(mutations[mut][0]), :] =  D[mutations[mut][0],:]
        count = count + len(mutations[mut][0])
        
    return new_D