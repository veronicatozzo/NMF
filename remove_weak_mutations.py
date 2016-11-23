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
    