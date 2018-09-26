# -*- coding: utf-8 -*-
"""
Created on Mon Nov 21 15:16:37 2016

@author: veronica
"""



from scipy import io
import remove_weak_mutations as rm_mut
import bootstrap 

file_name = "../alexandrov data/input/21_WTSI_BRCA_whole_genome_substitutions.mat"
data = io.loadmat(file_name, appendmat=False)
breast_cancers = data["originalGenomes"]


#remove the weak mutations (?)
res = rm_mut.remove_weak_mutatios(breast_cancers, 0.01)
mutational_catalogue = res["mutational_catalogue"]
removes_rows = res["removed_rows"]


number_of_signatures_to_try = 20
for number_of_signatures in range(2, number_of_signatures_to_try):

    #random bootstrap with replacement
    sampled_catalogue = bootstrap.boostrap(mutational_catalogue)    
    
    #decipher mutational signatures
    decipher_mutational_signatures(sampled_catalogue, number_of_signatures)
#try to decompose breast_cancers matrix with different number of patterns, it is necessary to repeat the experiment more times in order 
#to not consider the random initialization ---> random sampling????


#cluster the signatures and compute the stability for each different number of signatures


#plot the stability and the reconstruction error



