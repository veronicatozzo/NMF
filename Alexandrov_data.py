# -*- coding: utf-8 -*-
"""
Created on Mon Nov 21 15:16:37 2016

@author: veronica
"""



from scipy import io
import remove_weak_mutations as rm_mut
import decipher_mutational_signatures

file_name = "../alexandrov data/input/21_WTSI_BRCA_whole_genome_substitutions.mat"
data = io.loadmat(file_name, appendmat=False)
breast_cancers = data["originalGenomes"]


#remove the weak mutations (?)
res = rm_mut.remove_weak_mutatios(breast_cancers, 0.01)
mutational_catalogue = res["mutational_catalogue"]
removes_rows = res["removed_rows"]


#number_of_signatures_to_try = 20
#for number_of_signatures in range(2, number_of_signatures_to_try):

N = 4 #da 1 a min(K, G) -1
    
#decipher mutational signatures
decipher_mutational_signatures.decipher_mutational_signatures(mutational_catalogue, N)
#ottengo D, C, silhouette e errore per ogni N
#plotto tutte le silhouette e gli errori per fare model selection su N
