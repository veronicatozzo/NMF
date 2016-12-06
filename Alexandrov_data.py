# -*- coding: utf-8 -*-
"""
Created on Mon Nov 21 15:16:37 2016

@author: veronica
"""



from scipy import io
import remove_weak_mutations as rm_mut
from decipher_mutational_signatures import *
import data
import matplotlib.pyplot as pp
import model_selection

#file_name = "../alexandrov data/input/21_WTSI_BRCA_whole_genome_substitutions.mat"
#data = io.loadmat(file_name, appendmat=False)
#breast_cancers = data["originalGenomes"]

pp.close("all")
#remove the weak mutations (?)
#res = rm_mut.remove_weak_mutatios(breast_cancers, 0.01)
#mutational_catalogue = res["mutational_catalogue"]
#removes_rows = res["removed_rows"]

mutational_catalogue = data.get_synthetic_data(gaussian_noise = 0.5)

#number_of_signatures_to_try = 20
#for number_of_signatures in range(2, number_of_signatures_to_try):

N = 7 #da 1 a min(K, G) -1
    
#decipher mutational signatures
D, C, stability = model_selection.find_model(mutational_catalogue, min_N = 4, max_N = 7)

    
fig = pp.figure()
pp.plot(possible_Ns, stabilities, label = "Model stability")      
pp.plot(possible_Ns, reconstruction_errors, label = "Reconstruction error")
fig.suptitle('Model selection', fontsize=18)
pp.xlabel('Number of searched atoms', fontsize=16)
pp.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
              ncol=2, mode="expand", borderaxespad=0.)
pp.show() 

#for i in range(0, N):
#    fig = pp.figure()
#    fig.canvas.set_window_title('all algorithm'+ str(i))   
#    pp.plot(D[:, i])
#    pp.show()
    
#ottengo D, C, silhouette e errore per ogni N
#plotto tutte le silhouette e gli errori per fare model selection su N
