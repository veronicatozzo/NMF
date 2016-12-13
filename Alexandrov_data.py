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

#pp.close("all")
#breast_cancers, types = data.get_breast_cancer_data()
#res = rm_mut.remove_weak_mutatios(breast_cancers, 0.01)
#mutational_catalogue = res["mutational_catalogue"]
#removed_rows = res["removed_rows"]



#minimum = 4
#maximum = 10
##decipher mutational signatures
#Ds, Cs, stability, reconstruction_errors = model_selection.find_model(mutational_catalogue, min_N = minimum, max_N = maximum)
#
#    
#fig = pp.figure()
#x_axis = np.arange(minimum, maximum + 1)
#pp.scatter(x_axis, stability)
#pp.plot(x_axis, stability, label = "Model stability")      
#fig.suptitle('Model selection', fontsize=18)
#plt.xticks(np.arange(minimum, maximum + 1, 1.0))
#aux = np.concatenate((stability, reconstruction_errors))
#pp.axis([minimum - 0.2, maximum + 0.2, min(stability) - 0.5, max(stability) + 0.5])
#pp.xlabel('Number of searched atoms', fontsize=16)
#pp.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=2, mode="expand", borderaxespad=0.)
#pp.show() 
#
#fig = pp.figure()
#x_axis = np.arange(minimum, maximum + 1)
#pp.scatter(x_axis, reconstruction_errors)
#pp.plot(x_axis, reconstruction_errors, label = "Reconstruction error")    
#fig.suptitle('Model selection', fontsize=18)
#plt.xticks(np.arange(minimum, maximum + 1, 1.0))
#pp.axis([minimum - 0.2, maximum + 0.2, min(reconstruction_errors) - 0.5, max(reconstruction_errors) + 0.5])
#pp.xlabel('Number of searched atoms', fontsize=16)
#pp.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=2, mode="expand", borderaxespad=0.)
#pp.show()
#
#
#
pp.close("all")
D = Ds[5]
D = rm_mut.add_removed_rows(D, removed_rows)

#for i in range(0, 9):
#    fig = pp.figure()
#    fig.canvas.set_window_title('all algorithm'+ str(i))   
#    x = np.asarray(range(0, len(D[:,i])))
#    w = D[:, i]
#    pp.hist(x, bins = len(D[:,i]), weights =w )
#    #pp.plot(D[:,i])
#    pp.show()

new_D = rm_mut.ordering_for_types(D, types)
for i in range(0, 9):
    fig = pp.figure()
    fig.canvas.set_window_title('all algorithm ordered'+ str(i))   
    x = np.asarray(range(0, len(new_D[:,i])))
    w = new_D[:, i]
    pp.hist(x, bins = len(new_D[:,i]), weights =w )
    #pp.plot(D[:,i])
    pp.show()
