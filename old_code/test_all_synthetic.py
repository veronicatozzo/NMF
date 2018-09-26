# -*- coding: utf-8 -*-
"""
Created on Mon Dec 12 11:57:09 2016

@author: veronica
"""

from decipher_mutational_signatures import *
import data
import matplotlib.pyplot as pp
import model_selection

pp.close("all")
mutational_catalogue = data.get_synthetic_data(gaussian_noise = 0.5)

minimum = 4
maximum = 10
#decipher mutational signatures
Ds, Cs, stability, reconstruction_errors = model_selection.find_model(mutational_catalogue, min_N = minimum, max_N = maximum)

    
fig = pp.figure()
x_axis = np.arange(minimum, maximum + 1)
pp.scatter(x_axis, stability)
pp.plot(x_axis, stability, label = "Model stability")      
pp.scatter(x_axis, reconstruction_errors)
pp.plot(x_axis, reconstruction_errors, label = "Reconstruction error")    
fig.suptitle('Model selection', fontsize=18)
plt.xticks(np.arange(minimum, maximum + 1, 1.0))
aux = np.concatenate((stability, reconstruction_errors))
pp.axis([minimum - 0.2, maximum + 0.2, min(aux) - 0.5, max(aux) + 0.5])
pp.xlabel('Number of searched atoms', fontsize=16)
pp.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=2, mode="expand", borderaxespad=0.)
pp.show() 
