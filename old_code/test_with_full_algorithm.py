# -*- coding: utf-8 -*-
"""
Created on Fri Dec  2 16:43:50 2016

@author: veronica
"""

import numpy as np
import matplotlib.pyplot as pp
import decipher_mutational_signatures
import data

pp.close("all")

k = 7
V = data.get_synthetic_data(gaussian_noise = 0.5)
D, C = decipher_mutational_signatures.decipher_mutational_signatures(V, k)

print("reconstruction error", 0.5 * np.sum((V - np.dot(D, C))**2)/np.sum(V**2) )
for i in range(0, k):
    fig = pp.figure()
    fig.canvas.set_window_title('alexandrov algorithm'+ str(i))   
    pp.plot(D[:, i])
    pp.show()