# -*- coding: utf-8 -*-
"""
Created on Mon Nov 21 17:57:56 2016

@author: veronica
"""
import numpy as np


def boostrap(catalogue):
    selected = np.random.randint(0, high = catalogue.shape[1], size = (catalogue.shape[1]))
    return catalogue[:,selected] 