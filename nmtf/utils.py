from __future__ import division

import numpy as np
from itertools import combinations


def erdos_renyi(n, m):
    """
    n: int,
        Number of nodes
    m: int,
        Number of edges
    """
    network = np.zeros((n,n))
    comb = list(combinations(np.arange(0, n), 2))

    i = 0
    while i < m:
        v = int(np.random.uniform(low=0, high=len(comb)))
        if network[comb[v]]:
            continue
        network[comb[v]] = 1
        network[comb[v][::-1]] = 1
        i +=1

    return network

def dispersion_coefficient_Kim(X):
    coeff = 0
    for i in range(X.shape[0]):
        for j in range(X.shape[1]):
            coeff += 4*(X[i,j] - 0.5)**2
    return coeff/X.shape[0]**2


def dispersion_coefficients_Dognig(X, k):
    n = X.shape[0]
    non_diag = (np.ones(shape=X.shape) - np.identity(X.shape[0])).astype(bool)
    ravelled = X[np.where(non_diag)]
    eta = np.var(ravelled)
    eta /= ((n/k -1)/(n-1) - ((n/k -1)/(n-1))**2)

    aux = (X - 1/k)**2
    #print(aux)
    aux = aux - np.diag(np.diag(aux))
    v = np.sum(aux)
    v /= n*(n-1)*(1/k - 1/k**2)
    return eta, v

def connectivity_matrix(X):
    indices = np.argmax(X, axis=1)
    C = np.zeros_like(X)
    for r, i in enumerate(indices):
        C[r,i] = 1
    return C.dot(C.T)
