import pandas as pd
import numpy as np
from scipy.sparse import coo_matrix

def get_adjacency(file):
    data = pd.read_table(file,  sep=' ')
    nodes = data.iloc[:, 0].tolist() + data.iloc[:, 1].tolist()
    nodes = sorted(list(set(nodes)))
    nodes = [(i,nodes[i]) for i in range(len(nodes))]
    for i in range(len(nodes)):
        data = data.replace(nodes[i][1], nodes[i][0])
    M = coo_matrix((data.iloc[:,2], (data.iloc[:,0],data.iloc[:,1])), shape=(len(nodes), len(nodes))).todense()
    M = (M + M.T)/2
    M += np.eye(M.shape[0])
    return M
