import pandas as pd
import numpy as np
from scipy.sparse import coo_matrix


def get_adjacency(file):
    data = pd.read_table(file,  sep=' ')
    nodes = data.iloc[:, 0].tolist() + data.iloc[:, 1].tolist()
    nodes = sorted(list(set(nodes)))
    print(len(nodes))
    nodes = [(i,nodes[i]) for i in range(len(nodes))]
    for i in range(len(nodes)):
        data = data.replace(nodes[i][1], nodes[i][0])
    M = coo_matrix((data.iloc[:,2], (data.iloc[:,0],data.iloc[:,1])), shape=(len(nodes), len(nodes))).todense()
    M = (M + M.T)/2
    M += np.eye(M.shape[0])
    return M


def get_adjacency_csv(file):
    data = pd.read_csv(file, index_col=0)
    return data.values


def _read_enrichment_results(single_res, integrated_res):
    groups = [ "HER2", "LuminalA", "LuminalB", "TripleNegative",
                "Stage1", "Stage2", "Stage3", "Stage4"]
    network_type = ["AdaptiveLasso", "Lasso", "aracnea", "aracnem",
                    "c3net", "clr", "genenet", "Genie3", "mrnetb",
                    "mrnet", "wgcna", "integrated" ]
    bars = []
    for n in network_type:
        if n =='integrated':
            continue
        that_type = []
        for i in single_res:
            name = i[0].split('/')[-1].split('.')[-2]
            if name.startswith(n):
                if name[len(n)] == 'b':
                    continue
                that_type.append((name, i[2]))
        final_res = []
        for g in groups:
            for i, t in enumerate(that_type):
                if g in t[0]:
                    final_res.append(t[1])
        bars.append(final_res)

    integrated = []
    for g in groups:
        for i in integrated_res:
            if g in i[0]:
                integrated.append(i[2])
    bars.append(integrated)
    return bars, groups, network_type
