
import sys

import pickle as pkl
import os
os.environ["MKL_NUM_THREADS"] = "8"
os.environ["NUMEXPR_NUM_THREADS"] = "8"
os.environ["OMP_NUM_THREADS"] = "8"

import numpy as np 

from os import listdir
from os.path import isfile, join

from nmtf.nmtf import SSNMTF_CV
from nmtf.read import get_adjacency, get_adjacency_csv

parent_folder = "/cs/research/bioinf/bionet1/Veronica/NMF/network_inference/networks"
#folders = ["HER2Networks", "LuminalBNetworks", "Stage2Networks",
#           "Stage4Networks", "LuminalANetworks",  "Stage1Networks",
#           "Stage3Networks",  "TripleNegativeNetworks"]
#folders = ["Stage2Networks"]

#for i, fold in enumerate(folders):
#print("Analizing group "+fold+"...")
complete_path = parent_folder #xy+ fold
files = [join(complete_path, f) for f in listdir(complete_path)
             if isfile(join(complete_path, f))]
print("... for a total of %d networks.." % len(files))
print("Getting adjacency matrices..")
graphs = [get_adjacency_csv(f) for f in files]
for g in graphs:
	g[np.where(np.isnan(g))] = 0
	g[np.where(np.isinf(g))] = 0
print("Starting cross_validation")
est = SSNMTF_CV(ks=[10,20,30,50,100,200], mode='dognig', verbose=1, max_iter=50)
est.fit(graphs)

with open("cross_val_res_2.pkl", 'wb') as f:
    pkl.dump(est, f)
    print("Finished %d group"%i)
