

from __future__ import print_function
import sys, getopt
import os
os.environ["MKL_NUM_THREADS"] = "8"
os.environ["NUMEXPR_NUM_THREADS"] = "8"
os.environ["OMP_NUM_THREADS"] = "8"

import pickle as pkl
import numpy as np 

from os import listdir
from os.path import isfile, join

from nmtf.nmtf import SSNMTF_CV, SSNMTF
from nmtf.read import get_adjacency, get_adjacency_csv


parent_folder = "/cs/research/bioinf/bionet1/Veronica/NMF/network_inference/networks"
#folders = ["HER2Networks", "LuminalBNetworks", "Stage2Networks",
#           "Stage4Networks", "LuminalANetworks",  "Stage1Networks",
#           "Stage3Networks",  "TripleNegativeNetworks"]
#folders = ["Stage2Networks"]

#for i, fold in enumerate(folders):
#print("Analizing group "+fold+"...")
complete_path = parent_folder
files = [join(complete_path, f) for f in listdir(complete_path)
             if isfile(join(complete_path, f))]
print("... for a total of %d networks.." % len(files))
print("Getting adjacency matrices..")
graphs = [get_adjacency_csv(f) for f in files]
print("Starting cross_validation")
for g in graphs:
	g[np.where(np.isnan(g))] = 0
        g[np.where(np.logical_not(np.isfinite(g)))] = 0
est = SSNMTF(k=50, init='svd', verbose=1)
est.fit(graphs)

name = f.split("/")[-1].split('.')[-2]
with open("../../results_network_big.pkl", 'wb') as f:
    pkl.dump(est, f)
