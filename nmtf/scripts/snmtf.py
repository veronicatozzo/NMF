

from __future__ import print_function
import sys, getopt
import os
os.environ["MKL_NUM_THREADS"] = "8"
os.environ["NUMEXPR_NUM_THREADS"] = "8"
os.environ["OMP_NUM_THREADS"] = "8"

import pickle as pkl

from os import listdir
from os.path import isfile, join

from nmtf.nmtf import SSNMTF_CV, SSNMTF
from nmtf.read import get_adjacency


parent_folder = "/cs/research/bioinf/bionet1/Coexpression_Study/BCStages-SubtypesNetworks/"
folders = ["HER2Networks", "LuminalBNetworks", "Stage2Networks",
           "Stage4Networks", "LuminalANetworks",  "Stage1Networks",
           "Stage3Networks",  "TripleNegativeNetworks"]
#folders =  [ "Stage1Networks", "Stage2Networks", "Stage3Networks", "Stage4Networks"]
ks = [14, 14, 26, 11, 11, 23, 20, 11]
for i, fold in enumerate(folders):
    print("Analizing group "+fold+"...")
    complete_path = parent_folder + fold
    files = [join(complete_path, f) for f in listdir(complete_path)
             if isfile(join(complete_path, f))]
    for f in files:
        graphs = [get_adjacency(f)]
        est = SSNMTF(ks[i], init='svd', verbose=0)
        est.fit(graphs)

        name = f.split("/")[-1].split('.')[-2]
        with open("../../results_single2/"+name+fold+".pkl", 'wb') as f:
            pkl.dump(est, f)
