
import sys

import pickle as pkl
import os
os.environ["MKL_NUM_THREADS"] = "8"
os.environ["NUMEXPR_NUM_THREADS"] = "8"
os.environ["OMP_NUM_THREADS"] = "8"

from os import listdir
from os.path import isfile, join

from nmtf.nmtf import SSNMTF_CV
from nmtf.read import get_adjacency

parent_folder = "/cs/research/bioinf/bionet1/Coexpression_Study/BCStages-SubtypesNetworks/"
folders = ["HER2Networks", "LuminalBNetworks", "Stage2Networks",
           "Stage4Networks", "LuminalANetworks",  "Stage1Networks",
           "Stage3Networks",  "TripleNegativeNetworks"]

for i, fold in enumerate(folders):
    print("Analizing group "+fold+"...")
    complete_path = parent_folder + fold
    files = [join(complete_path, f) for f in listdir(complete_path)
             if isfile(join(complete_path, f))]
    print("... for a total of %d networks.." % len(files))
    print("Getting adjacency matrices..")
    graphs = [get_adjacency(f) for f in files]
    print("Starting cross_validation")
    est = SSNMTF_CV(mode='Dognig', verbose=1)
    est.fit(graphs)

    with open("../../results/cross_val_res"+folder+".pkl", 'wb') as f:
        pkl.dump(est, f)
    print("Finished %d group"%i)
