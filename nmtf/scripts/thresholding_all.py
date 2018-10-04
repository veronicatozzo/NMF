import sys

import pickle as pkl
import os
os.environ["MKL_NUM_THREADS"] = "8"
os.environ["NUMEXPR_NUM_THREADS"] = "8"
os.environ["OMP_NUM_THREADS"] = "8"

from os import listdir
from os.path import isfile, join

from nmtf.nmtf import SSNMTF_CV, SSNMTF
from nmtf.read import get_adjacency
from nmtf.nmtf import SSNMTF_CV, SSNMTF
from nmtf.read import get_adjacency
from nmtf.integration import integration_SSNMTF
from nmtf.thresholding import thresholding_generating_graphs

parent_folder = "/cs/research/bioinf/bionet1/Coexpression_Study/BCStages-SubtypesNetworks/"
folders = [ "LuminalBNetworks", "Stage2Networks",
           "Stage4Networks", "LuminalANetworks",  "Stage1Networks",
           "Stage3Networks",  "TripleNegativeNetworks"]

ks = [14, 26, 11, 11, 23, 20, 11]

for i, fold in enumerate(folders):
    print("Analizing group "+fold+"...")
    complete_path = parent_folder + fold
    files = [join(complete_path, f) for f in listdir(complete_path)
             if isfile(join(complete_path, f))]
    print("... for a total of %d networks.." % len(files))
    print("Getting adjacency matrices..")
    graphs = [get_adjacency(f) for f in files]

    est = SSNMTF(k=ks[i], init='svd', verbose=1)
    est.fit(graphs)
    print("fitted NMTF")
    with open("../../"+str(fold)+".pkl", 'wb') as f:
        pkl.dump(est, f)

    integrated = integration_SSNMTF(est.G_, est.S_, mode='mean')
    print("starting thresholding")
    res = thresholding_generating_graphs(integrated, min_v=0.01, max_v=0.99,
                    make_plot=False,
                      ax=None, label='', n_repetitions=10)
    with open("../../"+str(fold)+"_thresholding_results.pkl", 'wb') as f:
        pkl.dump(res, f)
