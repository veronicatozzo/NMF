from __future__ import print_function
import sys, getopt

import pickle as pkl

from os import listdir
from os.path import isfile, join

from nmtf.nmtf import SSNMTF_CV, SSNMTF
from nmtf.read import get_adjacency
from nmtf.integration import integration_SSNMTF
from nmtf.thresholding import thresholding_generating_graphs

def main(argv):
    k = argv[0]
    f = argv[1]
    outputfile = argv[2]

    files = [join(complete_path, f) for f in listdir(complete_path)
             if isfile(join(complete_path, f))]
    print("... for a total of %d networks.." % len(files))
    print("Getting adjacency matrices..")
    graphs = [get_adjacency(f) for f in files]
    est = SSNMTF(k=k, init='svd', verbose=1)
    est.fit(graphs)
    print("fitted NMTF")
    with open("../../"+str(f).split('/')[-1]+".pkl", 'wb') as f:
        pkl.dump(est, f)

    integrated = integration_SSNMTF(est.G_, est.S_, mode='mean')
    print(starting thresholding)
    res = thresholding_generating_graphs(integrated, min_v=0.01, max_v=0.99,
                    make_plot=False,
                      ax=None, label='', n_repetitions=10)
    with open("../../"+str(f).split('/')[-1]+"_thresholding_results.pkl", 'wb' as f):
        pkl.dump(res, f)

if __name__ == "__main__":
    main(sys.argv[1:])
