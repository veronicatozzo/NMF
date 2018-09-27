

from __future__ import print_function
import sys, getopt

import pickle as pkl

from os import listdir
from os.path import isfile, join

from nmtf.nmtf import SSNMTF_CV, SSNMTF
from nmtf.read import get_adjacency

def main(argv):
    k = int(argv[1])
    inputfile = argv[3]
    outputfile = argv[5]
    
    files = [join(inputfile, f) for f in listdir(inputfile)
             if isfile(join(inputfile, f))]

    graphs = [get_adjacency(f) for f in files]
    est = SSNMTF(k, init='svd', verbose=1)
    est.fit(graphs)

    with open(outputfile, 'wb') as f:
        pkl.dump(est, f)


if __name__ == "__main__":
    main(sys.argv[1:])
