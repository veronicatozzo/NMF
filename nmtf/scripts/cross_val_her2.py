
import sys
sys.path.append("/Users/veronica/Desktop/git_repos/NMF/")

import pickle as pkl

from os import listdir
from os.path import isfile, join

from nmtf.nmtf import SSNMTF_CV
from nmtf.read import get_adjacency

HER2_folder = "/cs/research/bioinf/bionet1/Coexpression_Study/BCStages-SubtypesNetworks/HER2Networks/"
files = [join(HER2_folder, f) for f in listdir(HER2_folder)
         if isfile(join(HER2_folder, f))]

graphs = [get_adjacency(f) for f in files]

est = SSNMTF_CV(verbose=1)
est.fit(graphs)

with open("Cross_validated_HER2.pkl", 'wb') as f:
    pkl.dump(est, f)
