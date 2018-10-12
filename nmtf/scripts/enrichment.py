

from __future__ import print_function
import sys, getopt
import os
os.environ["MKL_NUM_THREADS"] = "8"
os.environ["NUMEXPR_NUM_THREADS"] = "8"
os.environ["OMP_NUM_THREADS"] = "8"

import pickle as pkl
import pandas as pd
import numpy as np

from os import listdir
from os.path import isfile, join
from nmtf.utils import get_clusters
from scipy.special import binom
from scipy.stats import hypergeom
from nmtf.enrichment import enrichment_go, enrichment_kegg


go_annotations = pd.read_table("../../data/goa_human.gaf", sep='\t',
                               skiprows=30, header=None)
go_annotations = go_annotations.set_index(2)
go_annotations = go_annotations[(go_annotations[6] == 'EXP') |
                                (go_annotations[6] == 'IDA') |
                                (go_annotations[6] == 'IPI') |
                                (go_annotations[6] == 'IMP')]
go_annotations = go_annotations[go_annotations[8]=='P']
go_annotations.index = [str(s).lower() for s in go_annotations.index]

with open("../../data/genes_list_networks.pkl", 'rb') as f:
    genes = pkl.load(f)

intersection = list(set(genes).intersection(set(go_annotations.index)))
annotations = go_annotations.loc[intersection]


folder = "../../results_single2/"
files = [join(folder, f) for f in listdir(folder)
         if isfile(join(folder, f))]

results = []
for file in files:
    with open(file, 'rb') as f:
    u = pkl._Unpickler(f)
    u.encoding = 'latin1'
    p = u.load()
    clusters = get_clusters(p.G_)
    clusters_dim = []
    list_clusters = []
    for c in np.unique(clusters):
        indices_c = np.argwhere(clusters==c)
        genes_c = np.array(genes)[indices_c]
        list_clusters.append(genes_c)
        clusters_dim.append(len(genes_c))

    p_values, percentage = enrichment(list_clusters, clusters_dim, annotations)
    print("Done file "+file+" percentage "+str(percentage))
    results.append((file, p_values, percentage))

with open("../../results_enrichment.pkl", 'wb') as f:
    pkl.dump(results, f)
