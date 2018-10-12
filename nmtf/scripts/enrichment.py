

from __future__ import print_function, division
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


def enrichment_kegg(clusters, dims, pathways):
    N = len(np.unique(pathways.index))
    bonferroni_correction = len(np.unique(pathways['pathway'])*len(dims))

    p_values = []
    genes_to_count = []
    for i, c in enumerate(clusters):
        genes_at_least_one = list(set(pathways.index).intersection(
                                  set(list(c.ravel()))))
        n = len(genes_at_least_one)
        ann_c = pathways.loc[genes_at_least_one]
        for p in np.unique(pathways['pathway']):
            genes_annotated = pathways[pathways['pathway'] == p]
            K = genes_annotated.shape[0]
            k = ann_c[ann_c['pathway']==p].shape[0]
            pval = hypergeom.sf(k, N, K, n)
        #    print(k, N, K, n)
         #   print((i, pval, p))
            if pval < (0.05/bonferroni_correction):
                genes_to_count += list(np.unique(ann_c[ann_c['pathway']==p].index))
                p_values.append((i, pval, p))
    return p_values, len(set(genes_to_count))/N



def enrichment(clusters, dims, annotations):
    N = len(np.unique(annotations.index))
    bonferroni_correction = len(np.unique(annotations[4]))*len(dims)
    p_values = []
    genes_to_count = []
    for i, c in enumerate(clusters):
        genes_at_least_one = list(set(annotations.index).intersection(
                                  set(list(c.ravel()))))
        n = len(genes_at_least_one)
        ann_c = annotations.loc[genes_at_least_one]
        for a in np.unique(ann_c[4]):
            genes_annotated = annotations[annotations[4]==a]
            K = genes_annotated.shape[0]
            k = ann_c[ann_c[4]==a].shape[0]
            pval = hypergeom.sf(k-1, N, K, n)

            if pval < (0.05/bonferroni_correction):
                print(i, pval, a)
		genes_to_count += list(np.unique(ann_c[ann_c[4]==a].index))
                p_values.append((i, pval, a))
    print(len(set(genes_to_count)))
    return p_values, len(set(genes_to_count))/N

go=False
kegg=True

with open("../../data/genes_list_networks.pkl", 'rb') as f:
    genes = pkl.load(f)

if go:	
    go_annotations = pd.read_table("../../data/goa_human.gaf", sep='\t',
                               skiprows=30, header=None)
    go_annotations = go_annotations.set_index(2)
    go_annotations = go_annotations[(go_annotations[6] == 'EXP') |
                                (go_annotations[6] == 'IDA') |
                                (go_annotations[6] == 'IPI') |
                                (go_annotations[6] == 'IMP')]
    go_annotations = go_annotations[go_annotations[8]=='P']
    go_annotations.index = [str(s).lower() for s in go_annotations.index]


    intersection = list(set(genes).intersection(set(go_annotations.index)))
    annotations = go_annotations.loc[intersection]

if kegg:
    pathways = pd.read_csv("../../pathways_kegg.csv", index_col=0)

folder = "../../results_single2/" #"../../network_integrated/"
files = [join(folder, f) for f in listdir(folder)
         if isfile(join(folder, f))]

results = []
for file in files:
    if file.split('.')[-1] != 'pkl':
	continue
    with open(file, 'rb') as f:
   	 p = pkl.load(f)
    clusters = get_clusters(p.G_)
    clusters_dim = []
    list_clusters = []
    for c in np.unique(clusters):
        indices_c = np.argwhere(clusters==c)
        genes_c = np.array(genes)[indices_c]
        list_clusters.append(genes_c)
        clusters_dim.append(len(genes_c))
    if go:
        p_values, percentage = enrichment(list_clusters, clusters_dim, annotations)
    else:
        p_values, percentage = enrichment_kegg(list_clusters, clusters_dim, pathways)
    print("Done file "+file+" percentage "+str(percentage))
    results.append((file, p_values, percentage))

with open("../../results_enrichment_kegg.pkl", 'wb') as f:
    pkl.dump(results, f)
