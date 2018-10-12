
from __future__ import print_function, division
import sys, getopt

import pandas as pd
import numpy as np

from nmtf.utils import get_clusters
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


def enrichment_go(clusters, dims, annotations):
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
                genes_to_count += list(np.unique(ann_c[ann_c[4]==a].index))
                p_values.append((i, pval, a))
    return p_values, len(set(genes_to_count))/N
