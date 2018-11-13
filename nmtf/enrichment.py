
from __future__ import print_function, division
import sys, getopt

import pandas as pd
import numpy as np

from scipy.stats import hypergeom

def enrichment_go5(clusters, dims, annotations):
    N = len(np.unique(annotations.index))
    bonferroni_correction = len(np.unique(annotations['GOterm']))*len(dims)
   # print(len(np.unique(annotations['GOterm'])))
    p_values = []
    genes_to_count = []
    enriched_clusters = 0
    for i, c in enumerate(clusters):
        enriched=False
        genes_at_least_one = list(set(annotations.index).intersection(
                                  set(list(c.ravel()))))
        n = len(genes_at_least_one)
        ann_c = annotations.loc[genes_at_least_one]
        for a in np.unique(ann_c['GOterm']):
            genes_annotated = annotations[annotations['GOterm']==a]
            K = genes_annotated.shape[0]
            k = ann_c[ann_c['GOterm']==a].shape[0]
            #1 -  np.sum([binom(K, i)*binom(N-K, n-i)/binom(N, n) for i in range(k)])#
            pval = hypergeom.sf(k-1, N, K, n)
            #print(N, n, K, k)
            if pval < (0.05/bonferroni_correction):
                genes_to_count += list(np.unique(ann_c[ann_c['GOterm']==a].index))
                enriched=True
		 #  print(N, n, K, k)
               # print(pval)
                p_values.append((i, pval, a))
        if enriched:
	        enriched_clusters +=1
    return p_values, len(set(genes_to_count))/N, enriched_clusters/len(clusters)


def enrichment_kegg(clusters, dims, pathways):
    N = len(np.unique(pathways.index))
    bonferroni_correction = len(np.unique(pathways['pathway'])*len(dims))

    p_values = []
    genes_to_count = []
    enriched_clusters = 0
    for i, c in enumerate(clusters):
        enriched=False
        genes_at_least_one = list(set(pathways.index).intersection(
                                  set(list(c.ravel()))))
        n = len(genes_at_least_one)
        ann_c = pathways.loc[genes_at_least_one]
        for p in np.unique(pathways['pathway']):
            genes_annotated = pathways[pathways['pathway'] == p]
            K = genes_annotated.shape[0]
            k = ann_c[ann_c['pathway']==p].shape[0]
            pval = hypergeom.sf(k-1, N, K, n)
        #    print(k, N, K, n)
         #   print((i, pval, p))
            if pval < (0.05/bonferroni_correction):
                genes_to_count += list(np.unique(ann_c[ann_c['pathway']==p].index))
                p_values.append((i, pval, p))
                enriched=True
        if enriched:
	        enriched_clusters+=1
    return p_values, len(set(genes_to_count))/N, enriched_clusters/len(clusters)


def enrichment_go(clusters, dims, annotations):
    N = len(np.unique(annotations.index))
    bonferroni_correction = len(np.unique(annotations[4]))*len(dims)

    p_values = []
    genes_to_count = [] 
    enriched_clusters = 0
    for i, c in enumerate(clusters):
        enriched=False        
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
                enriched=True
        if enriched:
            enriched_clusters+=1
    return p_values, len(set(genes_to_count))/N, enriched_clusters/len(clusters)
