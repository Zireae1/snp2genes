#!/usr/bin/env python3

import os
import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform
import sys

# Set working directory
#os.chdir("/Users/vdubinkina/Desktop/Lab/Rotation/cbarrios/PGeSM/clean_code/r_pipeline/")

#species = "100003"

centroids_path=sys.argv[1]
genes_path=sys.argv[2]
annot_path=sys.argv[3]
stats_outpath=sys.argv[4]

# Load gene annotation files
flist = os.listdir(annot_path)

# Load centroids data
centroids = pd.read_csv(centroids_path, sep="\t", usecols=[0, 5])
centroids.columns = ["gene_id", "centroid_80"]

# Load all genes data
all_genes = pd.read_csv(genes_path, sep="\t")
all_genes.set_index("genome_id", inplace=True)
all_genes = (all_genes > 0).astype(int)  # Binarize

gene_list = all_genes.loc[:, all_genes.mean() > 0.1]
centroids = centroids[centroids["centroid_80"].isin(gene_list.columns)]

# Identify core genes
#core_indices = np.where(gene_list.mean() > 0.9)[0]

nets = []

for f in flist:
    print(f)
    f_path = annot_path + '/' + f
    tmp = pd.read_csv(f_path, sep="\t")
    
    # Add centroid information efficiently
    tmp = tmp.merge(centroids, on="gene_id", how="left")
    tmp = tmp.dropna(subset=["centroid_80"])
    
    # Construct adjacency matrix
    n = gene_list.shape[1]
    net = np.zeros((n, n), dtype=int)
    
    # Loop over unique contig IDs
    '''
    for contig in tmp['contig_id'].unique():
        # Filter the subset of tmp for this contig
        sub_tmp = tmp[tmp['contig_id'] == contig]
    
        # Iterate over consecutive rows in sub_tmp
    
        for i in range(len(sub_tmp) - 1):
            g1 = sub_tmp['centroid_80'].iloc[i]
            g2 = sub_tmp['centroid_80'].iloc[i + 1]
            
            if pd.isna(g1) or pd.isna(g2):
                continue  # Skip this pair if either g1 or g2 is NaN
            else:
                # Get the indices of these genes in gene_list columns
                if g1 in gene_list.columns and g2 in gene_list.columns:
                    id1 = gene_list.columns.get_loc(g1)
                    id2 = gene_list.columns.get_loc(g2)
            
                    # Set the corresponding adjacency matrix entries
                    net[id1, id2] = 1
                    net[id2, id1] = 1
    
    nets.append(net[np.ix_(core_indices, core_indices)][np.triu_indices(len(core_indices), k=1)])
'''
    for contig, sub_tmp in tmp.groupby("contig_id"):
        indices = [gene_list.columns.get_loc(g) for g in sub_tmp["centroid_80"].values if g in gene_list.columns]
        for i in range(len(indices) - 1):
            net[indices[i], indices[i + 1]] = 1
            net[indices[i + 1], indices[i]] = 1
    
    nets.append(net[np.triu_indices(n, k=1)])
# Convert nets into a matrix
nets_mat = np.vstack(nets)

distances = squareform(pdist(nets_mat, metric='euclidean')**2)
synt = 1 - (distances / nets_mat.sum(axis=1)[:, None])

# Save output
pd.DataFrame(synt, index=flist, columns=flist).to_csv(stats_outpath, sep="\t", index=False)

