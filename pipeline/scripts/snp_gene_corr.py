#!/usr/bin/env python3
import sys
import os
import pandas as pd
import numpy as np
from scipy.spatial.distance import cdist
from scipy import stats

# Set working directory
#os.chdir("/Users/vdubinkina/Desktop/Lab/Rotation/cbarrios/PGeSM/r_pipeline/")

from lib.logging_util import info

# Define a custom function to compute correlation ignoring NaNs
#def nan_corr(u, v):
#    return stats.spearmanr(u, v, nan_policy='omit')[0]

if __name__ == "__main__":
    snp_input = sys.argv[1]
    gene_input = sys.argv[2] 

    corr_outpath= sys.argv[3]

    info("Data load")
    # Read genes data
    #genes = pd.read_csv(f"data/subres_80_data/{species}.derep.hit.tsv", sep="\t", header=0)
    genes = pd.read_csv(gene_input, sep='\t', header=0)

    genes.set_index(genes.columns[0], inplace=True)
    genes.drop(genes.columns[0], axis=1, inplace=True)
    genes[genes > 0] = 1

    # Get list of gene files and process names
    #gene_list = os.listdir("data/snp_loc/")
    #gene_list = [gene.replace(".snp.coords.tsv", "") for gene in gene_list]

    # Read SNPs data
    #snps = pd.read_csv(f"data/{species}.snps.tsv", sep="\t", header=0)
    snps = pd.read_csv(snp_input, sep='\t', header=0)
    
    snps.index = snps.apply(lambda row: f"{row['contig']}_{row['pos']}", axis=1)
    snps.drop(['contig', 'pos'], axis=1, inplace=True)
    #snps.columns = snps.columns.str.replace("GUT_GENOME", "")
    snps = snps.loc[:, snps.columns.isin(genes.columns)]
    snps = snps.sort_index(axis = 1)
    genes = genes.sort_index(axis = 1)

    ### convert data to arrays
    snp_np=snps.to_numpy()
    gene_np=genes.to_numpy()
    
    info("Compute correlations")
    #correlation_distances = cdist(gene_np, snp_np, lambda u, v: nan_corr(u, v))
    
    ### use the same metric as for SNP ld decay
    snp_np[np.isnan(snp_np)] = 0.5
    
    correlation_distances = cdist(gene_np, snp_np, metric="correlation")
 
    info("Save correlations")
    corr=pd.DataFrame(correlation_distances)
    corr.index=genes.index
    corr.columns=snps.index

    corr.to_csv(corr_outpath, sep='\t', index=True, header=True)
