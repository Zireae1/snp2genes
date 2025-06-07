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
    corr_input = sys.argv[1]
    snp_loc_path = sys.argv[2]
    stats_outpath = sys.argv[3] 
    coord_outpath = sys.argv[4]
    # Read correlation data
    #corr_file_path = f"data/species/{species}/corr_cdist_snp_genes.tsv"
    corr = pd.read_csv(corr_input, sep='\t', index_col=0)

    # Initialize an empty DataFrame for counts
    counts = pd.DataFrame(columns=["gene", "snp", "dist", "corr"])

    # List all files in the snp_loc directory
    #snp_loc_path = f"data/species/{species}/snp_loc/"
    flist = os.listdir(snp_loc_path)

    for f in flist:
        print(f)
        coord_file_path = os.path.join(snp_loc_path, f)
        coord = pd.read_csv(coord_file_path, sep='\t')

        gene = f.replace(".snp.coords.tsv", "")
        sub_snps = coord['SNP_id'].unique()
        sub_snps = [snp for snp in sub_snps if snp in corr.columns]

        # Use list comprehension to gather results and append them to the DataFrame
        results = [
            {
                "gene": gene,
                "snp": snp,
                "dist": coord.loc[coord['SNP_id'] == snp, 'MinDist2'].mean(),
                "corr": corr.at[gene, snp] if gene in corr.index else None
            }
            for snp in sub_snps
        ]

        counts = counts.append(results, ignore_index=True)
    
    counts.to_csv(coord_outpath, sep='\t') 
    
    ### define bins
    xlim = np.array([0.5, 1e5])
    DISTANCE_BINS = np.unique(np.floor(np.logspace(*np.log10(xlim), num=51)).astype(int))
    LD_BINS = np.linspace(0, 1, num=51)

    ### compute LD stats
    ld_data = pd.DataFrame(dict(ld=(1 - counts['corr']) ** 2, pos_diss=counts['dist']))
    bin_stats = ld_data.groupby(pd.cut(ld_data.pos_diss, bins=DISTANCE_BINS, right=False)).ld.apply(lambda x: pd.Series(dict(ld90=x.quantile(0.9), count=len(x)))).unstack().assign(bin_left=DISTANCE_BINS[:-1])
    bin_stats.to_csv(stats_outpath, sep='\t')
