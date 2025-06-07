#!/usr/bin/env python3

import sys
import os
import pandas as pd
import numpy as np
from pathlib import Path
from lib.logging_util import info

if __name__ == "__main__":
    stats_input = sys.argv[1]
    master_input = sys.argv[2]
    stats_outpath = sys.argv[3] 
    
    # Read correlation and dist data
    dist = pd.read_csv(stats_input, sep='\t')
    input_path = Path(stats_input)
    # Assumes the path is like: data/species/{species}/snp_gene_counts.tsv
    species = input_path.parts[-2]  # get the second-to-last directory
    
    # Read master stats
    stats = pd.read_csv(master_input, sep="\t")
    
    species = str(species).strip()
    stats['species'] = stats['species'].astype(str).str.strip()

    dhalf = 10 ** stats.loc[stats['species'] == species, 'snp_dhalf'].values[0]
    sub_dist = dist[dist['dist'] < dhalf]
 
    sub_dist.to_csv(stats_outpath, sep='\t')
