#!/usr/bin/env python3

import pandas as pd
import subprocess
import pathlib
import sys
from lib.logging_util import info

if __name__ == "__main__":
    gene_inpath = sys.argv[1]
    min_genes = float(sys.argv[2])
    max_genes = float(sys.argv[3])
    gene_outpath = sys.argv[4]
   
    info("Loading data.")
    df_gene = pd.read_csv(gene_inpath,sep='\t', index_col=0, header=0)
    df_gene = df_gene.T
    info(df_gene.shape)
    
    # binarize gene data!!!
    df_gene[df_gene>0]=1
    
    ### delete core/rare genes
    rare_genes = df_gene.mean(axis=1) < min_genes
    rare_genes = rare_genes[rare_genes == True].index

    # remove rare genes from df_gene
    df_gene = df_gene.drop(rare_genes, axis=0)

    core_genes = df_gene.mean(axis=1) > max_genes
    core_genes = core_genes[core_genes == True].index

    # remove core genes from df_gene
    df_gene = df_gene.drop(core_genes, axis=0)
    info(df_gene.shape)
    ### Save
    info("Save dereplicated data tables")
    df_gene.to_csv(gene_outpath, sep='\t', index=True, header=True)
