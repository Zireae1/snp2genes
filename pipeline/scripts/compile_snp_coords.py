#!/usr/bin/env python3

import pandas as pd
import numpy as np
import pathlib
import scipy as sp
import scipy.cluster.hierarchy as sch
import sys
import re
import os.path
from lib.pandas_util import read_table_lz4
from lib.logging_util import info

if __name__ == "__main__":
    #### variable files
    path = sys.argv[1]
    outpath = sys.argv[2]
    snp_input = sys.argv[3]
    gene_input = sys.argv[4] 
    genes_info_path = sys.argv[5]
    ### replace:
    #path = '/Users/vdubinkina/Desktop/Lab/Rotation/cbarrios/PGeSM/nb/data/100122/'
    #outpath = path + 'snp_loc/'
    #snp_input = path + 'clean_data/sp-100122.denovo.derep.h.snp80.tsv'
    #gene_input = path + 'clean_data/sp-100122.denovo.derep.h.gene80.tsv'
    #snp_inpath = path + '/genome.norm.nucmer-spref.snps_table.filt-poly05-cvrg90.geno.tsv.lz4'
    rep_genomes_inpath = '/pollard/scratch/vdubinki/projects/snp2genes_clean/meta/species_to_reference_genome.tsv'
    #genes_info_path = path + '/genes_info.tsv'
    
    info("Data load")

    ### read snp file with names (i.e. global snp ids)
    #df_snp = read_table_lz4(snp_inpath)
    df_snp = pd.read_csv(snp_input, sep='\t',index_col=False, header=0)
    ### retrive rep_genome info to skip it
    rep_genome = pd.read_csv(rep_genomes_inpath, sep='\t', index_col=0, header=0)
    #rep_genome['representative'] = [re.sub("GUT_GENOME", "", x) for x in rep_genome['representative']]
    rep_genome=df_snp.columns.intersection(rep_genome['genome_id'])
    print(rep_genome)
    #rep_genome_id=re.sub("GUT_GENOME", "",rep_genome[0])
    # Create clean snp data
    snp_data = df_snp.iloc[:, 2:]

    ### leave only snp ids
    df_snp = df_snp['contig']+'_'+df_snp['pos'].astype(str)
    
    # Importing the gene data
    gene_data = pd.read_csv(gene_input, sep='\t', index_col=0, header=0)
    ### need to transpose it such that genomes are columns:
    #gene_data = gene_data.T
   
    # binarize gene data!!!
    #gene_data[gene_data>0]=1
    
    # make sure that they are ordered in the same way:
    snp_data = snp_data.sort_index(axis = 1)
    gene_data = gene_data.sort_index(axis = 1)
    
    ### retrieve centroids info
    genes_info = pd.read_csv(genes_info_path, sep='\t', index_col=False, header=0)
    ### leave only the map to centroid
    genes_info = genes_info[['gene_id', 'centroid_80']]
    ### leave only genes in selected accessory gene set
    genes_info = genes_info[genes_info['centroid_80'].isin(gene_data.index)]
    print(genes_info.shape)

    info("SNP compilation")
    for g in range(gene_data.shape[0]):
        print(g)
        info(".")
        test_gene = gene_data.index[g]
        ### create list of genomes where this gene is present
        test_gene_info = genes_info[genes_info['centroid_80']==test_gene]
        test_gene_info = test_gene_info.copy()
        ### now need to check if these genomes are in the derep subset...
        ### fix warning
        test_gene_info.loc[:,'genome_id'] = test_gene_info['gene_id'].str.split("_", n=2).str[:2].str.join("_")
        
        test_gene_info = test_gene_info[test_gene_info['genome_id'].isin(gene_data.columns)]
        #print(test_gene_info.shape)
    
        ### for this smaller set go and retrieve gene coords from the .gene files
        ### and snp coords from .snp files from mummer
        full_snp_coords = pd.DataFrame(columns=['P1', 'P2', 'C1', 'C2', 'gene_id2', 'MinDist2'])

        for i in range(test_gene_info.shape[0]): ### for every genome:
            #info("Genome")
            if(test_gene_info['genome_id'].values[i]==rep_genome):
                print("rep genome")
            else: 
                #info("else")
                ### load gene coordinates
                gene_coords_path = path + '/gene_annot/'+test_gene_info['genome_id'].values[i]+'.genes'
                gene_coords = pd.read_csv(gene_coords_path, sep='\t', index_col=False, header=0)

                ### make contig ID conversion
                ### NB!!!!! most dubious part
                gene_coords['contig_id'] = gene_coords['contig_id'].str.split("_").str[2].astype(int)-1
                ### find a gene:
                gene_coords = gene_coords[gene_coords['gene_id']==test_gene_info['gene_id'].values[i]]
                #print(gene_coords)

                ### grab SNP coords in the same genome
                snp_coords_path = path + '/genome/'+test_gene_info['genome_id'].values[i]+'/genome.norm.nucmer-spref.snps'
            
                if(os.path.isfile(snp_coords_path)):
                    snp_coords = pd.read_csv(snp_coords_path, skiprows=4, sep='\t', index_col=False, header=None)
                    snp_coords.columns = ['P1', 'Sub1', 'Sub2', 'P2', 'Buff', 'Dist', 'LenR', 'LenQ', 'FRM', 'Tags','C1', 'C2']
                    snp_coords = snp_coords[['P1', 'P2', 'C1', 'C2']]
            
                    #### check if GUT_GENOME or UHGG -> 2/1 
                    snp_coords['C2'] = [re.sub("C", "", x) for x in snp_coords['C2'].str.split("_").str[2]]                               
                    snp_coords['C2'] = snp_coords['C2'].astype(int)
                    ### grab contig for the first matched gene (but there might be muti-copy!!!)
                    snp_coords = snp_coords[snp_coords['C2']==gene_coords['contig_id'].values[0]]

                    ### add gene info to this table:
                    snp_coords['gene_id2'] = gene_coords['gene_id'].values[0]
                    ### calc distances:
                    snp_coords['MinDist2'] = [min(abs(x-gene_coords['start'].values[0]), abs(x-gene_coords['end'].values[0])) for x in snp_coords['P2']]

                    ### save into one giant table:
                    full_snp_coords = pd.concat([full_snp_coords,snp_coords], ignore_index=True)
                else:
                    print("No snp file!!!")

        print(full_snp_coords.shape)
        ### now fix it a bit:
        ### create a unique SNP id (using global reference)
        full_snp_coords['SNP_id'] = full_snp_coords['C1']+'_'+full_snp_coords['P1'].astype(str)
        print(full_snp_coords['SNP_id'].value_counts().shape)
        ### select only snps in the list
        full_snp_coords = full_snp_coords[full_snp_coords['SNP_id'].isin(df_snp)]
        print(full_snp_coords['SNP_id'].value_counts().shape)
        ### save this info:
        full_snp_coords_outpath = outpath +'/'+ test_gene +'.snp.coords.tsv'
        full_snp_coords.to_csv(full_snp_coords_outpath, sep='\t', index=True, header=True)


