#!/usr/bin/env python3
"""
Take a gene_info.txt.lz4 file and convert it into a genome-by-centroid matrix.

Assumes that input file is lz4 compressed.
Assumes that input file has a header, and is indexed by a column named gene_id.
Writes output as an xarray DataArray to NetCDF.
"""

import pandas as pd
import subprocess
import sys

if __name__ == "__main__":
    gene_path = sys.argv[1]
    genome_list = sys.argv[2]
    column = sys.argv[3]
    out_path = sys.argv[4]
    
    # grab full list of genomes
    genomes = pd.read_csv(genome_list, sep='\t', index_col=0, header=0)

    with subprocess.Popen(["lz4cat", gene_path], stdout=subprocess.PIPE) as f:
        reference_gene = pd.read_table(f.stdout, usecols=["gene_id", column])
    reference_gene = (
        reference_gene.assign(genome_id=lambda x: x.gene_id.str.split("_", n=2).str[:2].str.join("_"))[
            ["genome_id", column]
        ]
        .value_counts()
        .rename_axis(["genome_id", "gene_id"])
        .to_xarray()
        .fillna(0)
        .astype(int)
    )
    
    ### here need to filter down
    # Filter down reference_gene to include only genomes in genome_list
    valid_genomes = genomes.index  # Get the list of valid genome IDs from the genome_list file
    # Ensure only genome_ids present in reference_gene are used
    valid_genomes = valid_genomes.intersection(reference_gene["genome_id"].values)
    print(valid_genomes.shape)
    reference_gene = reference_gene.sel(genome_id=valid_genomes)
    print(reference_gene.shape)
    reference_gene.to_dataset(name="copy_number").to_netcdf(
        out_path, encoding=dict(copy_number=dict(zlib=True, complevel=5))
    )
