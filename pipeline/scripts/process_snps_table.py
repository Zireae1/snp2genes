#!/usr/bin/env python3

"""
{input.script} {params.poly} {params.cvrg} {params.align} {input.snps} {output}
"""

import sys
from lib.logging_util import info
from lib.pandas_util import read_table_lz4


if __name__ == "__main__":
    # Positions must have minimum fraction of genomes with a different base.
    poly_thresh = float(sys.argv[1])
    # Positions must have a minimum fraction of genomes with NOT "?"
    cvrg_thresh = float(sys.argv[2])
    inpath = sys.argv[3]

    data = read_table_lz4(inpath, index_col=["contig", "pos"])
    info(f"Input data has {data.shape[0]} positions across {data.shape[1]} genomes.")
    poly = data.isin(["A", "C", "G", "T"]).mean(1)
    cvrg = data.isin(["A", "C", "G", "T", "="]).mean(1)
    allele_ratios = data.apply(
        lambda x: x.value_counts(normalize=True), axis=1
    ).fillna(0)
    num_alleles = (allele_ratios[["=", "A", "C", "G", "T"]] > 0).sum(1)

    data = data[(poly >= poly_thresh) & (cvrg > cvrg_thresh) & (num_alleles == 2)]
    info(f"After filtering for data has {data.shape[0]} positions across {data.shape[1]} genomes.")

    binary_genotype = (data != "=").astype(float)
    binary_genotype[data == "?"] = float("nan")
    binary_genotype.to_csv(sys.stdout, sep="\t")
