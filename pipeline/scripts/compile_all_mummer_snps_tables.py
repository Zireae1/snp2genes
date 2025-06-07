#!/usr/bin/env python3

import pandas as pd
import sys
from tqdm import tqdm
import numpy as np
from lib.logging_util import info


if __name__ == "__main__":
    # For computational reasons, pre-filter positions to those with at least
    # this fraction of genomes showing a base mismatch [ACGT].
    min_mismatch_rate = float(sys.argv[1])  # 0.0
    arg_list = sys.argv[2:]

    info("START: Collecting input data.")
    genome_list = []
    mismatch_table = {}
    ref_base = []
    all_alignment = {}
    for arg in tqdm(arg_list):
        genome, mismatch_path, aligns_path = arg.split(":")
        genome_list.append(genome)

        # This parses the necessary columns from the MUMMER `show-snps -ClrT`
        # format.
        mismatch = pd.read_table(
            mismatch_path,
            skiprows=4,
            names=[
                "pos",
                "ref_base",
                "base",
                "_3",
                "_4",
                "_5",
                "_6",
                "_7",
                "_8",
                "_9",
                "contig",
                "_11",
            ],
        )[
            # NOTE: A reference base of "." means that there's an insertion in
            # the query genome.
            lambda x: (x.ref_base != ".")
        ]
        mismatch_table[genome] = (
            mismatch[["contig", "pos", "base"]].set_index(["contig", "pos"]).base
        )
        ref_base.append(mismatch[["contig", "pos", "ref_base"]])

        # This parses the necessary columns from the MUMMER `show-coords -rgT`
        # format.
        align = pd.read_table(
            aligns_path,
            skiprows=4,
            names=[
                "start_pos",
                "stop_pos",
                "_2",
                "_3",
                "_4",
                "_5",
                "identity",
                "contig",
                "_8",
            ],
        )
        all_alignment[genome] = align[["contig", "start_pos", "stop_pos", "identity"]]
    info("END: Collecting input data.")

    info("START: Collating reference positions.")
    ref_base = pd.concat(ref_base).drop_duplicates().set_index(["contig", "pos"])
    info("END: Collating reference positions.")

    info("START: Collating mismatches across all genomes.")
    # Populate a pre-allocated mismatch matrix by reindexing the mismatch
    # table inputs by the list of all "interesting" positions.
    position_list0 = ref_base.index
    all_mismatch = np.empty((len(position_list0), len(genome_list)), dtype="S1")
    for i, genome in tqdm(list(enumerate(genome_list))):
        all_mismatch[:, i] = mismatch_table[genome].reindex(position_list0).fillna("?")
    all_mismatch = pd.DataFrame(all_mismatch, index=position_list0, columns=genome_list)
    del mismatch_table  # No longer needed and potentially large.
    info("END: Collating mismatches across all genomes.")

    info("START: Collating alignment coverage matrix across all genomes.")
    # Pre-filter a list of positions by their rate of alternative bases.
    position_list1 = (
        all_mismatch.isin([b"A", b"C", b"G", b"T"])
        .mean(1)
        .sort_index()[lambda x: x > min_mismatch_rate]
        .index
    )
    # Iterate through each genome
    all_cov = {}
    for genome in tqdm(genome_list):
        # Construct a bit vector of positions.
        y = pd.Series(np.zeros(len(position_list1), dtype="bool"), index=position_list1)
        # Iterate through each valid alignment
        for _, x in all_alignment[genome].iterrows():
            # Flip the bits in the genome vector corresponding to all mismatch
            # positions contained within that alignment.
            y.loc[(x.contig, x.start_pos) : (x.contig, x.stop_pos)] = True
        all_cov[genome] = y
    all_cov = pd.DataFrame(all_cov, index=position_list1)
    # Unflip bits for any position with an indel (aligned but not covered).
    all_cov[(all_mismatch.loc[position_list1, genome_list] == b".")] = False
    info("END: Collating alignment coverage matrix across all genomes.")

    info("START: Merging coverage and mismatch data.")
    genotype = all_mismatch.loc[position_list1]
    genotype[all_cov.loc[position_list1] & (genotype == b"?")] = b"="
    info("END: Merging coverage and mismatch data.")

    info("START: Writing output.")
    print("contig", "pos", *genotype.columns, sep="\t")
    for (contig, pos), x in genotype.apply(lambda x: x.str.decode("utf-8")).iterrows():
        print(contig, pos, *x, sep="\t")
    info("END: Writing output.")

    # info("START: Filtering to bi-allelic positions.")
    # allele_ratios = genotype.apply(
    #     lambda x: x.value_counts(normalize=True), axis=1
    # ).fillna(0)
    # position_list2 = allele_ratios[
    #     lambda x: (x[[b"=", b"A", b"C", b"G", b"T"]] > 0).sum(1) == 2
    # ].index
    # info("END: Filtering to bi-allelic positions.")
    #
    # info("START: Binarizing.")
    # binary_genotype = (genotype.loc[position_list2] != b"=").astype(float)
    # binary_genotype[genotype.loc[position_list2] == b"?"] = np.nan
    # info("END: Binarizing.")
    #
    # info("START: Writing output.")
    # binary_genotype.stack().to_xarray().to_netcdf(outpath)
    # info("END: Writing output.")
