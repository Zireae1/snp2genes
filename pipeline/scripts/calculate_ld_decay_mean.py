#!/usr/bin/env python3
# "{input.script} {wildcards.nsnps} {wildcards.seed} {input.geno} {output}"

import pandas as pd
from scipy.spatial.distance import pdist
import sys
#from lib.pandas_util import read_table_lz4
from lib.pandas_util import read_table_tsv
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt


#DISTANCE_BINS = [2**i for i in range(24)]
#LD_BINS = np.linspace(0, 1, num=21)

xlim = np.array([0.5, 1e5])
DISTANCE_BINS = np.unique(np.floor(np.logspace(*np.log10(xlim), num=51)).astype(int))
LD_BINS = np.linspace(0, 1, num=51)


if __name__ == "__main__":
    geno_inpath = sys.argv[1]
    stats_outpath = sys.argv[2]
    #fig_outpath = sys.argv[3]

    # Load genotype data
    #geno = read_table_lz4(geno_inpath, index_col=["contig", "pos"])
    geno = read_table_tsv(geno_inpath, index_col=["contig", "pos"])

    # Map contigs to integers
    contig_list = list(geno.index.to_frame().contig.unique())
    contig_to_idx = pd.Series(range(len(contig_list)), index=contig_list)

    # Calculate compacted, pairwise arrays for all SNP pairs.
    contig_diff = pdist(
        geno.index.to_frame()[["contig"]].replace(contig_to_idx), metric="cityblock"
    )
    pos_diss = pdist(geno.index.to_frame()[["pos"]], metric="cityblock")
    corr_diss = pdist(geno.fillna(0.5), metric="correlation")

    # Compile a results table.
    # NOTE: This table is anonymized. It'll take additional work if we want to
    # keep track of _which_ SNPs are correlated.
    ld_data = pd.DataFrame(
        dict(ld=(1 - corr_diss) ** 2, pos_diss=pos_diss, same_contig=contig_diff == 0)
    )[lambda x: x.same_contig]#[lambda x: x.same_contig]


    bin_stats = ld_data.groupby(pd.cut(ld_data.pos_diss, bins=DISTANCE_BINS, right=False)).ld.apply(lambda x: pd.Series(dict(ldmean=x.mean(), count=len(x)))).unstack().assign(bin_left=DISTANCE_BINS[:-1])
    bin_stats.to_csv(stats_outpath, sep='\t')

    # LD-decay figure
    #fig, ax = plt.subplots()
    #ax.hist2d('pos_diss', 'ld', data=ld_data, bins=(DISTANCE_BINS, LD_BINS), norm=mpl.colors.LogNorm(), cmap='Greys')
    #ax.plot('bin_left', 'ld90', data=bin_stats)
    #ax.set_xscale('log')
    #fig.savefig(fig_outpath)
