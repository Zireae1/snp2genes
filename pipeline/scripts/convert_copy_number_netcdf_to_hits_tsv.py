#!/usr/bin/env python3

import xarray as xr
import sys

if __name__ == "__main__":
    inpath = sys.argv[1]
    outpath = sys.argv[2]

    data = xr.load_dataarray(inpath)
    data = data.to_pandas().astype(int)
    data.to_csv(outpath, sep="\t")
