#!/usr/bin/env python

"""
Jens Luebeck
UC San Diego, Bioinformatics & Systems Biology
jluebeck@ucsd.edu
"""

# Generate a tabular file containing BNX file molecules lengths in given size windows.

import argparse
import bisect
import os

import numpy as np


def get_mol_lens(bnxf):
    sizes = []
    with open(bnxf) as infile:
        for line in infile:
            if line.startswith("0"):
                sizes.append(float(line.rsplit()[2]))

    print("N molecules: " + str(len(sizes)))
    print("Median molecule length: " + str(np.median(sizes)))
    print("Mean molecule length: " + str(np.mean(sizes)))
    return sorted(sizes)


def bin_mol_sizes(ssizes, bs):
    maxsize = max(ssizes)
    print("Max molecule size: " + str(maxsize))
    bin_left_ends = np.arange(0, maxsize, bs)
    bincounts = [0]*len(bin_left_ends)
    for x in ssizes:
        ind = bisect.bisect(bin_left_ends, x) - 1
        bincounts[ind] += 1

    return zip(bin_left_ends, bincounts)


def write_bins(binsize_counts, bs, prefix):
    with open(prefix + "_binsize_counts.tsv", 'w') as outfile:
        outfile.write("bin_start\tnext_bin_start\tcount")
        for le, c in binsize_counts:
            oline = "\t".join([str(x) for x in [le, le+bs, c]]) + "\n"
            outfile.write(oline)

# Convert cmap file between versions. Currently supports 0.1 <--> 0.2 conversions.
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate a tabular file containing molecule counts in size windows "
                                                 "from BNX file")
    parser.add_argument("--bnx", "--input", help="Input CMAP file (file to be converted).", required=True)
    parser.add_argument("-s", help="Bin size (in bp) for reporting BNX molecule lengths (default 5000)", type=int,
                        default=5000)

    args = parser.parse_args()
    print("Reading BNX file")
    ssizes = get_mol_lens(args.bnx)
    print("Binning counts")
    binsize_counts = bin_mol_sizes(ssizes, args.s)
    prefix = os.path.splitext(os.path.basename(args.bnx))[0]
    write_bins(ssizes, args.s, prefix)
    print("Finished")
