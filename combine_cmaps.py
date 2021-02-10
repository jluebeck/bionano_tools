#!/usr/bin/env python

"""
Jens Luebeck
UC San Diego, Bioinformatics & Systems Biology
jluebeck@ucsd.edu
"""

# Combine CMAPs, creating new CMAP ids. Converted IDs will be printed for second file.
# GIVEN CMAPS MUST BE IN THE SAME FORMAT! (versioning must be same)

import os
import sys
import argparse


def parse_cmap(cmap_file):
    header_lines = []
    head_names = []
    line_dict_vector = []
    with open(cmap_file) as infile:
        for line in infile:
            if line.startswith("#"):
                header_lines.append(line)
                if line.startswith("#h"):
                    head_names = line.rstrip().lstrip("#h ").rsplit("\t")

            else:
                line_dict = dict(zip(head_names,line.rstrip().rsplit("\t")))
                line_dict_vector.append(line_dict)

    return header_lines, line_dict_vector, head_names


def write_cmap(header_lines, line_dict_vector1, line_dict_vector2, outname, head_names):
    with open(outname, 'w') as outfile:
        for line in header_lines:
            outfile.write(line)

        for line in line_dict_vector1:
            outfile.write("\t".join([line[k] for k in head_names]) + "\n")

        for line in line_dict_vector2:
            outfile.write("\t".join([line[k] for k in head_names]) + "\n")



def get_max_id(line_dict_vector):
    hitset = set()
    for x in line_dict_vector:
        # print(x.keys())
        hitset.add(int(x["CMapId"]))

    return max(hitset)


def remap_ids(line_dict_vector, max_id):
    remap = {}
    seen_ids = set()
    # print(line_dict_vector[0]["CMapId"])
    for x in line_dict_vector:
        cid = int(x["CMapId"])
        if cid in seen_ids:
            new_cid = remap[cid]

        else:
            max_id += 1
            new_cid = max_id
            remap[cid] = new_cid
            seen_ids.add(cid)

        x["CMapId"] = str(new_cid)

    # print(line_dict_vector[0]["CMapId"])
    return remap, line_dict_vector


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Combine two cmaps (merge in IDs). Must be same CMAP file version!")
    parser.add_argument("--cmap1", help="Input CMAP file 1", required=True)
    parser.add_argument("--cmap2", help="Input CMAP file 2", required=True)
    parser.add_argument("-o", help="output filename prefix", required=True)

    args = parser.parse_args()
    outfile = args.o
    #read cmap 1, count number of ids
    print("Reading CMAPs")
    header_lines1, line_dict_vector_1, head_names1 = parse_cmap(args.cmap1)
    max_id = get_max_id(line_dict_vector_1)
    header_lines2, line_dict_vector_2, head_names2 = parse_cmap(args.cmap2)
    remap, rldv_2 = remap_ids(line_dict_vector_2, max_id)
    print("Writing combined CMAP")
    write_cmap(header_lines1, line_dict_vector_1, line_dict_vector_2, outfile, head_names1)
    print("cmap2_original_id\tcmap2_merged_id")
    for old, new in sorted(remap.items()):
        print(str(old) + "\t" + str(new))

    print(outfile)



