#!/usr/bin/env python

"""
Jens Luebeck
UC San Diego, Bioinformatics & Systems Biology
jluebeck@ucsd.edu
"""

#Extract subset of xmap based on ID of ref, query or both.

import sys
import argparse


def extract_alignment_subset(xmapF, outfile_name, r_id_set, q_id_set, in_both):
    with open(xmapF) as infile, open(outfile_name, 'w') as outfile:
        for line in infile:
            if line.startswith("#"):
                outfile.write(line)
                if line.startswith("#h"):
                    head = line.rstrip().rsplit("\t")
                    head[0] = head[0].rsplit("#h ")[1]

            if not line.startswith("#"):
                fields = line.rsplit("\t")
                fDict = dict(zip(head, fields))
                queryID = fDict['QryContigID']
                refID = fDict['RefContigID']
                if in_both:
                    if refID in r_id_set and queryID in q_id_set:
                        outfile.write(line)

                else:
                    if refID in r_id_set or queryID in q_id_set:
                        outfile.write(line)


if __name__ == '__main__':
    #Parses the command line arguments
    parser = argparse.ArgumentParser(description="Subset of alignments from XMAP file")
    parser.add_argument("-i", "--input", help="name of input file",required=True)
    parser.add_argument("-o", "--outname", help="output XMAP filename",required=True)
    parser.add_argument("--r_id_list", help="Ref ID numbers", nargs='+',default=[])
    parser.add_argument("--q_id_list", help="Query ID numbers", nargs='+',default=[])
    parser.add_argument("--in_both", help="Elements in XMAP must appear in both --r_id_list AND --q_id_list to be "
                                          "extracted. By default xmap entry must be OR, not AND",action='store_true')
    args = parser.parse_args()

    if not args.q_id_list and not args.r_id_list:
        sys.stderr.write("--r_id_list or --q_id_list required\n")
        sys.exit(1)

    outfile_name = args.outname if args.outname.endswith(".xmap") else args.outname + ".xmap"
    q_id_set = set(args.q_id_list)
    r_id_set = set(args.r_id_list)

    extract_alignment_subset(args.input, outfile_name, r_id_set, q_id_set, args.in_both)
