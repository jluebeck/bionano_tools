#!/usr/bin/env python

"""
Jens Luebeck
UC San Diego, Bioinformatics & Systems Biology
jluebeck@ucsd.edu
"""

#Extract highest scoring alignment from .xmap file and write to file
#multiple returned if tied

import sys
import os
import argparse

def find_best_aln(xmapF,xmapF_out):
	bestScore = 0
	scoreToXmapID = {}
	with open(xmapF) as infile, open(xmapF_out,'w') as outfile:
		for line in infile:
			if line.startswith("#h"):
				head = line.rstrip().rsplit("\t")

			if not line.startswith("#"):
				fields = line.rsplit("\t")
				fDict = dict(zip(head,fields))
				currScore = float(fDict['Confidence'])
				if currScore > bestScore:
					bestScore = currScore

				if currScore not in scoreToXmapID:
					scoreToXmapID[currScore] = set()

				scoreToXmapID[currScore].add(line)

			else:
				outfile.write(line)

		if scoreToXmapID:
			for line in scoreToXmapID[bestScore]:
				outfile.write(line)


if __name__ == '__main__':
	#Parses the command line arguments
	parser = argparse.ArgumentParser(description="Extract best alignment from XMAP file")
	parser.add_argument("-i", "--input", help="name of input file",required=True)
	parser.add_argument("-o", "--outname", help="extracted file name (default to [inputname]_BESTALN.xmap")
	args = parser.parse_args()

	outname = args.outname
	if not outname:
		outname = args.input[:args.input.index(".xmap")] + "_BESTALN" + args.input[args.input.index(".xmap"):]

	find_best_aln(args.input,outname)
 