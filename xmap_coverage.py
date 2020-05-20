#!/usr/bin/env python

"""
Jens Luebeck
UC San Diego, Bioinformatics & Systems Biology
jluebeck@ucsd.edu
"""

import argparse
from collections import defaultdict

#Extract coverage in xmap file, total and per-reference ID
#TODO: Add capability to extract coverage for a region in a BED file.


#parses map file and extracts coverage (total, and per reference ID)
def getCov(xmapF, contigID = None):
	head = []
	refLens = {}
	qryLens = {}
	refQTotLen = defaultdict(float)

	with open(xmapF) as infile:
		for line in infile:
			if line.startswith("#h"):
				head = line.rstrip().rsplit("\t")
				head[0] = head[0].rsplit("#h ")[1]

			elif line and not line.startswith("#"):
				fields = line.rstrip().rsplit("\t")
				fDict = dict(zip(head, fields))
				refID = fDict['RefContigID']
				qryID = fDict['QryContigID']
				if contigID and refID != contigID:
					continue

				if refID not in refLens:
					refLens[refID] = float(fDict['RefLen'])

				if qryID not in qryLens:
					qryLens[qryID] = float(fDict['QryLen'])



				alnLen = float(fDict['RefEndPos']) - float(fDict['RefStartPos'])
				refQTotLen[refID]+=alnLen


	refCovs = defaultdict(float)
	totAlnLen = 0.
	totRefLen = 0.
	totQryLen = 0.
	for refID, refLen in refLens.items():
		refCovs[refID] = refQTotLen[refID]/refLen
		totAlnLen+=refQTotLen[refID]
		totRefLen+=refLen

	for qryID, qryLen in qryLens.items():
		totQryLen+=qryLen

	alnCov = totAlnLen/totRefLen
	LWCov = totQryLen/totRefLen

	return refCovs, alnCov, LWCov 


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="Extract aligned bionano molecules from contigs "
												 "at given reference genome positions")
	parser.add_argument("-x","--xmap",help="xmap file",required=True)
	parser.add_argument("-c","--contig",help="contig ID to get coverage from", default="")
	# parser.add_argument("--bed", help="BED file of reference regions toget coverage from")
	args = parser.parse_args()

	refCovs, alnCov, LWCov  = getCov(args.xmap, args.contig)
	sortedRefIDs = sorted(refCovs.keys(), key=lambda x: int(x))
	print("RefID\tCoverage of alignments")
	for refID in sortedRefIDs:
		print(refID + "\t" + str(refCovs[refID]))

	print("Total coverage of alignments: " + str(alnCov))
	print("Estimated reference coverage (Lander-Watterman): " + str(LWCov))

