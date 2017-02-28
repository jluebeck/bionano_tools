#!/usr/bin/env python

"""
Jens Luebeck
UC San Diego, Bioinformatics & Systems Biology
jluebeck@ucsd.edu
"""

#Make the seg marks figure 

import sys
import os
import argparse
from matplotlib import pyplot as plt
import matplotlib.patches as patches
from matplotlib import colors as mcolors
from intervaltree import Interval, IntervalTree


#parse cmap
def parse_cmap(cmapf):
    cmaps = {}
    contigCovs = {}
    with open(cmapf) as infile:
        for line in infile:
            if line.startswith("#h"):
                head = line.rstrip().rsplit()[1:]

            elif not line.startswith("#"):
                fields = line.rstrip().rsplit()
                fD = dict(zip(head,fields))
                if fD["CMapId"] not in cmaps:
                    cmaps[fD["CMapId"]] = [None]*int(fD["NumSites"])
                    contigCovs[fD["CMapId"]] = [None]*int(fD["NumSites"])

                if fD["LabelChannel"] == "1":
                    cmaps[fD["CMapId"]][int(fD["SiteID"])-1] = float(fD["Position"])
                    contigCovs[fD["CMapId"]][int(fD["SiteID"]) - 1] = float(fD["Coverage"])
    #print cmaps
    return cmaps,contigCovs

#parse xmap
def parse_xmap(xmapf):
    detailFields = ["QryContigID","RefContigID","Orientation","Confidence","QryLen","RefLen","QryStartPos","QryEndPos","RefStartPos","RefEndPos"]
    xmapAln = {}
    xmapPair = {}
    with open(xmapf) as infile:
        for line in infile:
            if line.startswith("#h"):
                head = line.rstrip().rsplit()[1:]

            elif not line.startswith("#"):
                fields = line.rstrip().rsplit()
                fD = dict(zip(head,fields))
                alnstring = ")" + fD["Alignment"] + "("
                xmapAln[fD["XmapEntryID"]] = alnstring.rsplit(")(")[1:-1]

                xmapPair[fD["XmapEntryID"]] = {x:fD[x] for x in detailFields}

    return xmapAln,xmapPair

def parse_cycles(fname="/home/jens/Dropbox/BafnaLab/bionano_analysis/kt48/KT48_amplicon1_structures.txt"):
    with open(fname) as infile:
        cycleLens = {}
        cycleD = {}
        for line in infile:
            if line.startswith("Segment"):
                fields = line.rstrip().split()
                start = int(fields[3])
                end = int(fields[4])
                cycleLens[fields[1]] = [end-start,start,end,fields[2]]
            elif line.startswith("Cycle"):
                fields = line.rstrip().rsplit(";")
                segs = fields[2].rsplit("=")[1].rsplit(",")
                segN = [x[:-1] for x in segs]
                dirs = [x[-1] for x in segs]
                cycleD[fields[0].rsplit("=")[1]] = (segN,dirs)

    return cycleD,cycleLens

def get_cycle_colors():
    colors = dict(mcolors.BASE_COLORS, **mcolors.CSS4_COLORS)
    by_hsv = sorted((tuple(mcolors.rgb_to_hsv(mcolors.to_rgba(color)[:3])), name)
                    for name, color in colors.items())
    sorted_names = [name for hsv, name in by_hsv]
    sorted_names = sorted_names[18:] + sorted_names[:18]
    return sorted_names

def plot_cyclebar(cycleD,cycleLens,cycleNumber):
    height = 0.15
    segN, dirs = cycleD[cycleNumber]
    #print segN
    #print dirs

    sorted_names = get_cycle_colors()
    spacing = len(sorted_names)/len(cycleLens)
    #print sorted_names[(2 - 1) * spacing]
    #print sorted_names[(12 - 1) * spacing]
    #print spacing
    currStart = 0
    for ind,x in enumerate(segN):
        ax.add_patch(
            patches.Rectangle(
                (currStart, -0.4),  # (x,y)
                cycleLens[x][0],  # width
                height,  # height
                facecolor=sorted_names[(int(x)-1)*spacing],
            )
        )
        #print (int(x) - 1) * spacing
        currStart+=cycleLens[x][0]

def parseGenes(chrom):
    print("Building interval tree for chrom " + chrom)
    t = IntervalTree()
    with open("/home/jens/Dropbox/BafnaLab/ref_genomes/UCSC/refGene.txt") as infile:
        for line in infile:
            fields = line.rsplit("\t")
            #gene = fields[-2]
            currChrom = fields[2]
            tstart = int(fields[4])
            tend = int(fields[5])
            if chrom == currChrom:
                t[tstart:tend] = fields

    return t

def plot_genes(cycleD,cycleLens,cycleNumber):
    superSegL,superSegPons,superSegDir = super_segs(cycleD,cycleLens,cycleNumber)
    #print superSegL
    #print superSegPons
    #print superSegDir
    currStart = 0
    chrIntTree = {}
    for seg,pTup,strand in zip(superSegL,superSegPons,superSegDir):
        relGenes = {}
        chrom = pTup[0]
        if chrom not in chrIntTree:
            chrIntTree[chrom] = parseGenes(chrom)
        overlappingT = chrIntTree[chrom][pTup[1]:pTup[2]]
        for i in overlappingT:
            #print i.data
            gene = i.data[-4].rsplit("-")[0]
            tstart = int(i.data[4])
            tend = int(i.data[5])
            if not (gene.startswith("LOC") or gene.startswith("LINC")):
                if gene not in relGenes:
                    relGenes[gene] = (tstart,tend)
                else:
                    oldTStart = relGenes[gene][0]
                    oldTEnd = relGenes[gene][1]
                    if tstart < oldTStart:
                        oldTStart = tstart
                    if tend > oldTEnd:
                        oldTEnd = tend
                    relGenes[gene] = (oldTStart,oldTEnd)


        for ind,i in enumerate(relGenes):
            truncStart = False
            truncEnd = False
            tstart,tend = relGenes[i]
            normStart = currStart + max(0,tstart-pTup[1])
            if tstart < pTup[1]:
                truncStart = True
            normEnd = currStart + min(pTup[2]-pTup[1],tend-pTup[1])
            if tend > pTup[2]:
                truncEnd = True

            offset = 0.05*ind
        # ax.add_patch(
        #     patches.Rectangle(
        #         (currStart, -0.4),  # (x,y)
        #         cycleLens[x][0],  # width
        #         height,  # height
        #         facecolor=sorted_names[(int(x)-1)*spacing],
        #     )
        # )
            plt.text(normStart - 10000, -0.55-offset, i,size=4)
            ax.plot([normStart, normEnd], [-0.55-offset, -0.55-offset], color="k",linewidth=2)
            if truncEnd:
                ax.scatter([normEnd],[-0.55-offset], color="r",s=0.3)

            #print (int(x) - 1) * spacing
        currStart+=(pTup[2]-pTup[1])


def super_segs(cycleD,cycleLens,cycleNumber):
    # make super-segs
    segN, dirs = cycleD[cycleNumber]
    prevEnd = -2
    prevChr = None
    superSegL = []
    superSegPons = []
    superSegDir = []
    superSegs = 0
    for ind, i in enumerate(segN):
        currLens = cycleLens[i]
        currStart, currEnd, chrom = currLens[1:4]
        currStart = int(currStart)
        currEnd = int(currEnd)
        if currStart != prevEnd+1 or prevChr != chrom:
            superSegs += 1
            superSegL.append(superSegs)
            superSegPons.append((chrom, currStart, currEnd))
            superSegDir.append(dirs[ind])

        else:
            superSegPons[-1] = (superSegPons[-1][0], superSegPons[-1][1], currEnd)

        prevEnd = currEnd
        prevChr = chrom

    print(len(segN),len(superSegL))
    return superSegL,superSegPons,superSegDir

#plot function
def plot_alignment(ref_cmaps,ref_contigConv,query_cmaps,query_contigCovs,xmapAln,xmapPairAll,xmapID):
    refBaseX = 0
    refBaseY = 0
    height = 0.5
    spaceBars = 0.9
    xmapPair = xmapPairAll[xmapID]
    #print xmapPair


    REnd = float(xmapPair["RefEndPos"])
    RStart = float(xmapPair["RefStartPos"])
    revOrient = False
    QEnd = float(xmapPair["QryEndPos"])
    QStart = float(xmapPair["QryStartPos"])
    RLen = float(xmapPair["RefLen"])
    QLen  = float(xmapPair["QryLen"])
    offset = QStart
    if xmapPair["Orientation"] == "-":
        revOrient = True
        QStart = float(xmapPair["QryEndPos"])
        QEnd = float(xmapPair["QryStartPos"])
        offset = QLen - QEnd


    #print RLen,QLen

    #REFERENCE CONTIG
    ax.add_patch(
        patches.Rectangle(
            (refBaseX, refBaseY),  # (x,y)
            RLen,  # width
            height,  # height
            facecolor="cornflowerblue",
        )
    )

    print ref_cmaps[xmapPair["RefContigID"]][:10]
    for pos in ref_cmaps[xmapPair["RefContigID"]]:
        ax.plot([pos, pos], [refBaseY, refBaseY+height], color="k",linewidth=0.3)

    #QUERY CONTIG
    ax.add_patch(
        patches.Rectangle(
            (RStart-offset, refBaseY+spaceBars),  # (x,y)
            QLen,  # width
            height,  # height
            facecolor="mediumspringgreen",
        )
    )

    #print ref_cmaps[xmapPair["RefContigID"]][:10]
    for pos in query_cmaps[xmapPair["QryContigID"]]:
        if revOrient:
            pos = QLen - pos
        pos = pos + (RStart-offset)
        ax.plot([pos, pos], [refBaseY+spaceBars, refBaseY+spaceBars+height], color="k",linewidth=0.3)

    for p in xmapAln[xmapID]:
        r,q = p.split(",")
        x1 = ref_cmaps[xmapPair["RefContigID"]][int(r)-1]
        x2 = query_cmaps[xmapPair["QryContigID"]][int(q)-1]
        if revOrient:
            x2 = QLen - x2
        x2 = x2 + (RStart - offset)
        ax.plot([x1,x2],[refBaseY+height,refBaseY+spaceBars],color='gray',alpha=0.8,linewidth=0.25)

    ax.set_xlim([-100000 + min(RStart,RStart-offset), max(REnd,QEnd)+offset+100000])
    ax.set_ylim([-2, 2])
    print x1,x2

if __name__ == '__main__':
    # Parses the command line arguments
    parser = argparse.ArgumentParser(description="Plot an alignment of BioNano XMAP")
    parser.add_argument("-r", "--ref", help="reference cmap", required=True)
    parser.add_argument("-q", "--query", help="query cmap", required=True)
    parser.add_argument("-x", "--xmap", help="xmap file", required=True)
    parser.add_argument("-o", "--outname", help="output filename prefix")
    # parser.add_argument("--each_query", help="extract the best alignment for each query CMAP aligned, default false",
    #                     action='store_true')
    # parser.add_argument("--each_ref", help="extract the best alignment for each ref CMAP aligned, default false",
    #                     action='store_true')
    args = parser.parse_args()

    outname = args.outname

    cycleNumber = "1"
    xmapID = "2"

    ref_cmaps,ref_contigCovs = parse_cmap(args.ref)
    query_cmaps,query_contigCovs = parse_cmap(args.query)
    xmapAln, xmapPair = parse_xmap(args.xmap)

    # print len(ref_cmaps),len(ref_contigCovs)
    # print len(query_cmaps),len(query_contigCovs)
    # print len(xmapAln)
    print xmapPair

    plt.figure(figsize=(16,8.5))
    ax =  plt.gca()

    plot_alignment(ref_cmaps, ref_contigCovs, query_cmaps, query_contigCovs, xmapAln, xmapPair, xmapID)
    cycleD, cycleLens = parse_cycles()
    plot_cyclebar(cycleD,cycleLens,cycleNumber)

    plot_genes(cycleD,cycleLens,cycleNumber)

    plt.savefig(args.outname + '.png', dpi=300)
    plt.close()
# if not outname:
    #     outname = args.input[:args.input.index(".xmap")] + "_BESTALN" + args.input[args.input.index(".xmap"):]


