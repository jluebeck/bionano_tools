#!/usr/bin/env python

"""
Jens Luebeck
UC San Diego, Bioinformatics & Systems Biology
jluebeck@ucsd.edu
"""

import os
import sys
import argparse

allowed_versions = ["0.1","0.2"]
map_1_head_names = ["#h CMapId","ContigLength","NumSites","SiteID","LabelChannel","Position","StdDev","Coverage","Occurrence","ChimQuality"]
map_1_head_types = ["#f int","float","int","int","int","float","float","float","float","float"]
map_2_head_names = map_1_head_names + ["SegDupL","SegDupR","FragileL","FragileR","OutlierFrac","ChimNorm","Mask"]
map_2_head_types =  map_1_head_types + ["float","float","float","float","float","float","Hex"]

def parse_cmap(cmap_file):
	header_lines = []
	head_names = []
	line_dict_vector = []
	field_names = []
	with open(cmap_file) as infile:
		for line in infile:
			if line.startswith("#"):
				header_lines.append(line)
				if line.startswith("#h"):
					head_names = line.rstrip().rsplit("\t")

			else:
				line_dict = dict(zip(head_names,line.rstrip().rsplit("\t")))
				line_dict_vector.append(line_dict)

	return header_lines,line_dict_vector

def write_cmap(header_lines,line_dict_vector,outname,map_head_names,map_head_types,version):
	with open(outname,'w') as outfile:
		#write header
		for line in header_lines:
			if line.startswith("# CMAP File Version:"):
				line_fields = line.rstrip().rsplit(":")
				line = line_fields[0] + ":\t" + version + "\n"

			elif line.startswith("#h"):
				break

			outfile.write(line)

		outfile.write("\t".join(map_head_names) + "\n")
		outfile.write("\t".join(map_head_types) + "\n")

		#write data
		for line_dict in line_dict_vector:
			outline = line_dict[map_head_names[0]]
			for ind, x in enumerate(map_head_names[1:]):
				outline+="\t"
				if x in line_dict:
					outline+=line_dict[x]

				#couldn't find the data field name, using a dummy value
				else:
					#find out what data type it is supposed to be
					dtype = map_head_types[ind+1]
					if dtype == "float":
						outline+="0.0"

					else:
						outline+="0"

			outline+="\n"
			outfile.write(outline)

def check_cmap_version(cmap_file):
	with open(cmap_file) as infile:
		for line in infile:
			if line.startswith("#"):
				if "CMAP File Version:" in line:
					version = line.rstrip().rsplit(":\t")[1]
					return version

	return "CMAP version not found."


#Convert cmap file between versions. Currently supports 0.1 <--> 0.2 conversions.
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert BioNano CMAP file between versions. Currently supports 0.1 <--> 0.2 CMAP conversions.")
    parser.add_argument("-i", "--input", help="Input CMAP file (file to be converted).", required=True)
    parser.add_argument("-o", "--output", help="Filename for output converted CMAP.")
    parser.add_argument("-v", help="CMAP version to convert input file to.",choices=allowed_versions,required=True)

    args = parser.parse_args()

    conv_version = args.v
    if not args.output:
    	args.output = os.path.splitext(args.input)[0] + "_v0_" + conv_version[-1] + ".cmap"

    curr_version = check_cmap_version(args.input)
    if not curr_version in allowed_versions:
    	sys.stderr.write("Error: Invalid CMAP version: " + curr_version + "\n")

    elif curr_version == conv_version:
    	sys.stdout.write("Input CMAP version (" + curr_version + ") already at desired version.\n")

    else:
    	sys.stdout.write("Converting input CMAP to version " + conv_version)
    	sys.stdout.write("\nConverted CMAP filename: " + args.output + "\n")
    	#do stuff
    	header_lines,line_dict_vector = parse_cmap(args.input)

    	if conv_version == "0.1":
    		#write 0.1
    		write_cmap(header_lines,line_dict_vector,args.output,map_1_head_names,map_1_head_types,conv_version)

    	else:
    		#write 0.2
    		write_cmap(header_lines,line_dict_vector,args.output,map_2_head_names,map_2_head_types,conv_version)


    sys.stdout.write("CMAP conversion finished.\n")