## Bionano tools
### Assorted scripts for working with Bionano data.

This collection of scripts does not require any additional libraries to be installed and supports both Python2 and Python3.


```
convert_cmap_version.py
Converts between CMAP versions. Outputs a new CMAP
Usage: python convert_cmap_version.py -i input.cmap -v [version to convert to (0.1/0.2)] [-o optional output CMAP filename]
```

```
extract_alignments_from_IDs.py
Extract subset of xmap based on ID of ref, query or both.
Usage: python extract_alignments_from_IDs.py -i input.xmap -o output_filename --r_id_list ref_ID1 refID2 ... refIDn --q_id_list qryID1 qryID2 ... qryIDn --in_both
```

```
extract_best_aln.py
Extracts the single best alignment[s] from an XMAP. Writes a new XMAP.
Usage: python extract_best_aln.py -i input.xmap [-o (optional output XMAP filename)] [--each_query (flag to extract best alignment for each query)] [--each_ref (flag to extract best alignment for each reference entry)]
```

```
extract_cmap.py
Extract CMAP(s) given a list of IDs. Writes a new CMAP.
Usage: python extract_cmap.py -i input.cmap -l ID1 ID2 ... IDn [-o optional output CMAP filename]
```

```
extract_contig_aln_mols.py
Extract all molecules aligned to a contig. Writes a BNX with relevant molecules.
Usage: python extact_contig_aln_mols.py -m pileup_map_file.MAP -c contig_id --bnx molecule_bnx_file.BNX [--coords (coordinates on contig, formatted as start-end), default all molecules to entire contig] [-o optional output BNX filename]
```

```
xmap_coverage.py
Compute coverage of reference alignments in XMAP file, for all refIDs present and in total
Usage: python xmap_coverage -x input.xmap [-c refID to get coverage for]
```

Coming soon: ``convert_bnx_version.py`` & tools for working with BNX molecules
