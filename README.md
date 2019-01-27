## BioNano tools
### Assorted scripts for working with BioNano data.

```
convert_cmap_version.py
Converts between CMAP versions. Outputs a new CMAP
Usage: python convert_cmap_version.py -i input.cmap -v [version to convert to (0.1/0.2)] [-o optional output CMAP filename]
```

```
extract_best_aln.py
Extracts the single best alignment[s] from an XMAP. Writes a new CMAP.
Usage: python extract_best_aln.py -i input.xmap [-o (optional output XMAP filename)] [--each_query (flag to extract best alignment for each query)] [--each_ref (flag to extract best alignment for each reference entry)]
```

```
extract_cmap.py
Extract a single CMAP given ID. Writes a new CMAP.
Usage: python extract_cmap.py -i input.cmap -c contig_id [-o optional output CMAP filename]
```

```
extract_contig_aln_mols.py
Extract all molecules aligned to a contig. Writes a BNX with relevant molecules.
Usage: python extact_contig_aln_mols.py -m pileup_map_file.MAP -c contig_id --bnx molecule_bnx_file.BNX [--coords (coordinates on contig, formatted as start-end), default all molecules to entire contig] [-o optional output BNX filename]
```

Coming soon: ``convert_bnx_version.py`` & tools for working with BNX molecules
