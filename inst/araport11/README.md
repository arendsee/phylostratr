# Araport11 blast results subset 

## `araport11_subset9.tab`

This is a subset of the results of the following BLAST against the RefSeq release 84 database.

``` sh
#!/bin/bash
# perfomrs NR blast (blastp)
infile="$1"
outfile="$(basename "${infile%.*}").out"
database="/work/GIF/databases/ncbi_nr/nr"
module load ncbi-blast
blastp \
 -query "${infile}" \
 -db "${database}" \
 -out "${outfile}" \
 -num_threads 16 \
 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids salltitles qcovs"
```

The columns `qseqid`, `evalue`, `bitscore`, and `staxid` have been extracted
and a header has been added. A data were subset with the regular expression
`AT.G..9`, which extracted 1/10 loci.
