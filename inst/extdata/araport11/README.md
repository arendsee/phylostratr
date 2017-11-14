# Araport11 blast results subset 

The `*.faa` files are protein FASTA files from the Araport11 annotation of
*Arabidopsis thaliana*.

For these example files, only a subset of the proteins are used. The sequences
were filtered either on the pattern `AT.G..9` (1/10 of the dataset) or on
`AT.G.99` (1/100).

 1. `at_ATxGxx9xx.faa` - 1/10 of the Araport11 proteins
 2. `at_ATxGx99xx.faa` - 1/100 of the Araport11 proteins
 3. `at_reverse_ATxGxx9xx.faa` - 1/10 of the proteins with the protein sequences reversed
 4. `at_reverse_ATxGx99xx.faa` - 1/100 of the proteins with the protein sequences reversed

The `araport11.Rda` contains the full BLAST result or `Arabidopsis thaliana`
against select uniprot relatives.
