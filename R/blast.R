#' Build a blast database for one species
#'
#' @param taxid A single NCBI taxon id
#' @param indir The directory in which to find the protein FASTA sequence
#' (should contain a file '<taxid>.faa')
#' @param outdir The directory in which the blast database should be stored
make_blast_database <- function(taxid, indir, outdir){
  fastafile <- file.path(indir, paste0(taxid, ".faa"))
  out <- file.path(outdir, taxid)
  system2('makeblastdb', args=c('-dbtype', 'prot', '-in', fastafile, '-title', taxid, '-out', out))
}

#' BLAST query protein FASTA file against a target species 
#'
#' @param query_fastafile A protein FASTA file for the focal species
#' @param target_taxid The target NCBI taxon ID
#' @param blastdb The directory in which the blast database should be made
#' @param blastresult The output TAB-delimited result file
#' @param nthreads Number of threads
#' @param seg Whether to mask the query protein sequences
run_blastp <- function(query_fastafile, target_taxid, blastdb, blastresult, nthreads=1, seg=FALSE){
  system2(
    'blastp',
    stdout=blastresult,
    args=c(
      '-db', file.path(blastdb, target_taxid),
      '-query', query_fastafile,
      '-outfmt', '"6 qseqid evalue score"',
      '-num_threads', nthreads,
      '-seg', if(seg) {'yes'} else {'no'}
    )
  )
  blastresult
}
