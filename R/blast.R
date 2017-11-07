#' Build a blast database for one species
#'
#' @param taxid A single NCBI taxon id
#' @param indir The directory in which to find the protein FASTA sequence
#' (should contain a file '<taxid>.faa')
#' @param outdir The directory in which the blast database should be stored
#' @return path to blast database (without extension)
make_blast_database <- function(
  taxid,
  indir  = 'uniprot-seqs',
  outdir = 'blastdb'
){
  fastafile <- file.path(indir, paste0(taxid, ".faa"))
  out <- file.path(outdir, taxid)
  if(!file.exists(paste0(out, ".psq"))){
    message(sprintf("%s: making blast database ...", taxid))
    system2('makeblastdb', args=c('-dbtype', 'prot', '-in', fastafile, '-title', taxid, '-out', out))
  }
  out
}

#' BLAST query protein FASTA file against a target species 
#'
#' @param query_fastafile A protein FASTA file for the focal species
#' @param target_taxid The target NCBI taxon ID
#' @param blastdb The directory in which the blast database should be made
#' @param blastresult The output TAB-delimited result file
#' @param nthreads Number of threads
#' @param seg Whether to mask the query protein sequences
#' @return The path to the output file
run_blastp <- function(
  query_fastafile,
  target_taxid,
  blastdb     = 'blastdb',
  blastresult = paste0(target_taxid, ".tab"),
  nthreads    = 1,
  seg         = FALSE
){
  if(file.exists(blastresult)){
    message(sprintf("Skipping %s", target_taxid))
  } else {
    message(sprintf("%s: blasting ...", target_taxid))
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
  }
  blastresult
}

#' Blast uncle-level strata 
#'
#' Takes a uncle-level strata, blasts each leaf, returns a flat strata with blast result filenames
#'
#' @param query The protein FASTA file for the focal species
#' @param strata Uncle level strata
#' @param makedb_args Additional arguments passed to \code{make_blast_database} 
#' @param blast_args Additional arguments passed to \code{run_blastp}
strata_blast <- function(query, strata, makedb_args=list(), blast_args=list()){
  lapply(
    strata,
    function(uncles){
      taxids <- unname(unlist(uncles))
      vapply(
        FUN.VALUE="",
        taxids,
        function(taxid){
          blastdb <- do.call(make_blast_database, args=append(taxid, makedb_args))
          blast_args <- append(list(query, taxid, blastdb=blastdb), blast_args)
          do.call(run_blastp, args=blast_args)
        }
      )
    }
  )
}
