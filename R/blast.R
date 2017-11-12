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

# produce an empty blast result
empty_blast_result <- function(){
  data.frame(
    qseqid = character(0),
    staxid = integer(0),
    evalue = numeric(0),
    score  = numeric(0)
  )
}

#' BLAST query protein FASTA file against a target species 
#'
#' @param query_fastafile A protein FASTA file for the focal species
#' @param target_taxid The target NCBI taxon ID
#' @param blastdb A path to a blast database (as returned from \code{make_blast_database})
#' @param blastresult The output TAB-delimited result file
#' @param nthreads Number of threads
#' @param seg Whether to mask the query protein sequences
#' @return The path to the output file
run_blastp <- function(
  query_fastafile,
  target_taxid,
  blastdb,
  blastresult  = paste0(target_taxid, ".tab"),
  nthreads     = 1,
  seg          = FALSE
){
  if(file.exists(blastresult)){
    message(sprintf("Skipping %s", target_taxid))
  } else {
    message(sprintf("%s: blasting ...", target_taxid))
    system2(
      'blastp',
      stdout=blastresult,
      args=c(
        '-db', blastdb,
        '-query', query_fastafile,
        '-outfmt', '"6 qseqid evalue score"',
        '-num_threads', nthreads,
        '-seg', if(seg) {'yes'} else {'no'}
      )
    )
    # Add the subject taxon ID, name and order columns, write with header
    readr::read_tsv(
      blastresult,
      col_names = c('qseqid', 'evalue', 'score'),
      col_types = 'cdd'
    ) %>%
      dplyr::mutate(staxid = target_taxid) %>%
      dplyr::select(qseqid, staxid, evalue, score) %>%
      readr::write_tsv(blastresult)
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

#' Load each blast result and filter out the best hit against each query gene
#'
#' @param blast_strata The output of strata_blast
#' @return A list of lists of filtered blast results
strata_besthits <- function(blast_strata){
  lapply(blast_strata,
    function(blast_stratum) {
      if(length(blast_stratum) == 0){
        list(empty_blast_result())
      } else {
        lapply(
          blast_stratum,         
          function(blast_file){
            if(length(blast_file) > 0){
              readr::read_tsv(blast_file, col_types='cidd') %>% get_max_hit
            } else {
              empty_blast_result()
            }
          }
        )
      }
    }
  )
}

#' Build a single data.frame with an MRCA column from stratified blast results
#'
#' @param besthits_strata The output of \code{strata_besthits}
merge_besthits <- function(besthits_strata){
  strata <- names(besthits_strata)
  ps <- 1:length(strata)
  lapply(
    ps,
    function(i) {
      do.call(rbind, besthits_strata[[i]]) %>% {
        if(nrow(.) > 0){
          .$mrca <- strata[i]
          .$ps <- i
        } else {
          .$mrca <- integer(0)
          .$ps <- integer(0)
        }
        as.data.frame(.)
      }
    }
  ) %>%
    do.call(what=rbind) %>%
    {
      d <- .

      mrca_map <- d %>%
        dplyr::select(staxid, mrca, ps) %>%
        dplyr::distinct()

      d %>%
        dplyr::select(-mrca, -ps) %>%
        tidyr::complete(qseqid, staxid) %>%
        merge(mrca_map, by='staxid')
    }
}
