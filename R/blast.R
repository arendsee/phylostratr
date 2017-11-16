# This is probably a bit unstable, and certainly hacky. If NCBI decides to
# change the warning message, this will break.
#
# Example:
#
# $ blastdbcmd -info -db blastdb/b.faa
# Database: blastdb/b.faa
#         87 sequences; 26,545 total residues
#
# Date: Nov 15, 2017  4:43 PM     Longest sequence: 2,287 residues
#
# Volumes:
#         /path/to/blastdb/b.faa
#
# A failure looks like:
#
# $ blastdbcmd -info -db waldo.faa
# BLAST Database error: No alias or index file found for protein database [waldo.faa] in search path ...
.blastdb_exists <- function(db){
  suppressWarnings(
    msg <- system2(
      'blastdbcmd',
      stdout = TRUE,
      stderr = TRUE,
      args   = c('-info', '-dbtype', 'prot', '-db', db)
    )
  )
  any(grepl('([0-9,]+) sequences; ([0-9,]+) total residues', msg, perl=TRUE))
}

#' Build a blast database for one species
#'
#' @param fastafile The path to a protein FASTA file
#' @param blastdb The directory in which the blast database should be stored
#' @return path to blast database (without extension)
make_blast_database <- function(
  fastafile,
  blastdb = 'blastdb',
  verbose = FALSE
){
  out <- file.path(blastdb, basename(fastafile))
  if(!.blastdb_exists(out)){
    maybe_message("%s: making blast database ...", verbose, basename(fastafile))
    dbmsg <- system2(
      'makeblastdb',
      stderr = TRUE,
      stdout = TRUE,
      args   = c('-dbtype', 'prot', '-in', fastafile, '-out', out)
    )
    maybe_message(dbmsg, verbose)
    if(!.blastdb_exists(out)){
      stop(sprintf("Failed to make blast database %s", out))
    }
  }
  out
}

#' BLAST query protein FASTA file against a subject species 
#'
#' @param query_fastafile A protein FASTA file for the focal species
#' @param subject_taxid The subject NCBI taxon ID
#' @param blastdb A path to a blast database (as returned from \code{make_blast_database})
#' @param blastresult The output TAB-delimited result file
#' @param nthreads Number of threads
#' @param seg Whether to mask the query protein sequences
#' @return The path to the tabular BLAST result output 
run_blastp <- function(
  query_fastafile,
  subject_taxid,
  blastdb,
  blastresult  = paste0(subject_taxid, ".tab"),
  nthreads     = 1,
  seg          = FALSE,
  verbose      = TRUE
){
  if(file.exists(blastresult)){
    maybe_message("Skipping %s", verbose, subject_taxid)
  } else {
    maybe_message("%s: blasting ...", verbose, subject_taxid)
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
      dplyr::mutate(staxid = subject_taxid) %>%
      dplyr::select(qseqid, staxid, evalue, score) %>%
      readr::write_tsv(blastresult)
  }
  blastresult
}

#' Blast strata 
#'
#' @param query The protein FASTA file for the focal species
#' @param strata named list of phylostrata, where each element is a vector of cousin protein FASTA filenames
#' @param makedb_args Additional arguments passed to \code{make_blast_database} 
#' @param blast_args Additional arguments passed to \code{run_blastp}
#' @return named list of phylostrata, where each element is a vector of blast result filenames 
strata_blast <- function(query, strata, makedb_args=list(), blast_args=list()){
  lapply(
    strata,
    function(cousins){
      vapply(
        FUN.VALUE = "",
        names(cousins),
        function(cousin_name){
          fastafile <- cousins[[cousin_name]]
          blastdb <- do.call(make_blast_database, args=append(fastafile, makedb_args))
          blast_args <- append(list(query, cousin_name, blastdb=blastdb), blast_args)
          do.call(run_blastp, args=blast_args)
        }
      )
    }
  )
}

#' Load each blast result and filter out the best hit against each query gene
#'
#' @param blast_strata A named list of phylostrata where each element consists
#' of a possibly empty list of filenames. The filenames are raw BLAST results.
#' @return A list of lists of data.frames, where each data.frame is a filtered blast result
strata_besthits <- function(blast_strata){
  # produce an empty blast result
  empty_blast_result <- data.frame(
    qseqid = character(0),
    staxid = integer(0),
    evalue = numeric(0),
    score  = numeric(0)
  )
  lapply(blast_strata,
    function(blast_stratum) {
      if(length(blast_stratum) == 0){
        list(empty_blast_result)
      } else {
        lapply(
          blast_stratum,         
          function(blast_file){
            if(length(blast_file) > 0){
              readr::read_tsv(blast_file, col_types='ccdd') %>% get_max_hit
            } else {
              empty_blast_result
            }
          }
        )
      }
    }
  )
}

#' Build a single data.frame with an MRCA column from stratified blast results
#'
#' @param besthits_strata A named list of lists of dataframes. The names are
#' the phylostrata names. The data.frames are filtered BLAST results.
#' @return A single dataframe holding the top hits of each focal gene against
#' each subject species.
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
