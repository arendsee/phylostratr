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
#' @param verbose Print progress messages
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

#' Read a blast result of the form expected by phylostratr
#'
#' @param x filename
#' @param with_taxid If TRUE, expect the staxid column to be in the table
#' @export
read_blast <- function(x, with_taxid=TRUE){
  col_types = readr::cols(
    qseqid = readr::col_character(),
    sseqid = readr::col_character(),
    qstart = readr::col_integer(),
    qend   = readr::col_integer(),
    sstart = readr::col_integer(),
    send   = readr::col_integer(),
    evalue = readr::col_double(),
    score  = readr::col_double()
  )
  if(with_taxid){
    col_types[['staxid']] <- readr::col_character()
  }
  readr::read_tsv(x, col_names=names(col_types$cols), col_types=col_types)
}

#' BLAST query protein FASTA file against a subject species 
#'
#' @param query_fastafile A protein FASTA file for the focal species
#' @param subject_taxid The subject NCBI taxon ID
#' @param blastdb A path to a blast database (as returned from \code{make_blast_database})
#' @param blastresult The output TAB-delimited result file
#' @param nthreads Number of threads
#' @param seg Whether to mask the query protein sequences
#' @param verbose Print progress messages
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
        '-outfmt', '"6 qseqid sseqid qstart qend sstart send evalue score"',
        '-num_threads', nthreads,
        '-seg', if(seg) {'yes'} else {'no'}
      )
    )
    # Add the subject taxon ID, name and order columns, write with header
    read_blast(blastresult, with_taxid=FALSE) %>%
      dplyr::mutate(staxid = as.character(subject_taxid)) %>%
      readr::write_tsv(path=blastresult)
  }
  blastresult
}

#' Blast strata 
#'
#' @param strata Strata object where the 'faa' vector is included in the data slot
#' @param makedb_args Additional arguments passed to \code{make_blast_database} 
#' @param blast_args Additional arguments passed to \code{run_blastp}
#' @return named list of phylostrata, where each element is a vector of blast result filenames 
#' @export
strata_blast <- function(
  strata,
  makedb_args=list(),
  blast_args=list()
){
  is_valid_strata(strata, required='faa')

  query <- strata@data$faa[[strata@focal_species]]

  strata@data$blast_result <- lapply(names(strata@data$faa), function(taxid){
    fastafile <- strata@data$faa[[taxid]]
    blastdb <- do.call(make_blast_database, args=append(fastafile, makedb_args))
    blast_args <- append(list(query, taxid, blastdb=blastdb), blast_args)
    do.call(run_blastp, args=blast_args)
  })
  names(strata@data$blast_result) <- names(strata@data$faa)
  strata
}

#' Get the "besthit" table for a given species
#'
#' @param strata Strata object
#' @param taxid A species name or taxonomy ID
#' @return data.frame of best hits (\code{get_max_hit})
#' @export
#' @examples
#' \dontrun{
#' strata@data$besthits[["Unicorn"]] <- get_besthit(strata, "Unicorn")
#' }
get_besthit <- function(strata, taxid){
  is_valid_strata(strata, required=c('faa', 'blast_result'))
  blast_file <- strata@data$blast_result[[taxid]]
  readr::read_tsv(blast_file) %>% get_max_hit
}

#' Load each blast result and filter out the best hit against each query gene
#'
#' @param strata Strata object with 'blast_result' vector in data
#' of a possibly empty list of filenames. The filenames are raw BLAST results.
#' @return Strata object with 'besthits' field in data slot. This field holds a
#' data.frame for each target species, where each data.frame is a filtered
#' blast result
#' @export
strata_besthits <- function(strata){
  is_valid_strata(strata, required=c('faa', 'blast_result'))

  # produce an empty blast result
  empty_blast_result <- data.frame(
    qseqid = character(0),
    staxid = character(0),
    evalue = numeric(0),
    score  = numeric(0)
  )
  taxa <- names(strata@data$blast_result)
  strata@data$besthit <- lapply(taxa, get_besthit, strata=strata)
  names(strata@data$besthit) <- taxa
  strata
}

#' Build a single data.frame with an MRCA column from stratified blast results
#'
#' @param strata A Strata object with a list of dataframe as the data$besthit slot.
#' @return A single dataframe holding the top hits of each focal gene against
#' each subject species.
#' @export
merge_besthits <- function(strata){
  is_valid_strata(strata, required='besthit')

  strata_fold(strata, function(s){
    do.call(what=rbind, s@data$besthit)
  }) %>%
  tuplify %>%
  lapply(function(x){
    x$value$mrca <- x$name
    x$value$ps   <- x$position
    x$value
  }) %>%
  do.call(what=rbind) %>%
  dplyr::mutate(ps = as.integer(.data$ps))
}
