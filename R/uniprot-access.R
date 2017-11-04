#' Get the uniprot ids downstream of a node
#'
#' @param taxid Ancestral node of clade of interest
#' @param reference_only Should only reference proteomes be considered?
#' @return A numeric vector of NCBI taxon ids listing all species in this clade
#' for which Uniprot has a complete proteome.
#' @export
uniprot_downstream_ids <- function(taxid, reference_only=FALSE){
  ref_str <- if(reference_only){
    'reference:yes'
  } else {
    'complete:yes'
  }
  url_str <- glue::glue(
    'http://www.uniprot.org/taxonomy/?query=ancestor:{taxid}+{ref_str}&format=list'
  )
  con <- curl::curl(url_str)
  children <- readLines(con) %>% as.integer
  close(con)
  children
}

#' Download UniProt proteomes for given taxa
#'
#' Each proteom is written to a file with the name '<taxid>.faa', for example,
#' Arabidopsis thaliana, which has the id 3702, will be written the '3702.faa'.
#'
#' @param taxa Vector of NCBI taxonomy ids
#' @param keep_isoforms Should all isoforms of each protein be included?
#' @param dir Directory in which to write all FASTA files
#' @param dryrun If TRUE, do not download genomes, just print the URLs and destination files
#' @export
#' @examples
#' \dontrun{
#' # uniprot_downstream_ids(3701) %>% uniprot_retrieve_genomes
#' }
uniprot_retrieve_genomes <- function(taxa, keep_isoforms=FALSE, dir='uniprot-seqs', dryrun=FALSE){
  if(!dir.exists(dir) && !dryrun){
    dir.create(dir, recursive=TRUE)
  }
  inc_str <- if(keep_isoforms){ 'include=yes' } else { 'include=no' }
  for_str <- 'format=fasta'
  for(taxid in taxa){
    fastafile <- file.path(dir, paste0(taxid, '.faa'))
    if(file.exists(fastafile)){
      message(sprintf("Skipping %s - already retrieved", taxid))
    } else {
      url_str <- glue::glue(
        "http://www.uniprot.org/uniprot/?query=organism:{taxid}&{for_str}&{inc_str}"
      )
      if(dryrun){
        message(sprintf("Checking %s ...", taxid))
        message(sprintf("  url: %s", url_str))
        message(sprintf("  destination: %s", fastafile))
      } else {
        message(sprintf("Retrieving %s ...", taxid))
        curl::curl_download(url_str, fastafile)
      }
    }
  }
}

#' Given a focal taxid, find the uniprot descendents of a taxon's uncles 
#'
#' @param taxid The focal species NCBI taxon id
#' @param ... Additional arguments sent to \code{uniprot_downstream_ids}
#' @return list of lists of id vectors
uniprot_cousins <- function(taxid, ...){
  us <- uncles(taxid)
  for(ancestor in names(us)){
    taxa <- us[[ancestor]]
    us[[ancestor]] <- lapply(taxa, uniprot_downstream_ids)
    names(us[[ancestor]]) <- taxa
  }
  us
}

#' Retrive the sequences from the uniprot cousins
#'
#' cousin_sets A named list of vectors of species taxon ids. Each vector
#' in the list bears the name of an uncle, all at the same level in the tree.
#' There may be multiple uncles in multifurcating nodes of the tree (this is
#' very common in the NCBI common tree).
#'
#' @param cousins The output of \code{uniprot_cousins}
#' @param dir Directory in which to write sequences
#' @param prefix The phylostratum-level prefix
#' @param dryrun If TRUE, do not download files or create directories
#' @param ... Additional arguments passed to \code{uniprot_retrieve_genomes}
#' @return Nothing, this function is run for its effects
#' @examples
#' \dontrun{
#' uniprot_cousins(3702) %>%
#'   lapply(take_first) %>%
#'   uniprot_cousin_genomes
#' 
#' cfilter <- make_do_if_over(3, take_first)
#' uniprot_cousins(3702) %>%
#'   lapply(cfilter) %>%
#'   uniprot_cousin_genomes
#' }
uniprot_cousin_genomes <- function(cousins, dir='strata', prefix='ps_', dryrun=FALSE, ...){
  strata <- names(cousins)
  for(ps in seq_along(cousins)){
    stratum_str <- strata[ps]
    for(node_id in names(cousins[[ps]])){
      psdir <- file.path(dir, paste0(prefix, ps), node_id)
      if(dryrun){
        message(sprintf("Retrieving descendents of '%s' ...", node_id))
      } else {
        message(sprintf("Checking descendents of '%s' ...", node_id))
        if(!dir.exists(psdir) && length(cousins[[ps]][[node_id]] > 0)){
          dir.create(psdir, recursive=TRUE)
        }
      }
      for(taxid in cousins[[ps]][[node_id]]){
        message(node_id)
        uniprot_retrieve_genomes(node_id, dir=psdir, dryrun=dryrun, ...)
      }
    }
  }
}

#' Transform list of taxids to list of taxid summaries
#'
#' @param taxids Vector of NCBI taxonomy ids
#' @param number of records to retrive at a time (the default should be fine)
#' @export
taxid2summary <- function(taxids, chunk_size=250){
  id_summaries <- as.list(taxids) %>% lapply(function(x) NULL)
  N <- length(taxids)
  for(i in (0:(N %/% chunk_size))){
    i <- i*chunk_size + 1
    j <- min(N, i + chunk_size - 1)
    # If retrieving only one record, the return type is [fields]
    if(i == j){
      id_summaries[[i]] <- rentrez::entrez_summary(id=taxids[i:j], db='taxonomy')
    # Else the return type is [[fields]]
    } else {
      id_summaries[i:j] <- rentrez::entrez_summary(id=taxids[i:j], db='taxonomy')
    }
    if(j < N)
      Sys.sleep(0.3) # so as not to piss off NCBI
  }
  names(id_summaries) <- taxids
  id_summaries
}

#' Get the scientific names for a list of taxids
#'
#' @param taxids A list of NCBI taxonomy ids
#' @param ... Arguments passed to \code{taxid2summary}
#' @export
taxid2name <- function(taxids, ...){
  id_summaries <- taxid2summary(taxids, ...)
  vapply(FUN.VALUE=character(1), id_summaries, function(x) x$scientificname)
}

#' Select representatives for a strata
#'
#' @export
take_first <- function(cousin_sets){
  lapply(cousin_sets, head, 1) 
}

#' Make a filter
#'
#' @param n The number of taxa that must be present before resorting to
#' \code{fun}
#' @param fun The filter function to use for strata with more than \code{n}
#' representatives
#' @export
#' @return A filter that can be used in \code{uniprot_cousin_genomes}
make_do_if_over <- function(n=3, fun=take_first){
  function(cousin_sets){
    ncousins <- sum(vapply(FUN.VALUE=integer(1), cousin_sets, length))
    if(ncousins > n){
      fun(cousin_sets)
    } else {
      cousin_sets
    }
  }
}
