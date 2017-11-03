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
#' @export
#' @examples
#' \dontrun{
#' # uniprot_downstream_ids(3701) %>% uniprot_retrieve_genomes
#' }
uniprot_retrieve_genomes <- function(taxa, keep_isoforms=FALSE, dir='uniprot-seqs'){
  if(!dir.exists(dir)){
    dir.create(dir, recursive=TRUE)
  }
  inc_str <- if(keep_isoforms){ 'include=yes' } else { 'include=no' }
  for_str <- 'format=fasta'
  for(taxid in taxa){
    fastafile <- file.path(dir, paste0(taxid, '.faa'))
    if(file.exists(fastafile)){
      message(sprintf("Skipping %s - already retrieved", taxid))
    } else {
      message(sprintf("Retrieving %s ...", taxid))
      url_str <- glue::glue(
        "http://www.uniprot.org/uniprot/?query=organism:{taxid}&{for_str}&{inc_str}"
      )
      curl::curl_download(url_str, fastafile)
    }
  }
}

#' Given a focal taxid, find the uniprot uncle representatives
#'
#' @param taxid The focal species NCBI taxon id
#' @param ... Additional arguments sent to \code{uniprot_downstream_ids}
#' @return list of lists of id vectors
uniprot_uncle_ids <- function(taxid, ...){
  us <- uncles(taxid)
  for(ancestor in names(us)){
    taxa <- us[[ancestor]]
    us[[ancestor]] <- lapply(taxa, uniprot_downstream_ids)
    names(us[[ancestor]]) <- taxa
  }
  us
}

#' Retrive the sequences from the uniprot uncle ids
#'
#' @param uncle_id The output of \code{uniprot_uncle_ids}
#' @param dir Directory in which to write sequences
#' @param prefix The phylostratum-level prefix
uniprot_uncle_genomes <- function(uncle_list, dir='strata', prefix='ps_', ...){
  strata <- names(uncle_id)
  for(ps in seq_along(uncle_id)){
    stratum_str <- strata[ps]
    for(node_id in names(uncle_id[[ps]])){
      psdir <- file.path(dir, paste0(prefix, ps), node_id)
      if(!dir.exists(psdir)){
        dir.create(psdir, recursive=TRUE)
      }
      message(sprintf("Retrieving descendents of '%s' ...", node_id))
      for(taxid in uncle_id[[ps]][[node_id]]){
        uniprot_retrieve_genome(node_id, dir=psdir, ...)
      }
    }
  }
}
