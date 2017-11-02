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
