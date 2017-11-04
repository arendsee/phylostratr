#' Get the lineage of a given species
#'
#' @param taxid a single NCBI taxon
#' @param byname Get ancestor names, instead of NCBI ids
#' @export
ancestors <- function(taxid, byname=FALSE){
  taxize::classification(taxid, db='ncbi')[[1]]$id 
}

#' Get the species level cousins of a taxid
#'
#' @param taxid a single NCBI taxon
#' @export
cousins <- function(taxid){
  taxize::downstream(ancestors(taxid)[-1], downto='species', db='ncbi')
}

#' Get all uncles of a species
#'
#' @param taxid a single NCBI taxon
#' @export
uncles <- function(taxid){
  # FIXME: cannot find root (taxize issue #639)
  #        so I remove the first index (root)
  lineage <- ancestors(taxid)[-1]
  children <- taxize::children(lineage, db='ncbi') %>% {.[-length(.)]}
  lineage <- lineage[-1]
  lapply(
    seq_along(children),
    function(i) {
      setdiff(children[[i]]$childtaxa_id, lineage[i]) %>% as.integer
    }
  ) %>%
    magrittr::set_names(names(children))
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
