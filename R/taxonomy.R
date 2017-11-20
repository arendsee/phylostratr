#' Get the lineage of a given species
#'
#' @param taxid a single NCBI taxon
#' @param byname Get ancestor names, instead of NCBI ids
#' @export
ncbi_ancestors <- function(taxid, byname=FALSE){
  taxize::classification(taxid, db='ncbi')[[1]]$id 
}

#' Get the species level cousins of a taxid
#'
#' @param taxid a single NCBI taxon
#' @export
ncbi_cousins <- function(taxid){
  taxize::downstream(ncbi_ancestors(taxid)[-1], downto='species', db='ncbi')
}

# turn a flat list into a nested tree, with one element per level
.unflatten <- function(x){
  if(length(x) == 1){
    list(
      name = x[1]
    )
  } else if(length(x) > 0){
    list(
      name = x[1],
      .unflatten(x[-1])
    )
  }
}

#' Get all uncles of a species
#'
#' @param taxid a single NCBI taxon
#' @export
ncbi_uncles <- function(taxid){

  # FIXME: cannot find root (taxize issue #639)
  #        so I remove the first index (root)
  tree <- ncbi_ancestors(taxid)[-1] %>%
    .unflatten %>%
    data.tree::FromListSimple()

  tree$Do(function(node) {
    children <- setdiff(
      taxize::children(node$name, db='ncbi')[[1]]$childtaxa_id,
      sapply(node$children, function(n) n$name)
    )
    for(child in children){
      node$AddChild(child)
    }
  })

  tree
}

#' Transform list of taxids to list of taxid summaries
#'
#' @param taxids Vector of NCBI taxonomy ids
#' @param chunk_size number of records to retrive at a time (the default should be fine)
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
  sapply(
    # FUN.VALUE=character(1),
    id_summaries,
    function(x) {
      if('scientificname' %in% names(x)){
        x$scientificname
      } else {
        NA_character_
      }
    }
  )
}
