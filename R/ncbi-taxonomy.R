#' Build a phylogenetic tree from a list of NCBI taxon IDs
#'
#' @param taxa vactor of NCBI taxon IDs
#' @return phylo object
#' @export
#' @examples
#' \dontrun{
#' ncbi_tree(c(9606, 9598, 9593))
#' }
ncbi_tree <- function(taxa){
  lineages_to_phylo(taxizedb::classification(taxa, db='ncbi'))
}

#' Get the lineage of a given species
#'
#' @param taxid a single NCBI taxon
#' @param byname Get ancestor names, instead of NCBI ids
#' @export
ncbi_ancestors <- function(taxid, byname=FALSE){
  taxizedb::classification(taxid, db='ncbi')[[1]]$id 
}

#' Get the species level cousins of a taxid
#'
#' @param taxid a single NCBI taxon
#' @export
ncbi_cousins <- function(taxid){
  taxizedb::downstream(ncbi_ancestors(taxid)[-1], downto='species', db='ncbi')
}

#' Get all ancestral sisters of a taxon
#'
#' @param taxid a single NCBI taxon
#' @return phylo object with all NCBI sisters of all ancestors
#' @export
ncbi_aunts <- function(taxid){
  lin <- ncbi_ancestors(taxid)[-1] %>%
    head(-1) %>% tail(-1) %>%
    taxizedb::children()
  lapply(names(lin), function(x){
      kids <- lin[[x]]$childtaxa_id
      matrix(c(rep(x, length(kids)), kids), ncol=2)
    }) %>%
    do.call(what=rbind) %>%
    edgelist_to_phylo
}
