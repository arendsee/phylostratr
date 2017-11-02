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

check_noise <- function(noise){
  # FIXME: stub
}

add_mrca_and_ps <- function(d, qtaxid){
  check_hit_table(d)
  # FIXME: stub
}

find_revenants <- function(d){
  check_hit_table(d, has_mrca=TRUE, has_ps=TRUE)
  # FIXME: stub
}

find_strange <- function(d){
  check_hit_table(d, has_mrca=TRUE, has_ps=TRUE)
  # FIXME: stub
}

constrict <- function(d, on){
  check_hit_table(d, has_mrca=TRUE, has_ps=TRUE)
  stopifnot(on %in% d$mrca)
  # FIXME: stub
}

make_noise <- function(f, r){
  check_hit_table(f)
  check_hit_table(r)
  # FIXME: stub
}

add_pvalue <- function(d, noise){
  check_hit_table(d)
  check_noise(noise)
  # FIXME: stub
}

predict_phylostrata <- function(d, noise){
  # FIXME: stub
}
