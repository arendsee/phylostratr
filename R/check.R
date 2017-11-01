check_hit_table <- function(d, has_mrca=FALSE, has_ps=FALSE){
  is_good <-
    is.data.frame(d) &&
    all(c('qseqid', 'evalue', 'bitscore', 'staxid') %in% names(d)) &&
    (!has_mrca || 'mrca' %in% names(d)) &&
    (!has_ps   || 'ps'   %in% names(d))
  if(nrow(d) == 0){
    warning("The hit table is empty")
  }
  if(!is_good){
    stop("Bad hit table")
  }
}

check_taxids <- function(d){
  check_hit_table(d)
  taxids <- unique(d$staxids)
  if(!all(taxize::itis_taxrank(taxids) == 'Species')){
    stop("Expected all subject taxid ids to refer to species, but they don't")
  }
}
