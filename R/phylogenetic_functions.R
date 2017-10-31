check_blast_table <- function(d, has_mrca=FALSE, has_ps=FALSE){
  is_good <-
    is.data.frame(d) &&
    all(c('qseqid', 'evalue', 'bitscore', 'staxid') %in% names(d)) &&
    (!has_mrca || 'mrca' %in% names(d)) &&
    (!has_ps   || 'ps'   %in% names(d))

  if(!is_good){
    stop("Bad blast table")
  }
}

check_taxids <- function(d){
  check_blast_table(d)
  taxids <- unique(d$staxids)
  if(!all(taxize::itis_taxrank(taxids) == 'Species')){
    stop("Expected all subject taxid ids to refer to species, but they don't")
  }
  
}

check_noise <- function(noise){
  # FIXME: stub
}

add_mrca_and_ps <- function(d, qtaxid){
  check_blast_table(d)
  # FIXME: stub
}

find_revenants <- function(d){
  check_blast_table(d, has_mrca=TRUE, has_ps=TRUE)
  # FIXME: stub
}

find_strange <- function(d){
  check_blast_table(d, has_mrca=TRUE, has_ps=TRUE)
  # FIXME: stub
}

constrict <- function(d, on){
  check_blast_table(d, has_mrca=TRUE, has_ps=TRUE)
  stopifnot(on %in% d$mrca)
  # FIXME: stub
}

make_noise <- function(f, r){
  check_blast_table(f)
  check_blast_table(r)
  # FIXME: stub
}

add_pvalue <- function(d, noise){
  check_blast_table(d)
  check_noise(noise)
  # FIXME: stub
}

predict_phylostrata <- function(d, noise){
  # FIXME: stub
}
