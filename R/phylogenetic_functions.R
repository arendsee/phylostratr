backbone <- function(taxid, byname=FALSE){
  taxize::classification(taxid, db='ncbi')[[1]]$id 
}

cousins <- function(taxid){
  taxize::downstream(backbone(taxid)[-1], downto='species', db='ncbi')
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
