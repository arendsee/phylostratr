find_revenants <- function(d, classifier=classify_by_evalue(1e-5)){
  check_hit_table(d, has_mrca=TRUE, has_ps=TRUE)
  ds <- split(d, d$qseqid)
  is_revenant <- vapply(
    ds,
    FUN.VALUE=TRUE,
    function(x){
      x$has_hit <- classifier(x)
      ps_has_hit <- x %>%
        dplyr::filter(has_hit) %$% ps %>% min
      ps_no_hit <- x %>%
        dplyr::group_by(ps) %>%
        dplyr::filter(!any(has_hit)) %$% ps %>% max
      ps_no_hit > ps_has_hit
    }
  )
  names(ds)[is_revenant]
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
