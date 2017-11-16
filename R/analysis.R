#' Find genes that with matches that skip strata
#'
#' @param d data.frame of maximal hits
#' @param classifier function of a blast result that infers whether each
#' subject species contains a homolog
#' @return data.frame with columns qseqid and ps, where ps represents the
#' strata that are skipped
find_revenants <- function(d, classifier=classify_by_evalue(1e-5)){
  check_hit_table(d, has_mrca=TRUE, has_ps=TRUE)
  focal_ps <- max(d$ps)
  d %>%
    dplyr::mutate(has_hit = classifier(d)) %>%
    dplyr::select(staxid, qseqid, ps, has_hit) %>%
    dplyr::group_by(qseqid) %>%
    dplyr::mutate(basal_ps = min(c(focal_ps, ps[has_hit]))) %>%
    dplyr::group_by(qseqid, ps) %>%
    dplyr::mutate(mrca_has_hit = any(has_hit)) %>%
    dplyr::group_by(qseqid) %>%
    dplyr::filter((ps > basal_ps) & (!mrca_has_hit)) %>%
    dplyr::ungroup() %>%
    dplyr::select(qseqid, ps, basal_ps) %>%
    dplyr::distinct()
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
