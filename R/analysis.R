#' Find genes that with matches that skip strata
#'
#' @param d data.frame of maximal hits
#' @param classifier function of a blast result that infers whether each
#' subject species contains a homolog
#' @return data.frame with columns [qseqid, ps, mrca, basal_ps, n].  Where
#' qseqid is the focal species taxond IDs, ps is the phylostratum level (where
#' 1 is generally cellular_organisms) of strata that are skipped, mrca is the
#' taxon ID of the ancestor, basal_ps is the oldest stratum with a match to
#' qseqid, and n is the total number of revenants for seqid.
#' @export
find_revenants <- function(d, classifier=classify_by_evalue(1e-5)){
  check_hit_table(d, has_mrca=TRUE, has_ps=TRUE)
  focal_ps <- max(d$ps)
  d %>%
    dplyr::mutate(has_hit = classifier(d)) %>%
    dplyr::select(.data$staxid, .data$qseqid, .data$ps, .data$mrca, .data$has_hit) %>%
    dplyr::group_by(.data$qseqid) %>%
    dplyr::mutate(basal_ps = min(c(focal_ps, .data$ps[.data$has_hit]))) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(.data$qseqid, .data$ps) %>%
    dplyr::mutate(mrca_has_hit = any(.data$has_hit)) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(.data$qseqid) %>%
    dplyr::filter((.data$ps > .data$basal_ps) & (!.data$mrca_has_hit)) %>%
    dplyr::ungroup() %>%
    dplyr::select(.data$qseqid, .data$ps, .data$mrca, .data$basal_ps) %>%
    dplyr::distinct() %>%
    dplyr::add_count(.data$qseqid) %>%
    dplyr::arrange(.data$qseqid, .data$ps)
}

find_skip_runs <- function(d, classifier=classify_by_evalue(1e-5)){
  has_hit <- classifier(d) 
  d <- d[has_hit, ]

  k = max(d$ps) + 1 

  x <- d %>%
    dplyr::select(.data$qseqid, .data$ps) %>%
    dplyr::arrange(.data$qseqid, .data$ps) %>%
    dplyr::distinct() %>%
    dplyr::group_by(.data$qseqid) %>%
    dplyr::transmute(
      a = .data$ps[-1],
      b = seq_along(.data$ps)
    )
}

#' Find the focal sequences that are not represented in the given strata
#'
#' @param d data.frame of maximal hits
#' @param on The strata levels (NOT taxon IDs) on which to constrict
#' @param revenants A table of the form created by \code{find_revenants}. If
#' NULL, a new table is created.
#' @param ... Additional arguments sent to \code{find_revenants} if a new
#' revenant table is being created (i.e. 'revenants=NULL').
#' @return A subset of d where no gene has a match to the strata specified in 'on'
#' @export
constrict <- function(d, on, revenants=NULL, ...){
  check_hit_table(d)
  if(is.null(revenants)){
    revenants <- find_revenants(d, ...)
  }
  seqids <- revenants %>%
    dplyr::group_by(.data$qseqid) %>%
    dplyr::mutate(constrict = all(on %in% .data$ps)) %>%
    dplyr::ungroup() %>%
    dplyr::filter(.data$constrict)$qseqid %>% unique
  d[d$qseqid %in% seqids, ]
}

find_strange <- function(d){
  check_hit_table(d, has_mrca=TRUE, has_ps=TRUE)
  # FIXME: stub
}

make_noise <- function(f, r){
  check_hit_table(f)
  check_hit_table(r)
  # FIXME: stub
}
