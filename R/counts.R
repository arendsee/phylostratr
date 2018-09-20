# Run a couple comparisons between methods
# Internal
run_comparison <- function(results){
  list(
    b1 = stratify(results, classify_by_adjusted_pvalue(1e-1))[, c('qseqid', 'mrca_name')],
    b2 = stratify(results, classify_by_adjusted_pvalue(1e-2))[, c('qseqid', 'mrca_name')],
    b5 = stratify(results, classify_by_adjusted_pvalue(1e-5))[, c('qseqid', 'mrca_name')]
  )
}

make_matrix_from_two_strata <- function(d){
  # for now, only use the first two columns
  d <- d[, c(1,2)]
  labels <- names(d)
  names(d) <- c("a", "b")

  mat <- dplyr::group_by(d, .data$a, .data$b) %>%
    dplyr::count() %>%
    reshape2::acast(a ~ b, fill=0)

  CountMatrix(x=mat, ylab=labels[1], xlab=labels[2])
}

#' Build a matrix relating 2014 study to 2018
#' 
#' With 2014 classifications on the y axis and 2018 on the x axis.
#'
#' @param strata Strata object
#' @param results data.frame of highest scoring hits of each gene against each target
#' @param classifier function of 'results' that returns a logical vector
#' @return matrix
#' @export
make_TIPS_comparison_matrix <- function(strata, results, classifier){
  # Get the strata from the current analysis
  a2018 <- stratify(hittable=results, classifier=classifier)
  # Get and clean the (Arendsee et al., 2014) data (stored in phylostratr package)
  a2014 <- readr::read_tsv(system.file('extdata', 'arendsee2014_strata.tab', package='phylostratr'))
  # Add NCBI taxonomy ID to tips data
  a2014$mrca <- taxizedb::name2taxid(a2014$name)
  # match column naming conventions
  a2014 <- dplyr::select(a2014, .data$qseqid, .data$mrca, .data$ps, mrca_name=.data$name)
  # underscores to spaces 
  a2014$mrca_name <- sub('_', ' ', a2014$mrca_name)
  # get a lineage map for Arabidopsis
  strata_map <- taxizedb::classification('3702')[[1]]
  # factor the backbone
  a2014$mrca_name <- droplevels(factor(a2014$mrca_name, levels=strata_map$name))

  ps <- standardize_strata(list(
      a2014 = a2014,
      a2018 = a2018
  ))
  ps <- lapply(ps, function(x) x[c('qseqid', 'mrca_name')])

  species2018 <- strata %>% strata_fold(function(s) taxizedb::taxid2name(s@tree$tip.label)) %>%
    { names(.) <- taxizedb::taxid2name(names(.)); . } %>%
    reshape2::melt() %>%
    dplyr::select(phylostratum = .data$L1, species = .data$value)

  species2014 <- readr::read_tsv(system.file('case-studies', 'tips-ages.tab', package='phylostratr'))

  d <- merge(ps[[1]], ps[[2]], by="qseqid", suffixes=paste0(".", names(ps)))

  make_matrix_from_two_strata(d[, -1])
}

#' Get comparisons of inferences across scoring systems
#' 
#' @param m data.frame of mrca
#' @return data.frame of stratum labels
#' @export
make_significance_matrices <- function(m){
  if(ncol(m) < 2){
    stop("Expected a data.frame with 2 or more columns")
  }
  gs <- list()
  for(i in 1:(ncol(m)-1)){
    gss <- list()
    for(j in (i+1):ncol(m)){
      gss[[names(m)[j]]] <- make_matrix_from_two_strata(m[, c(i, j)])
    }
    gs[[names(m)[i]]] <- gss
  }
  gs
}

#' Get comparisons of inferences across scoring systems
#' 
#' @param d data.frame containing a column of group labels and a column of indices
#' @param labels character vector of labels
#' @return matrix
#' @export
make_jump_matrix <- function(d, labels=NULL){
  # expect two columns
  stopifnot(ncol(d) == 2)
  names(d) <- c("group", "index")
  stopifnot(is.numeric(d$index))

  mat <- d %>%
    dplyr::select(.data$group, .data$index) %>%
    dplyr::distinct() %>%
    dplyr::arrange(.data$index) %>%
    dplyr::group_by(.data$group) %>%
    dplyr::transmute(
      from = c(NA, head(.data$index, -1)),
      to   = c(NA, tail(.data$index, -1))
    ) %>%
    dplyr::ungroup() %>%
    dplyr::filter(!is.na(.data$from) & !is.na(.data$to)) %>%
    dplyr::select(.data$from, .data$to) %>%
    dplyr::group_by(.data$from, .data$to) %>%
    dplyr::count() %>%
    dplyr::ungroup() %>%
    reshape2::acast(from ~ to, fill=0) %>% t

  if(!is.null(labels)){
    dimnames(mat) <- list(labels[as.integer(dimnames(mat)[[1]])], labels[as.integer(dimnames(mat)[[2]])])
  }

  CountMatrix(x=mat)
}
