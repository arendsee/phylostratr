#' Make a count matrix from two strata
#'
#' @param d data.frame with two columns holding stratum labels
#' @return CountMatrix object
#' @export
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
