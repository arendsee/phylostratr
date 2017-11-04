#' Select representatives for a strata
#'
#' @export
take_first <- function(cousin_sets){
  lapply(cousin_sets, head, 1) 
}

#' Make a filter
#'
#' @param n The number of taxa that must be present before resorting to
#' \code{fun}
#' @param fun The filter function to use for strata with more than \code{n}
#' representatives
#' @export
#' @return A filter that can be used in \code{uniprot_cousin_genomes}
make_do_if_over <- function(n=3, fun=take_first){
  function(cousin_sets){
    ncousins <- sum(vapply(FUN.VALUE=integer(1), cousin_sets, length))
    if(ncousins > n){
      fun(cousin_sets)
    } else {
      cousin_sets
    }
  }
}
