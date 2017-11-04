#' Select representatives for a strata
#'
#' @param x A list of uncles at a specific stratum, each of which contains a
#' list of descendent taxids
#' @export
take_first <- function(x){
  lapply(x, head, 1) 
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

#' Substitute names for taxids in strata
#'
#' @param x strata list, as made by \code{uniprot_cousins}
#' @export
as_named_strata <- function(x){
  scinames <- taxid2name(c(unname(unlist(x)), names(x), unlist(sapply(x, names))))
  backbone <- unname(scinames[names(x)])
  x <- lapply(x,
        function(y) {
          uncles <- unname(scinames[names(y)])
          if(length(y) > 0)
            y <- lapply(y, function(z) unname(scinames[as.character(z)]))
          names(y) <- uncles
          y
        }
      )
  names(x) <- backbone
  x
}

#' Print a strata list
#'
#' @param x A named or unnamed strata
#' @export
print_strata <- function(x){
  for(stratum in names(x)){
    cat(sprintf('%s\n', stratum))
    if(length(unlist(x[[stratum]])) == 0){
      cat("  NO REPRESENTATIVE\n")
    } else {
      for(uncle in names(x[[stratum]])){
        # If the uncle has no included children, ignore it
        if(length(x[[stratum]][[uncle]]) > 0)
          cat(sprintf('  %s\n', uncle))
        for(species in x[[stratum]][[uncle]]){
          cat(sprintf('    %s\n', species))
        }
      }
    }
  }
}
