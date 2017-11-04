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

#' Map a function over a specific stratum
#'
#' @param cousin_sets A nested list of cousins
#' @param stratum_id Stratum NCBI taxonomy id
#' @param fun A function of the uncle list
#' @export
do_on <- function(cousin_sets, stratum_id, fun){
  stratum_id <- as.character(stratum_id)
  cousin_sets[[stratum_id]] <- fun(cousin_sets[[stratum_id]]) 
  cousin_sets
}

#' Add an id to a representative list
#'
#' @param cousin_sets A nested list of cousins
#' @param stratum_id Stratum NCBI taxonomy id
#' @param uncle_id Uncle NCBI taxonomy id
#' @param new_id The ID to be added
#' @export
add_to <- function(cousin_sets, stratum_id, uncle_id, new_id){
  stratum_id <- as.character(stratum_id)
  uncle_id <- as.character(uncle_id)
  cousin_sets[[stratum_id]][[uncle_id]] <- append(cousin_sets[[stratum_id]][[uncle_id]], new_id)
  cousin_sets
}

#' Clear all representatives of a given uncle
#'
#' @param cousin_sets A nested list of cousins
#' @param stratum_id Stratum NCBI taxonomy id
#' @param uncle_id Uncle NCBI taxonomy id
#' @export
clear_uncle <- function(cousin_sets, stratum_id, uncle_id){
  stratum_id <- as.character(stratum_id)
  uncle_id <- as.character(uncle_id)
  cousin_sets[[stratum_id]][[uncle_id]] <- integer(0)
  cousin_sets
}

#' Clear all representatives of a given stratum
#'
#' @param cousin_sets A nested list of cousins
#' @param stratum_id Stratum NCBI taxonomy id
#' @export
clear_stratum <- function(cousin_sets, stratum_id){
  cousin_sets[[as.character(stratum_id)]] <- vapply(FUN.VALUE=integer(0), function(x) integer(0))
  cousin_sets
}

#' Add representatives to the strata
#'
#' @param cousin_sets A nested list of cousins
#' @param scinames scientific names of representatives to add
#' @param taxids NCBI ids of representatives to add
#' @export
add_representative <- function(cousin_sets, scinames=NULL, taxids=NULL){
  if(!is.null(scinames)){
    taxids <- taxize::get_uid(scinames, verbose=FALSE) %>% as.character
  }
  lineages <- taxize::classification(taxids, db='ncbi')

  taxids <- as.character(taxids)

  for(taxid in names(lineages)){
    lineage <- lineages[[taxid]]
    stratum <- lineage$id[lineage$id %in% names(cousin_sets)] %>% tail(1)
    uncle <- lineage$id[which(lineage$id %in% stratum) + 1]
    cousin_sets <- add_to(cousin_sets, stratum, uncle, taxid)
  }

  cousin_sets
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
