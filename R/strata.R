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

#' Get a map from taxid id to MRCA and stratum level
#'
#' @param strata A three-level list. The first level has one element for each
#' node in the focal species lineage (and is named accordingly). The second
#' level has one element for each 'uncle'. The third level is a possibly empty
#' vector of taxon IDs.
#' @export
strata2mrca <- function(strata){
  lapply(
    seq_along(strata),
    function(i) {
      taxa <- unlist(strata[[i]])
      data.frame(
        taxid = taxa,
        mrca  = rep(names(strata)[i], length(taxa)),
        ps    = rep(i, length(taxa))
      )
    }
  ) %>% do.call(what=rbind) %>% { rownames(.) <- NULL; . }
}


#' Infer homology inference based on a hard e-value threshold
#'
#' @param theshold An e-value threshold
#' @return A function of a dataframe that has the column 'evalue', the function
#' returns a logical vector
classify_by_evalue <- function(threshold){
  function(x){
    !is.na(x$evalue) & (x$evalue < threshold)
  }
}

#' Get an ordered factor mapping MRCA taxon IDs (as vector names) to names
#'
#' @param hittable A data.frame with the columns 'ps' and 'mrca'
#' @export
get_mrca_names <- function(hittable){
  hittable %>%
    dplyr::select(ps, mrca) %>%
    dplyr::distinct() %>%
    dplyr::arrange(ps) %>%
    { .$mrca } %>%
    taxid2name %>%
    { factor(., levels=unname(.)) }
}

#' Get the phylostratum of each query gene
#' 
#' @param hittable data.frame with the columns 'qseqid', 'mrca', 'ps', and
#' whichever columns arerequired by 'classifier'
#' @param classifier A function to infer homology between the query gene and a specific subject
#' @param strata_name A factor of MRCA names (scientific names by default, but
#' you can set your own labels) with levels ordered by phylostrata. This must
#' be a named factor, with names corresponding to taxon IDs. It is used to make
#' the \code{mrca_name} column. It will be used mostly in plots and reports;
#' internally the taxon ids are used.
#' @export
stratify <- function(
    hittable,
    classifier   = classify_by_evalue(1e-5),
    strata_names = get_mrca_names(hittable)
  ){
  hittable %>%
    dplyr::mutate(has_homolog = classifier(hittable)) %>%
    dplyr::filter(has_homolog) %>%
    dplyr::select(qseqid, mrca, ps) %>%
    dplyr::group_by(qseqid) %>%
    dplyr::filter(ps == min(ps)) %>%
    dplyr::ungroup() %>%
    {
      strata <- .
      rbind(
        strata,
        hittable[!(hittable$qseqid %in% strata$qseqid), ] %>%
        dplyr::select(qseqid, mrca, ps) %>%
        dplyr::group_by(qseqid) %>%
        dplyr::filter(ps == max(ps)) %>%
        dplyr::ungroup()
      )
    } %>%
    dplyr::distinct() %>%
    dplyr::mutate(mrca_name = strata_names[mrca])
}

