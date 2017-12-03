#' Get a diverse subset of species from tree
#'
#' @param tree phylo object
#' @param n number of species to keep (if greater then the number of species in
#' the tree, the tree is returned unchanged)
#' @param weights Numeric vector with length equal to the number of species in
#' the tree. A weight of 1 will have no influence, lower than one means the
#' species is less likely to be selected.
#' @param collapse Should nodes with a single descendent be collapse?
#' @return phylo object
#' @export
#' @examples
#' data(atree)
#' 
#' # default
#' diverse_subtree(atree, 4)
#'
#' # do not to include the 't1' species
#' diverse_subtree(atree, 4, weights=c(0,1,1,1,1,1,1,1,1,1))
diverse_subtree <- function(tree, n, weights=rep(1, nleafs(tree)), collapse=FALSE){
  if(n < 1){
    stop('Must select at least one species')
  }
  if(length(weights) != nleafs(tree)){
    stop('The number of weights must equal the number of species')
  }

  n <- min(n, nleafs(tree))

  lins <- lapply(1:nleafs(tree), lineage, tree=tree)

  for(i in 1:n){
    if(i == 1){
      # start by taking the species with the highest phylogeny independent weight
      chosen_taxa <- head(which.max(weights), 1)
      # initial list of visited taxa
      seen <- lins[[chosen_taxa]]
    } else {
      # then scale weights to penalize species that are close to the chosen species
      scaled.weights <- weights *
        sapply(lins, function(x){
          1 - length(intersect(x, seen)) / length(x)
        })
      new_taxon <- head(which.max(scaled.weights), 1)
      seen <- union(seen, lins[[new_taxon]])
      chosen_taxa <- c(chosen_taxa, new_taxon)
    }
  }

  subtree(tree, chosen_taxa, collapse=collapse)

}

#' Apply f to each outgroup branch ascending node id
#'
#' @param strata Strata object
#' @param f A function of a phylo object that return a phylo object
#' @param id A name or index
#' @param ... Additional arguments passed to f
#' @return Strata object
strata_apply <- function(strata, f, id=strata@focal_id, ...){
  lin <- lineage(strata@tree, id)[-1]

  outgroups <- lapply(lin, function(ancestor){
    outgroup <- tree_names(strata@tree)[sisters(strata@tree, ancestor)] %>%
      subtree(tree=strata@tree, collapse=FALSE, descend=TRUE) %>%
      f(...)
  })
  names(outgroups) <- tree_names(strata@tree)[lin]
  outgroups <- outgroups[!sapply(outgroups, is.null)]

  final <- lapply(outgroups, tree2edgelist) %>%
    do.call(what=rbind) %>%
    unique %>%
    edgelist_to_phylo %>%
    ape::collapse.singles()

  new_data <- lapply(strata@data, function(x){
    w <- rep(NA, nleafs(final))
    names(w) <- final$tip.label
    common <- intersect(final$tip.label, strata@tree$tip.label)
    w[common] <- strata@data[common]
    w
  })

  new_strata <- strata

  new_strata@data <- new_data
  new_strata@tree <- final

  new_strata
}

#' Map a function over a specific stratum
#'
#' @param cousin_sets A nested list of cousins
#' @param stratum_id Stratum NCBI taxonomy id
#' @param fun A function of the aunt list
#' @export
do_on <- function(cousin_sets, stratum_id, fun){
  ### FIXME: convert to tree
  # stratum_id <- as.character(stratum_id)
  # cousin_sets[[stratum_id]] <- fun(cousin_sets[[stratum_id]])
  # cousin_sets
}

#' Add an id to a representative list
#'
#' @param cousin_sets A nested list of cousins
#' @param stratum_id Stratum NCBI taxonomy id
#' @param aunt_id aunt NCBI taxonomy id
#' @param new_id The ID to be added
#' @export
add_to <- function(cousin_sets, stratum_id, aunt_id, new_id){
  ### FIXME: convert to tree
  # stratum_id <- as.character(stratum_id)
  # aunt_id <- as.character(aunt_id)
  # cousin_sets[[stratum_id]][[aunt_id]] <- append(cousin_sets[[stratum_id]][[aunt_id]], new_id)
  # cousin_sets
}

#' Clear all representatives of a given aunt
#'
#' @param cousin_sets A nested list of cousins
#' @param stratum_id Stratum NCBI taxonomy id
#' @param aunt_id aunt NCBI taxonomy id
#' @export
clear_aunt <- function(cousin_sets, stratum_id, aunt_id){
  ### FIXME: convert to tree
  # stratum_id <- as.character(stratum_id)
  # aunt_id <- as.character(aunt_id)
  # cousin_sets[[stratum_id]][[aunt_id]] <- integer(0)
  # cousin_sets
}

#' Clear all representatives of a given stratum
#'
#' @param cousin_sets A nested list of cousins
#' @param stratum_id Stratum NCBI taxonomy id
#' @export
clear_stratum <- function(cousin_sets, stratum_id){
  ### FIXME: convert to tree
  # cousin_sets[[as.character(stratum_id)]] <- vapply(FUN.VALUE=integer(0), function(x) integer(0))
  # cousin_sets
}

#' Add representatives to the strata
#'
#' @param cousin_sets A nested list of cousins
#' @param scinames scientific names of representatives to add
#' @param taxids NCBI ids of representatives to add
#' @export
add_representative <- function(cousin_sets, scinames=NULL, taxids=NULL){
  ### FIXME: convert to tree
  # if(!is.null(scinames)){
  #   taxids <- taxize::get_uid(scinames, verbose=FALSE) %>% as.character
  # }
  # lineages <- taxize::classification(taxids, db='ncbi')
  #
  # taxids <- as.character(taxids)
  #
  # for(taxid in names(lineages)){
  #   lineage <- lineages[[taxid]]
  #   stratum <- lineage$id[lineage$id %in% names(cousin_sets)] %>% tail(1)
  #   aunt <- lineage$id[which(lineage$id %in% stratum) + 1]
  #   cousin_sets <- add_to(cousin_sets, stratum, aunt, taxid)
  # }
  #
  # cousin_sets
}

#' Print a strata list
#'
#' @param x A named or unnamed strata
#' @export
print_strata <- function(x){
  ### FIXME: convert to tree
  # for(stratum in names(x)){
  #   cat(sprintf('%s\n', stratum))
  #   if(length(unlist(x[[stratum]])) == 0){
  #     cat("  NO REPRESENTATIVE\n")
  #   } else {
  #     for(aunt in names(x[[stratum]])){
  #       # If the aunt has no included children, ignore it
  #       if(length(x[[stratum]][[aunt]]) > 0)
  #         cat(sprintf('  %s\n', aunt))
  #       for(species in x[[stratum]][[aunt]]){
  #         cat(sprintf('    %s\n', species))
  #       }
  #     }
  #   }
  # }
}

#' Get a map from taxid id to MRCA and stratum level
#'
#' @param strata A three-level list. The first level has one element for each
#' node in the focal species lineage (and is named accordingly). The second
#' level has one element for each 'aunt'. The third level is a possibly empty
#' vector of taxon IDs.
#' @export
strata2mrca <- function(strata){
  ### FIXME: convert to tree
  # lapply(
  #   seq_along(strata),
  #   function(i) {
  #     taxa <- unlist(strata[[i]])
  #     data.frame(
  #       taxid = taxa,
  #       mrca  = rep(names(strata)[i], length(taxa)),
  #       ps    = rep(i, length(taxa))
  #     )
  #   }
  # ) %>% do.call(what=rbind) %>% { rownames(.) <- NULL; . }
}


#' Infer homology inference based on a hard e-value threshold
#'
#' @param threshold An e-value threshold
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
    dplyr::select(.data$ps, .data$mrca) %>%
    dplyr::distinct() %>%
    dplyr::arrange(.data$ps) %>%
    { .$mrca } %>%
    taxid2name %>%
    { factor(., levels=unname(.)) }
}

#' Get the phylostratum of each query gene
#' 
#' @param hittable data.frame with the columns 'qseqid', 'mrca', 'ps', and
#' whichever columns arerequired by 'classifier'
#' @param classifier A function to infer homology between the query gene and a specific subject
#' @param strata_names A factor of MRCA names (scientific names by default, but
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
  hittable[classifier(hittable), ] %>%
    dplyr::select(.data$qseqid, .data$mrca, .data$ps) %>%
    dplyr::group_by(.data$qseqid) %>%
    dplyr::filter(.data$ps == min(.data$ps)) %>%
    dplyr::ungroup() %>%
    {
      strata <- .
      rbind(
        strata,
        hittable[!(hittable$qseqid %in% strata$qseqid), ] %>%
        dplyr::select(.data$qseqid, .data$mrca, .data$ps) %>%
        dplyr::group_by(.data$qseqid) %>%
        dplyr::filter(.data$ps == max(.data$ps)) %>%
        dplyr::ungroup()
      )
    } %>%
    dplyr::distinct() %>%
    dplyr::mutate(mrca_name = strata_names[.data$mrca])
}

