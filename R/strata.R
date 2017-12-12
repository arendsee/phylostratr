#' Assert a Strata object is valid, die on failure
#'
#' @param strata Strata object
is_valid_strata <- function(strata){
  if(!(strata@focal_species %in% strata@tree$tip.label)){
    stop("Invalid Strata object, the focal species '", strata@focal_species, "' is not found in the tree")
  }
  for(field in names(strata@data)){
    if(!all(names(strata@data[[field]]) %in% strata@tree$tip.label)){
      stop("Invalid Strata object, data field '", field, "' contains species labels that are not in the tree")
    }
  }
}

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
#' diverse_subtree(atree, 4, weights=c(t1=0))
diverse_subtree <- function(tree, n, weights=NULL, collapse=FALSE){
  if(n < 1){
    stop('Must select at least one species')
  }

  if(is.null(weights)){
    weights <- rep(1, length(tree$tip.label))
  } else {
    weights <- weights[tree$tip.label]
    weights[is.na(weights)] <- 1
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
#' @param ... Additional arguments passed to f
#' @return Strata object
strata_apply <- function(strata, f, ...){
  is_valid_strata(strata)

  lin <- lineage(strata@tree, strata@focal_species, type='name')[-1]

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

#' Add an id to a representative list
#'
#' @param strata A Strata object
#' @param taxid NCBI taxon ID of the species to add
#' @return Strata object
#' @export
add_taxid <- function(strata, taxid){
  if(length(taxid) > 1)
    stop("Currently only one taxid is handled at a time in this 'add_taxid'")
  if(any(taxid %in% tree_names(strata@tree)))
    return(strata)

  # Get new node with its lineage
  new_tree <- taxizedb::classification(taxid)[[1]]$id %>%
    lineage_to_ancestor_tree
  # Find the most recent common node
  max_id <- max(which(tree_names(new_tree) %in% tree_names(strata@tree)))
  # Find the location at which this should attach to the old tree
  common_node_index <- which(tree_names(strata@tree) == tree_names(new_tree)[max_id])
  # Name of parent
  pid <- tree_names(strata@tree)[parent(strata@tree, common_node_index, type='index')]

  # All species descending from the common node need to be reclassified since
  # their lineages may have been collapsed. To resolve the relationship between
  # these species, I need to recover the possibly collapsed nodes and rebuild
  # the tree.
  all_tree <- c(taxid, tree_names(strata@tree)[descendents(strata@tree, common_node_index)]) %>%
    taxizedb::classification() %>%
    lineages_to_phylo

  # Then remove the root of the subtree
  strata@tree <- prune(strata@tree, common_node_index, type='index')

  # Finally, bind the subtree to the stub
  where <- which(tree_names(strata@tree) == pid)
  strata@tree <- ape::bind.tree(strata@tree, all_tree, where=where)

  strata
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
    taxizedb::taxid2name() %>%
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

#' Convert strata data, tip, and node names to and from NCBI taxonomy IDs
#'
#' @param strata Strata object
#' @param target What to convert, may be 'tip', 'node', or 'all'
#' @param to What to convert to, may be 'id' or 'name'
#' @return Strata object with new names
strata_convert <- function(strata, target='tip', to='id'){
  FUN <- switch(
    to,
    id = taxizedb::name2taxid,
    name = taxizedb::taxid2name,
    stop("Illegal 'to' value: must be either 'id' or 'name'")
  )
  if(target == 'tip' || target == 'all'){
    strata@tree$tip.label <- FUN(strata@tree$tip.label)
    for(item in names(strata@data)){
      names(strata@data[[item]]) <- FUN(names(strata@data[[item]]))
    }
    strata@focal_species <- FUN(strata@focal_species)
  }
  if(target == 'node' || target == 'all'){
    strata@tree$node.label <- FUN(strata@tree$node.label)
  }
  strata
}

#' Sort strata relative to focal species
#'
#' @param strata Strata object
#' @return Strata object with reorded tips
#' @export
sort_strata <- function(strata){
  is_valid_strata(strata)
  strata@tree <- make_tree_relative_to(strata@tree, strata@focal_species)
  strata
}
