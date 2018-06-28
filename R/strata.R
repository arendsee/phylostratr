#' Assert a Strata object is valid, die on failure
#'
#' @param required A character vector of required data fields
#' @param check_focal logical, if TRUE, then presence of the focal species in
#' the tree is checked
#' @param strata Strata object
is_valid_strata <- function(strata, required=NULL, check_focal=TRUE){
  if(class(strata) != "Strata"){
    stop("Expected Strata object, got '", class(strata), "'")
  }
  if(check_focal && !(strata@focal_species %in% strata@tree$tip.label)){
    stop("Invalid Strata object, the focal species '", strata@focal_species, "' is not found in the tree")
  }
  for(field in names(strata@data)){
    if(!all(names(strata@data[[field]]) %in% strata@tree$tip.label)){
      stop("Invalid Strata object, data field '", field, "' contains species labels that are not in the tree")
    }
  }
  if(!is.null(required)){
    missing <- setdiff(required, names(strata@data))
    if(length(missing) > 0){
      msg <- "Required data field(s) [%s] are missing from Strata object"
      msg <- sprintf(msg, paste0(missing, collapse=", "))
      stop(msg)
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
#' @param FUN The diversification algorithm to use
#' @param ... Additional arguments passed to the diversification algorithm function
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
diverse_subtree <- function(tree, n, weights=NULL, collapse=FALSE, FUN=.algo1, ...){

  if(class(tree) == 'Strata'){
    tree <- tree@tree
  }

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

  lineages <- lapply(1:nleafs(tree), lineage, type="index", x=tree)

  chosen <- FUN(lineages, n, weights, ...)

  subtree(tree, chosen, collapse=collapse)

}

.algo1 <- function(lineages, n, weights, ...){
  for(i in 1:n){
    if(i == 1){
      # start by taking the species with the highest phylogeny independent weight
      chosen <- head(which.max(weights), 1)
      # initial list of visited taxa
      seen <- lineages[[chosen]]
    } else {
      # then scale weights to penalize species that are close to the chosen species
      scaled.weights <- weights *
        sapply(lineages, function(x){
          1 - length(intersect(x, seen)) / length(x)
        })
      new_taxon <- head(which.max(scaled.weights), 1)
      seen <- union(seen, lineages[[new_taxon]])
      chosen <- c(chosen, new_taxon)
    }
  }
  chosen
}

.algo2 <- function(lineages, n, weights, ...){
  k <- rep(0, max(unlist(lineages))) 
  for(i in 1:n){

    w <- if(i == 1){
      1
    } else {
      # Calculate diversity weights: the mean number of times each ancestral
      # node has been passed through.
      w <- sapply(lineages, function(x) mean(k[x]))
      # ignore observed elements
      w[chosen] <- Inf
      w
    }

    # Divide initial weight by the diversity score
    chosen_id <- which.max(weights / w)

    # number of times each node has been passed through
    k[lineages[[chosen_id]][-1]] <- k[lineages[[chosen_id]][-1]] + 1
    chosen <- append(chosen, chosen_id)
  }
  chosen
}


#' Apply f to each outgroup branch ascending node id
#'
#' @param strata Strata object
#' @param f A function of a phylo object that return a phylo object
#' @param ... Additional arguments passed to f
#' @return Strata object
#' @export
strata_apply <- function(strata, f, ...){
  outgroups <- strata_fold(strata, f, ...)

  merge_trees <- function(trees){
    lapply(trees, tree2edgelist) %>%
      do.call(what=rbind) %>%
      unique %>%
      edgelist_to_phylo %>%
      ape::collapse.singles()
  }

  merge_data <- function(x,y){
    if(identical(x, list())){
      return(y)
    }
    stopifnot(identical(names(x), names(y)))
    z <- list()
    for(field in names(x)){
      z[[field]] <- c(x[[field]], y[[field]]) 
    }
    z
  }

  if(all(sapply(outgroups, class) == 'phylo')){
    strata@tree <- merge_trees(outgroups)
    strata@data <- lapply(strata@data, function(x){
      x[strata@tree$tip.label]
    })
    strata
  } else if(all(sapply(outgroups, class) == 'Strata')) {
    strata@tree <- lapply(outgroups, function(x) x@tree) %>% merge_trees 
    strata@data <- lapply(outgroups, function(x) x@data) %>% Reduce(f=merge_data, init=list()) 
    strata
  } else {
    stop("'f' in 'strata_apply' must map to a tree of Strata object")
  }
}

#' Apply f to each outgroup branch ascending node id, returning a list 
#'
#' @param strata Strata object
#' @param f A function of a phylo object that may return anything 
#' @param ... Additional arguments passed to f
#' @return A named list, with names corresponding to phylostrata 
#' @export
strata_fold <- function(strata, f, ...){
  is_valid_strata(strata)

  # Get list of internal IDs
  lin <- lineage(strata@tree, strata@focal_species, type='name')[-1]

  outgroups <- lapply(lin, function(ancestor){
    sisters(strata@tree, ancestor, type='index') %>%
      subtree(x=strata, collapse=FALSE, descend=TRUE, type='index') %>%
      f(...)
  })

  names(outgroups) <- tree_names(strata@tree)[parent(strata@tree, lin, type='index')]
  outgroups <- outgroups[!sapply(outgroups, is.null)]
  outgroups[[strata@focal_species]] <- f(subtree(strata, strata@focal_species, type='name'), ...)

  outgroups
}


#' Add an id to a representative list
#'
#' @param x A Strata or phylo object (all IDs MUST be NCBI taxonomy IDs)
#' @param taxa NCBI taxon ids to add 
#' @param ... Additional arguments (not used currently) 
#' @return Strata or phylo object
#' @export
#' @name add_taxa 
add_taxa <- function(x, taxa, ...){
  UseMethod('add_taxa', x)
}

#' @rdname add_taxa
#' @export
add_taxa.Strata <- function(x, taxa, ...){
  # TODO: assert that all IDs are NCBI taxonomy IDs
  x@tree <- unique(c(taxa, tree_names(x@tree))) %>%
    taxizedb::classification() %>%
    lineages_to_phylo
  x 
}

#' @rdname add_taxa
#' @export
add_taxa.phylo <- function(x, taxa, ...){
  unique(c(taxa, tree_names(x))) %>%
    taxizedb::classification() %>%
    lineages_to_phylo
}

#' Infer homology inference based on a hard e-value threshold
#'
#' @param threshold An e-value threshold
#' @param ... Unused
#' @return A function of a dataframe that has the column 'evalue', the function
#' @export
classify_by_evalue <- function(threshold, ...){
  function(x){
    !is.na(x$evalue) & (x$evalue < threshold)
  }
}

.evalue2pvalue <- function(x) { 1 - exp(-1 * x) }

#' Infer homology inference based on a hard p-value threshold
#'
#' @param threshold An e-value threshold
#' @param ... Unused
#' @return A function of a dataframe that has the column 'evalue', the function
#' @export
classify_by_pvalue <- function(threshold, ...){
  function(x){
    p <- .evalue2pvalue(x$evalue)
    !is.na(p) & (p < threshold)
  }
}

#' Classify homologs under the independence assumption
#'
#' @param threshold The target alpha (e.g. 0.05)
#' @param ... Additional arguments passed to p.adjust
#' @export
classify_assuming_iid <- function(threshold, ...){
  function(x){
    exp_n <- length(unique(x$qseqid))

    m <- dplyr::group_by(x, qseqid, mrca) %>%
    dplyr::mutate(
      pval.adj = p.adjust(.evalue2pvalue(evalue), ...)
    ) %>%
    dplyr::select(qseqid, pval.adj, mrca)

    z <- dplyr::group_by(m, .data$qseqid, .data$mrca) %>%
      dplyr::summarize(pval = min(pval.adj)) %>%
      # cast and melt: this completes the data 
      reshape2::acast(qseqid ~ mrca, value.var="pval", fill=1) %>%
      apply(1, p.adjust, ...) %>% t

    # Add in any rows that were missing, this can happedn when, for reasons
    # most mysterious, the focal gene does not even match itself. 

    p <- z[as.matrix(x[, c('qseqid', 'mrca')])]

    if(exp_n != nrow(z)){
      stop("classify_assuming_iid is broken - complain to the maintainer")
    }

    { !is.na(p) & p < threshold }
  }
}

#' Infer homology based on p-value with Fisher's method 
#'
#' @param threshold P-value threshold
#' @export
classify_by_fisher <- function(threshold){
  fishers_method <- function(p){
      q <- -2 * sum(log(p))
      df <- 2 * length(p) # degrees of freedom
      pchisq(q, df, lower.tail=FALSE)
  }
  function(x){
    x$pval <- ifelse(is.na(x), 1, .evalue2pvalue(x$evalue))
    dplyr::group_by(x, .data$qseqid, .data$mrca) %>%
      dplyr::mutate(pval = fishers_method(.data$pval)) %>%
      { .$pval < threshold }
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
    # convert to name if possible
    partial_id_to_name %>%
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
  orphan_ps <- max(hittable$ps)
  orphan_mrca <- levels(strata_names)[length(strata_names)]
  hittable[classifier(hittable), ] %>%
    dplyr::select(.data$qseqid, .data$mrca, .data$ps) %>%
    dplyr::group_by(.data$qseqid) %>%
    dplyr::filter(.data$ps == min(.data$ps)) %>%
    dplyr::ungroup() %>%
    {
      strata <- .
      rbind(
        strata,
        # add in queries that have no hit against anything
        hittable[!(hittable$qseqid %in% strata$qseqid), ] %>%
        dplyr::select(.data$qseqid, .data$mrca, .data$ps) %>%
        # take just one entry for each query
        dplyr::group_by(.data$qseqid) %>%
        head(1) %>%
        dplyr::ungroup() %>%
        # assign it orphan status
        dplyr::mutate(ps = orphan_ps, mrca = orphan_mrca)
      )
    } %>%
    dplyr::distinct() %>%
    dplyr::mutate(mrca_name = strata_names[.data$ps])
}

#' Given a strata table, remove given strata 
#'
#' The genes in a stratum that is removed will be reassigned to the parent
#' stratum. If the stratum is the root stratum, the genes will instead be
#' reassigned to the child stratum.
#' 
#' Note that removing strata will result in non-monophyletic branches.
#'
#' @param strata_table data.frame with the factored column 'mrca_name'
#' @param strata_names character vector of strata to drop
#' @export
prune_phylostrata <- function(strata_table, strata_names){
  if(!('mrca_name' %in% names(strata_table))){
    stop("strata_table must have the column mrca_name")
  }
  if(!is.factor(strata_table$mrca_name)){
    stop("strata_table$mrca_name must be a factor ordered by descending age")
  }

  mrca_levels <- levels(strata_table$mrca_name)
  id2name_map <- dplyr::distinct(strata_table[, c('mrca', 'mrca_name')])

  indices <- as.integer(strata_table$mrca_name)

  indices_to_drop <- intersect(strata_names, mrca_levels) %>%
    {which(mrca_levels %in% .)} %>%
    sort(decreasing=TRUE)

  for(i in indices_to_drop){
    indices <- ifelse(indices >= i, indices-1, indices)
    mrca_levels <- mrca_levels[-i]
  }

  # TODO: I am making the mistake of storing the same information several
  # times: ps, mrca, and mrca_name. Then everytime I change one I have to
  # change them all. I should add a phylostratigrph object that handles all
  # this for me.
  strata_table$mrca_name <- ifelse(indices > 0, mrca_levels[indices], mrca_levels[1]) 
  strata_table$mrca <- NULL
  strata_table <- merge(strata_table, id2name_map)
  strata_table$mrca_name <- factor(strata_table$mrca_name, levels=mrca_levels)
  strata_table$ps <- as.integer(strata_table$mrca_name)
  strata_table

}

#' Standardize two or more phylostrata tables 
#'
#' Each table is subset with the union of sequence IDs. Also all strata that
#' are not shared between all input tables are dropped, with the genes
#' reassigned to older strata as needed. 
#'
#' @param strata_tables phylostrata tables, all must the column 'mrca_name'
#' @return list of standardized phylostrata tables 
#' @export
standardize_strata <- function(strata_tables){

  # Get intersection of all sequence IDs
  common_ids <- Reduce(
    f = function(x,y){intersect(x$qseqid, y$qseqid)},
    x = strata_tables[-1],
    init = strata_tables[[1]]
  )

  # Get the phylostrata common to all studies
  common_strata <- Reduce(
    f = function(x,y){intersect(levels(x$mrca_name), levels(y$mrca_name))},
    x = strata_tables[-1],
    init = strata_tables[[1]]
  )

  tuplify(strata_tables) %>%
    lapply(function(x){
      d <- x$value
      d$group <- x$name
      d <- subset(d, d$qseqid %in% common_ids)
      d <- prune_phylostrata(d, setdiff(levels(d$mrca_name), common_strata))
      x$value <- d
      x
    }) %>%
    untuplify
}

#' Convert strata data, tip, and node names to and from NCBI taxonomy IDs
#'
#' @param strata Strata object
#' @param target What to convert, may be 'tip', 'node', or 'all'
#' @param to What to convert to, may be 'id' or 'name'
#' @param FUN a custom function to convert names
#' @return Strata object with new names
#' @export
strata_convert <- function(strata, target='tip', to='id', FUN=NULL){
  if(is.null(FUN)){
    FUN <- switch(
      to,
      id = taxizedb::name2taxid,
      name = taxizedb::taxid2name,
      stop("Illegal 'to' value: must be either 'id' or 'name'")
    )
  }
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

#' Convert vector of mixed IDs and names to a vector of names
#'
#' @param x character vector
#' @return Vector of names
#' @export
partial_id_to_name <- function(x){
  is_id <- grepl('^[0-9]+$', x, perl=TRUE)
  x[is_id] <- taxizedb::taxid2name(x[is_id])
  x
}

#' Get a phylostrata map from a Strata object
#'
#' @param strata Strata object
#' @return data.frame with columns 'species', 'mrca', an 'ps', where 'mrca' is
#' a factor of phylostrata with levels ordered by age
#' @export
get_phylostrata_map <- function(strata){
  map <- strata_fold(strata, leafs) %>%
    tuplify %>%
    lapply(function(x){
      tibble::data_frame(
        species = x$value,
        mrca = rep(x$name, length(x$value)),
        ps = x$position 
      )
    }) %>%
    do.call(what=rbind)

  mrca_levels <- map %>%
    dplyr::select(.data$mrca, .data$ps) %>%
    dplyr::distinct() %>%
    dplyr::arrange(.data$ps) %>% {.$mrca}

  map$mrca <- factor(map$mrca, levels=mrca_levels)

  map
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

#' Rename a species in a Strata object
#'
#' @param strata Strata object
#' @param old_name character vector of names that are currently in the Strata object
#' @param new_name character vector of names to replace the old_name vector
#' @return Strata object with the same number of species
#' @examples
#' \dontrun{
#' data(saccharomyces)
#' x <- rename_species(saccharomyces, "Saccharomyces_arboricola", "Sa")
#' is_valid_strata(x)
#'
#' x <- rename_species(saccharomyces,
#'   c("Saccharomyces_arboricola", "Saccharomyces_cerevisiae"), c("Sa", "Sc"))
#' is_valid_strata(x)
#' }
rename_species <- function(strata, old_name, new_name){
  is_valid_strata(strata)
  if(length(old_name) != length(new_name)){
    stop("old_name must have the same length as the new_name")
  }
  if(any(!(old_name %in% strata@tree$tip.label))){
    name <- which(!(old_name %in% strata@tree$tip.label))
    stop(old_name[name]," is not in the tree")
  }
  for (i in 1:length(old_name)){
    old <- which(strata@tree$tip.label %in% old_name[i])
    strata@tree$tip.label[old] <- new_name[i]
    for (j in names(strata@data)){
      names(strata@data[[j]])[old] <- new_name[i]
    }
  }
  strata
}
