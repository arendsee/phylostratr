#' Convert a tree to an edge list
#'
#' @param tree phylo object
#' @export
tree2edgelist <- function(tree){
  if(is.null(tree$node.label)){
    stop("Cannot convert this tree to an edge list since it lacks node labels")
  }
  edge <- tree$edge
  matrix(c(
      tree_names(tree)[edge[ , 1]],
      tree_names(tree)[edge[ , 2]]
    ),
    ncol=2,
    byrow=FALSE
  )
}

#' Make a tree from the union of lineages
#'
#' @param lineages list of lineages, where each lineage is a data.frame with
#' the columns id and name
#' @param clean Should leafs with descendents be removed? This can occur when
#' both a species and a descendent subspecies are in the lineage set.
#' @return phylo object with node names and no branch lengths
#' @export
lineages_to_phylo <- function(lineages, clean=FALSE){
  to_edge <- function(xs){    # LTC comment: this makes a data table with all the unique connections between e.g., genus (col 1) and sp (col 2), order (col 1) to family (col 2) etc
    # build edge list based on taxonomy IDs
    lapply(xs, function(x){
      matrix(c(head(x$id, -1), tail(x$id, -1)), ncol=2)
    }) %>%
    do.call(what=rbind) %>%
    unique
  }
  edges <- to_edge(lineages)
  if(clean){
    bad <- names(lineages)[names(lineages) %in% edges[, 1]]
    edges <- Filter(x=lineages, f=function(x){
      !any(head(x$id, -1) %in% bad)      
    }) %>% to_edge
  }
  edgelist_to_phylo(edges) %>% ape::collapse.singles()
}

#' Make a tree of ancestors from a lineage
#'
#' @param lineage A list of ancestors, ranked from old to young
#' @return phylo object with 0 length branches
#' @export
lineage_to_ancestor_tree <- function(lineage){
  edgelist_to_phylo(matrix(c(head(lineage, -1), tail(lineage, -1)), ncol=2)) 
}

#' Make a phylo object from an edgelist
#'
#' @param edgelist A matrix of node names
#' @return phylo object
#' @export
#' @examples
#' edgelist <- matrix(c('A', 'A', 'B', 'B', 'B', 'C', 'D', 'E', 'A', 'F'), ncol=2)
#' tree <- edgelist_to_phylo(edgelist)
#' # A bug in ape causes a segfault here
#' # plot(tree, show.node.label=TRUE)
edgelist_to_phylo <- function(edgelist){
  from <- edgelist[,1]
  to <- edgelist[,2]
  ids <- unique(c(edgelist[,1], edgelist[,2]))

  tip.label <- setdiff(ids, from)
  node.label <- unique(from)

  # make a map from taxonomy ID to internal 1:n ids
  idmap <- 1:(length(tip.label) + length(node.label))
  names(idmap) <- c(tip.label, node.label)

  # make a phylo object
  tree <- list(
    edges      = matrix(c(idmap[as.character(from)], idmap[as.character(to)]), ncol=2),
    tip.label  = unname(tip.label),
    node.label = unname(node.label),
    Nnode      = length(node.label)
  )
  class(tree) <- 'phylo'

  tree
}
