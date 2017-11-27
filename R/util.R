#' Get the maximum hit 
#'
#' The resulting data.frame will be complete, with every (qseqid, staxid) pair represented
#'
#' @param d A dataframe with the columns 'qseqid', 'staxid', and 'score'
#' @export
get_max_hit <- function(d){
  d %>%
    dplyr::group_by(.data$qseqid, .data$staxid) %>%
    dplyr::filter(.data$score == max(.data$score)) %>%
    # the filter step can lead to multiple hits with equal score, I want just
    # one, so have to pass through distinct
    dplyr::ungroup() %>%
    dplyr::distinct(.data$qseqid, .data$staxid, .keep_all=TRUE)
}

maybe_message <- function(msg, verbose=TRUE, ...){
  if(verbose)
    message(sprintf(msg, ...))
}

#' Make a tree from the union of lineages
#'
#' @param lineages list of lineages, where each lineage is a data.frame with
#' the columns id and name
#' @return phylo object with node names and no branch lengths
lineages_to_phylo <- function(lineages){
  # get species names
  tip.label <- unique(sapply(lineages, function(x) tail(x$name, 1)))

  # get species taxon IDs
  tip.id <- unique(sapply(lineages, function(x) tail(x$id, 1)))

  # build edge list based on taxonomy IDs
  edges <- unique(do.call(what=rbind,
    lapply(lineages, function(x){
      matrix(c(head(x$id, -1), tail(x$id, -1)), ncol=2)
    })
  ))

  # get node IDs
  node.id <- setdiff(c(edges[,1], edges[,2]), tip.id)

  # make a map from taxonomy ID to internal 1:n ids
  idmap <- 1:(length(tip.label) + length(node.id))
  names(idmap) <- c(tip.id, node.id)

  # make a phylo object
  tree <- list(
    edges      = matrix(c(idmap[edges[,1]], idmap[edges[,2]]), ncol=2),
    tip.label  = unname(tip.label),
    node.label = unname(node.id),
    Nnode      = length(node.id)
  )
  class(tree) <- 'phylo'

  tree
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
#' # A bug in ape causes a segfault here, but the French bastards in charge don't
#' # seem to take the issues seriously
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
    edges      = matrix(c(idmap[from], idmap[to]), ncol=2),
    tip.label  = unname(tip.label),
    node.label = unname(node.label),
    Nnode      = length(node.label)
  )
  class(tree) <- 'phylo'

  tree
}

#' Convert node/tip name to id, check size
#'
#' @param tree phylo object
#' @param id integer id or a node label
#' @param len require id to be of this length (NULL for no requirement)
#' @return on success, return a valid integrel id
vet_phyid <- function(tree, id, len=NULL){
  if(length(id) == 0){
    return(integer(0)) 
  }
  if(is.character(id)){
    names <- tree_names(tree)
    if(all(id %in% names)){
      id <- match(id, names)
    } else if(!all(grepl("^[0-9]+$", id, perl=TRUE))) {
      stop("Could not find nodes with the requested labels")
    }
  }
  if(!is.null(len) && length(id) != len)
    stop("Expected ", len, " id(s), got ", length(id))
  id <- as.integer(id)
  if(max(id, na.rm=TRUE) > tree_size(tree))
    stop("Cannot get node ", max(id, na.rm=TRUE), " for tree of size ", tree_size(tree))
  if(min(id, na.rm=TRUE) < 1)
    stop("Invalid id, node ids must be greater than 0")
  id
}

tree_names <- function(tree){
  c(tree$tip.label, tree$node.label)
}

ntip <- function(tree){
  length(tree$tip.label)
}

#' Find the ancestors of a node
#'
#' @examples
#' set.seed(42)
#' a <- rtree(10)
#' a$tip.label <- paste0("t", 1:ntip(a))
#' a$node.label <- paste0("n", 1:a$Nnode + ntip(a))
lineage <- function(tree, id){
  id <- vet_phyid(tree, id, len=1)
  if(is_root(tree, id)){
    id
  } else {
    c(lineage(tree, parent(tree, id)), id)
  }
}

#' Get the number of species and ancestors in the tree
#'
#' @param tree phylo object
#' @return integer
tree_size <- function(tree){
  ntip(tree) + tree$Nnode
}

#' Vectorized node parent getter
#'
#' @param tree phylo object
#' @param id vector of ids or names
#' @return vector of parent ids (NA if root)
parent <- function(tree, id=1:tree_size(tree)){
  id <- vet_phyid(tree, id)
  tree$edge[match(id, tree$edge[, 2]), 1]
}

#' Get the immediate children of a node
#'
#' @param tree phylo object
#' @param id vector of ids or names
#' @return vector of children, (integer(0) if there are no children)
children <- function(tree, id){
  id <- vet_phyid(tree, id, len=1)
  tree$edge[tree$edge[,1] == id, 2]
}

#' Vectorized tree root finder
#'
#' @param tree phylo object
#' @param id vector of ids or names
#' @return logical vector
is_root <- function(tree, id=1:tree_size(tree)){
  is.na(parent(tree, id))
}

#' Get the root of a tree
#'
#' There should normally be only one, but if there are more, this function will
#' find all of them.
#'
#' @param tree phylo object
get_root <- function(tree){
  which(is_root(tree))
}

#' Get all nodes (leafs and ancestors) descending from, and including, a set of ids
#'
#' @param tree phylo object
#' @param id vector of ids or names
#' @return integer vector of leaf and node ids
descendent_nodes <- function(tree, id){
  id <- vet_phyid(tree, id)
  children <- tree$edge[tree$edge[,1] %in% id, 2]
  if(length(children) > 0)
    c(id, descendent_nodes(tree, children))
  else
    id
}

#' Get all leafs descending from set of ids
#'
#' This is a simple wrapper around \code{descendent_nodes}.
#'
#' @param tree phylo object
#' @param id vector of ids or names
#' @return integer vector of leaf ids
descendents <- function(tree, id){
  descendent_nodes(tree, id) %>%
    { .[. <= length(tree$tip.label)] }
}

subset_phylo <- function(tree, id){
  id <- vet_phyid(tree, id)
  new_edge <- tree$edge[tree$edge[, 1] %in% id, ]
  new_Nnode <- length(unique(new_edge[, 1])) 
  new_size <- length(unique(c(new_edge[, 1], new_edge[, 2])))
  leaf_ids <- intersect(1:ntip(tree), id)
  node_ids <- intersect((ntip(tree)+1):(tree_size(tree)), id)
  idmap <- 1:new_size
  names(idmap) <- c(leaf_ids, node_ids)
  new_edge[, 1] <- idmap[as.character(new_edge[, 1])]
  new_edge[, 2] <- idmap[as.character(new_edge[, 2])]
  new_tip.label <- tree$tip.label[leaf_ids]
  new_node.label <- tree$node.label[node_ids - ntip(tree)]
  new_tree <- list(
    edge        = new_edge,
    tip.label   = new_tip.label,
    edge.length = tree$edge.length[tree$edge[, 1] %in% id],
    Nnode       = new_Nnode,
    node.label  = new_node.label
  )
  class(new_tree) <- 'phylo'
  if(sum(is_root(new_tree)) != 1){
    stop("Malformed tree: must have exactly 1 root")
  }
  new_tree
}

#' Get a subtree rooted at a given node
#'
#' @param tree phylo object
#' @param id vector of ids or names
#' @return phylo object
subtree <- function(tree, id){
  id <- vet_phyid(tree, id, len=1)
  nodes <- descendent_nodes(tree, id)
  subset_phylo(tree, nodes)
}

#' Get list of sister trees for a given node
#'
#' @param tree phylo object
#' @param id vector of ids or names
#' @return list of phylo objects
sister_trees <- function(tree, id){
  id <- vet_phyid(tree, id, len=1)
  sister_ids <- setdiff(children(tree, parent(tree, id)), id)
  sister_trees <- lapply(sister_ids, subtree, tree=tree)
  if(length(sister_trees) > 0){
    names(sister_trees) <- sister_ids
  }
  sister_trees
}

# prune <- function(tree, id){
#   id <- vet_phyid(tree, id)
# }
#
# strata_apply <- function(strata, f){
#
# }
#
# strata_join <- function(trees){
#
# }
