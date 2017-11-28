#' If nodes are not named, give them default names
#'
#' In phylostratr, all nodes in a tree need to be named, since the names are
#' used to merge trees.
#'
#' @param tree phylo object
#' @param default_names node names to use if they are missing. Must be a vector
#' with a number of elements equal to the number of nodes in the tree.
#' @return phylo object
#' @export
set_node_names <- function(tree, default_names=paste0("n", nodes(tree))){
  if(is.null(tree$node.label)){
    tree$node.label <- default_node_names
  } else {
    tree$node.label <- ifelse(
      is.na(tree$node.label),
      default_node_names,
      tree$node.label
    )
  }
  tree
}

#' Convert node/tip name to id, check size
#'
#' @param tree phylo object
#' @param id integer id or a node label
#' @param len require id to be of this length (NULL for no requirement)
#' @return on success, return a valid integrel id
clean_phyid <- function(tree, id, len=NULL){
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

#' Get the names of all leafs and nodes in a tree
#'
#' @param tree phylo object
#' @export
tree_names <- function(tree){
  c(tree$tip.label, tree$node.label)
}

#' Get leaf indices
#'
#' @param tree phylo object
#' @export
leafs <- function(tree){
  1:length(tree$tip.label)
}

#' Get node indices
#'
#' @param tree phylo object
#' @export
nodes <- function(tree){
  (length(tree$tip.label)+1):tree_size(tree)
}

#' Get the number of leafs
#'
#' @param tree phylo object
#' @export
nleafs <- function(tree){
  length(leafs(tree))
}

#' Find the ancestors of a node
#'
#' @examples
#' data(atree)
#' lineage(atree, 't1')
#' lineage(atree, 1)
#' @export
lineage <- function(tree, id){
  id <- clean_phyid(tree, id, len=1)
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
#' @export
tree_size <- function(tree){
  max(node_ids(tree))
}

#' Vectorized node parent getter
#'
#' @param tree phylo object
#' @param id vector of ids or names
#' @return vector of parent ids (NA if root)
#' @export
parent <- function(tree, id=1:tree_size(tree)){
  id <- clean_phyid(tree, id)
  tree$edge[match(id, tree$edge[, 2]), 1]
}

#' Get the immediate children of a node
#'
#' @param tree phylo object
#' @param id vector of ids or names
#' @return vector of children, (integer(0) if there are no children)
#' @export
children <- function(tree, id){
  id <- clean_phyid(tree, id, len=1)
  tree$edge[tree$edge[,1] == id, 2]
}

#' Vectorized tree root finder
#'
#' @param tree phylo object
#' @param id vector of ids or names
#' @return logical vector
#' @export
is_root <- function(tree, id=1:tree_size(tree)){
  is.na(parent(tree, id))
}

#' Get the root of a tree
#'
#' There should normally be only one, but if there are more, this function will
#' find all of them.
#'
#' @param tree phylo object
#' @return index of the root node (or nodes if there are more than one)
#' @export
get_root <- function(tree){
  which(is_root(tree))
}

#' Get all nodes (leafs and ancestors) descending from, and including, a set of ids
#'
#' @param tree phylo object
#' @param id vector of ids or names
#' @return integer vector of leaf and node ids
#' @export
descendent_nodes <- function(tree, id){
  id <- clean_phyid(tree, id)
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
#' @export
descendents <- function(tree, id){
  descendent_nodes(tree, id) %>%
    intersect(leafs(tree))
}

#' Select a subset of nodes from a tree
#'
#' @param tree phylo object
#' @param id vector of ids or names
#' @param collapse logical collapse edges around nodes with a single descendent
#' @param descend logical include all descendents of each id
#' @return phylo object
#' @examples
#' data(atree)
#' subset_phylo(atree, c('n15', 'n19'))
#' subset_phylo(atree, c('t7', 't4', 't1'))
#' @export
subset_phylo <- function(tree, id, collapse=TRUE, descend=TRUE){
  id <- clean_phyid(tree, id)
  if(descend){
    id <- descendent_nodes(tree, id)
  }
  leaf_ids <- intersect(leafs(tree), id)
  id <- lapply(leaf_ids, lineage, tree=tree) %>%
    do.call(what=c) %>% unique
  new_edge <- tree$edge[(tree$edge[, 1] %in% id) & (tree$edge[, 2] %in% id), ]
  new_Nnode <- length(unique(new_edge[, 1])) 
  new_size <- length(unique(c(new_edge[, 1], new_edge[, 2])))
  leaf_ids <- intersect(leafs(tree), id)
  node_ids <- intersect(nodes(tree), id)
  idmap <- 1:new_size
  names(idmap) <- c(leaf_ids, node_ids)
  new_edge[, 1] <- idmap[as.character(new_edge[, 1])]
  new_edge[, 2] <- idmap[as.character(new_edge[, 2])]
  new_tip.label <- tree$tip.label[leaf_ids]
  new_node.label <- tree$node.label[node_ids - nleafs(tree)]
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
  if(collapse){
    # collapse the edges of nodes that have only a single descendent 
    new_tree <- collapse.singles(new_tree)
  }
  new_tree
}

#' Get a subtree rooted at a given node
#'
#' @param tree phylo object
#' @param id vector of ids or names
#' @return phylo object
#' @export
subtree <- function(tree, id){
  id <- clean_phyid(tree, id, len=1)
  subset_phylo(tree, descendent_nodes(tree, id))
}

#' Get list of sister trees for a given node
#'
#' @param tree phylo object
#' @param id vector of ids or names
#' @return list of phylo objects
#' @export
sister_trees <- function(tree, id){
  id <- clean_phyid(tree, id, len=1)
  sister_ids <- setdiff(children(tree, parent(tree, id)), id)
  sister_trees <- lapply(sister_ids, subtree, tree=tree)
  if(length(sister_trees) > 0){
    names(sister_trees) <- sister_ids
  }
  sister_trees
}

#' Remove the specified indices, and their descendents, from the tree
#'
#' @param tree phylo object
#' @param id vector of indices or names 
#' @export
prune <- function(tree, id){
  id <- clean_phyid(tree, id)
  id <- descendents(tree, id)
  subset_phylo(tree, setdiff(leafs(tree), id))
}

#' Merge fully named subtrees according to a reference tree
#'
#' @param tree phylo object
#' @param subtrees A list of phylo objects. All nodes and leafs of each tree
#' must be named. All names must be from 'tree'.
#' @return phylo object
#' @examples
#' data(atree)
#' subtrees <- list(
#'   subtree(atree, 'n19'),
#'   subtree(atree, 'n18'),
#'   subtree(atree, 'n15')
#' )
#' merge_phylo(atree, subtrees)
#' @export
merge_phylo <- function(tree, subtrees){
  lapply(subtrees, tree_names) %>% do.call(what='c') %>% unique %>%
    subset_phylo(tree=tree)
}



#' Map indices from a to b based off tip labels
#'
#' @param a phylo object
#' @param b phylo object
#' @return matrix mapping indices from `a` to `b`
map_ids <- function(a, b){
  # TODO: assert a and b have the same topology
  # start with map from the tips of `a` to those of `b`
  idmaps <- list(matrix(
      c(
        clean_phyid(a, a$tip.label),
        clean_phyid(b, a$tip.label)
      ), ncol=2))
  # iteratively find the parents of nodes, until root is reached
  while(TRUE){
    last_a <- tail(idmaps, n=1)[[1]][, 1]
    last_b <- tail(idmaps, n=1)[[1]][, 2]
    if(!all(is_root(a, last_a))){
      idmaps <- append(
        idmaps,
        list(unique(matrix(c(
            parent(a, last_a),
            parent(b, last_b)
          ), ncol=2)))
      )
    } else {
      break
    }
  }
  # merge the tables
  idmap <- do.call(rbind, idmaps)
  # get the unique rows that are not at root
  idmap <- unique(idmap[!is.na(idmap[, 1]), ])
  idmap
}

#' Map nodes names from a to b
#'
#' The trees a and by must have the same topology and same tip labels.
#'
#' @param a phylo object with node labels
#' @param b phylo object
map_node_label <- function(a, b){
  # get map of indices from a to b
  idmap <- map_ids(a, b)
  # get the number of tips
  offset <- length(a$tip.label)
  # get the node indices, these are above the tip indices
  idmap <- idmap[idmap[, 1] > offset, ]
  # rebase
  idmap[, 1] <- idmap[, 1] - offset
  idmap[, 2] <- idmap[, 2] - offset
  # map names from a to b
  b$node.label[idmap[, 2]] <- a$node.label[idmap[, 1]]
  b
}