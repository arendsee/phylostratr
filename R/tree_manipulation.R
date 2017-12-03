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
    tree$node.label <- default_names
  } else {
    tree$node.label <- ifelse(
      is.na(tree$node.label),
      default_names,
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
#' @param type ['auto', 'index', or 'name'], is 'index' then the id is expected
#' to be an index in the phylo object, if 'name' it is expected to be a node or
#' tip name, if 'auto' the value is assumed to be a name if it is a character
#' string.
#' @return on success, return a valid integrel id
clean_phyid <- function(tree, id, len=NULL, type='auto'){
  is_integer <- function(i){
    if(is.character(i)){
      all(grepl("^[0-9]+$", id, perl=TRUE))
    } else if(is.numeric(i)){
      !all((i %% 1) > 1e-6)
    } else {
      FALSE
    }
  }
  if(!is.null(len) && length(id) != len) {
    stop("Expected ", len, " id(s), got ", length(id))
  }
  if(length(id) == 0){
    return(integer(0)) 
  }
  index <- if(type == 'auto'){
    # If all ids match node/tip names, assume they are names ...
    if(all(id %in% tree_names(tree))){
      match(id, tree_names(tree))
    # if all the ids are integers, assume they are indices
    } else if(is_integer(id)) {
      as.integer(id)
    # otherwise, die screaming
    } else {
      stop("Could not interpret ids")
    }
  } else if(type == 'index'){
    if(is_integer(id)) {
      as.integer(id)
    } else {
      stop("Could not find requested ids")
    }
  } else if(type == 'name'){
      if(all(id %in% tree_names(tree))){
        match(id, tree_names(tree))
      } else {
        stop("Could not find nodes with requested names")
      }
  } else {
    stop('Unsupported type')
  }
  if(max(index, na.rm=TRUE) > tree_size(tree)){
    stop("Cannot get node ", max(id, na.rm=TRUE), " for tree of size ", tree_size(tree))
  }
  if(min(index, na.rm=TRUE) < 1){
    stop("Invalid id, node ids must be greater than 0")
  }
  index
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
#' @param byname logical Return leaf labels, rather than indices
#' @export
leafs <- function(tree, byname=FALSE){
  ids <- 1:length(tree$tip.label)
  if(byname){
    tree$tip.label[ids]
  } else {
    ids
  }
}

#' Get node indices
#'
#' @param tree phylo object
#' @param rebase logical Make indices base 1
#' @param byname logical Return node labels, rather than indices
#' @export
nodes <- function(tree, rebase=FALSE, byname=FALSE){
  ids <- (length(tree$tip.label)+1):(nleafs(tree) + tree$Nnode)
  ids_base1 <- ids - nleafs(tree)
  if(byname){
    if(is.null(tree$node.label))
      stop("Cannot get nodes by name since tree has no node labels")
    tree$node.label[ids_base1]
  } else if(rebase){
    ids_base1
  } else {
    ids
  }
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
#' @param tree phylo object
#' @param id A single name or index
#' @param use_name Should names be used rather than internal ids
#' @examples
#' data(atree)
#' lineage(atree, 't1')
#' lineage(atree, 1)
#' @export
lineage <- function(tree, id, use_name=FALSE){
  id <- clean_phyid(tree, id, len=1)
  id <- if(is_root(tree, id)){
    id
  } else {
    c(lineage(tree, parent(tree, id)), id)
  }
  if(use_name){
    tree_names(tree)[id]
  } else {
    id
  }
}

#' Get the number of species and ancestors in the tree
#'
#' @param tree phylo object
#' @return integer
#' @export
tree_size <- function(tree){
  max(nodes(tree))
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
  sort(tree$edge[tree$edge[,1] == id, 2])
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
  ids <-
    if(length(children) > 0)
      c(id, descendent_nodes(tree, children))
    else
      id
  sort(ids)
}

#' Get all leafs descending from a set of ids
#'
#' This is a simple wrapper around \code{descendent_nodes}.
#'
#' @param tree phylo object
#' @param id vector of ids or names
#' @return integer vector of leaf ids
#' @export
descendents <- function(tree, id){
  descendent_nodes(tree, id) %>%
    intersect(leafs(tree)) %>%
    sort
}

#' Select a subset of nodes from a tree
#'
#' @param tree phylo object
#' @param id vector of ids or names
#' @param collapse logical collapse edges around nodes with a single descendent
#' @param descend logical include all descendents of each id
#' @param type id type: ['name', 'id', 'auto']
#' @return phylo object
#' @examples
#' data(atree)
#' subtree(atree, c('n15', 'n19'))
#' subtree(atree, c('t7', 't4', 't1'))
#' @export
subtree <- function(tree, id, collapse=TRUE, descend=TRUE, type='auto'){
  id <- clean_phyid(tree, id, type=type)
  if(descend){
    id <- descendent_nodes(tree, id)
  }
  leaf_ids <- intersect(leafs(tree), id)
  id <- lapply(leaf_ids, lineage, tree=tree) %>%
    do.call(what=c) %>% unique
  new_edge <- tree$edge[(tree$edge[, 1] %in% id) & (tree$edge[, 2] %in% id), , drop=FALSE]
  new_Nnode <- length(unique(new_edge[, 1])) 
  new_size <- length(unique(c(new_edge[, 1], new_edge[, 2])))
  leaf_ids <- intersect(leafs(tree), id)
  node_ids <- intersect(nodes(tree), id)
  idmap <- seq_along(c(leaf_ids, node_ids))
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
  if(nleafs(new_tree) == 1 && collapse){
    # ape::collapse.singles dies one single tip input. So I need to implement
    # my own handling for this edge case. I keep the parent of the single node
    # to keep the ape from screaming later. 
    pid <- new_tree$edge[new_edge[, 2] == 1, 1]
    new_tree$node.label <- new_tree$node.label[pid - 1]
    new_tree$edge <- new_edge[new_edge[, 2] == 1, , drop=FALSE]
    new_tree$Nnode <- 1
  } else if(new_size > 1 && collapse){
    # collapse the edges of nodes that have only a single descendent 
    new_tree <- ape::collapse.singles(new_tree)
  }
  new_tree
}

#' Get the sisters of a node
#'
#' @param tree phylo object
#' @param id vector of ids or names
#' @return integer vector of indices
#' @export
sisters <- function(tree, id){
  id <- clean_phyid(tree, id, len=1)
  if(is_root(tree, id)){
    integer(0)
  } else {
    setdiff(children(tree, parent(tree, id)), id)
  }
}

#' Get list of sister trees for a given node
#'
#' @param tree phylo object
#' @param id vector of ids or names
#' @return list of phylo objects
#' @export
sister_trees <- function(tree, id){
  sister_ids <- sisters(tree, id)
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
#' @param type id type: ['name', 'id', 'auto']
#' @export
prune <- function(tree, id, type='auto'){
  id <- clean_phyid(tree, id, type=type)
  id <- descendents(tree, id)
  subtree(tree, setdiff(leafs(tree), id))
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
    subtree(tree=tree)
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
