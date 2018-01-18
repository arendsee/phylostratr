# Strata / phylo compatibility function
.fmap <- function(x, f, fout=.identity, ...){
  fout(x, f(x@tree, ...))
}
.strata <- function(x, tree){
  new_data <- lapply(x@data, function(datum){
    datum[tree$tip.label]
  })
  x@tree <- tree 
  x@data <- new_data
  x
}
.name <- function(x, indices){
  tree_names(x)[indices] 
}
.identity <- function(x, result) result

#' If nodes are not named, give them default names
#'
#' In phylostratr, all nodes in a tree need to be named, since the names are
#' used to merge trees.
#'
#' @param x phylo or Strata object
#' @param default_names node names to use if they are missing. Must be a vector
#' with a number of elements equal to the number of nodes in the tree.
#' @param ... Arguments passed to \code{set_node_names.phylo}
#' @return phylo or Strata object
#' @export
#' @name set_node_names
set_node_names <- function(x, ...){
  UseMethod('set_node_names', x)
}

#' @rdname set_node_names
#' @export
set_node_names.Strata <- function(x, ...){
  .fmap(x, set_node_names, fout=.strata, ...)
}

#' @rdname set_node_names
#' @export
set_node_names.phylo <- function(x, default_names=paste0("n", nodes(x)), ...){
  if(is.null(x$node.label)){
    x$node.label <- default_names
  } else {
    x$node.label <- ifelse(
      is.na(x$node.label),
      default_names,
      x$node.label
    )
  }
  x
}

#' Convert node/tip name to id, check size
#'
#' @param x phylo or Strata object
#' @param id integer id or a node label
#' @param len require id to be of this length (NULL for no requirement)
#' @param type ['auto', 'index', or 'name'], is 'index' then the id is expected
#' to be an index in the phylo object, if 'name' it is expected to be a node or
#' tip name, if 'auto' the value is assumed to be a name if it is a character
#' string.
#' @return on success, return a valid integrel id
clean_phyid <- function(x, id, len=NULL, type='auto'){
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
    id_name <- if(all(id %in% tree_names(x))){
      match(id, tree_names(x))
    }
    # if all the ids are integers, assume they are indices
    id_index <- if(is_integer(id)) {
      as.integer(id)
    }
    if(!(is.null(id_name) || is.null(id_index))){
      stop("Cannot automatically resolve tree id.",
           "It could be either a phylo object index or a taxon name.",
           "Please specify the type.")
    } else if(!is.null(id_name)){
      id_name
    } else if(!is.null(id_index)){
      id_index
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
      if(all(id %in% tree_names(x))){
        match(id, tree_names(x))
      } else {
        stop("Could not find nodes with requested names")
      }
  } else {
    stop('Unsupported type')
  }
  if(max(index, na.rm=TRUE) > tree_size(x)){
    stop("Cannot get node ", max(id, na.rm=TRUE), " for tree of size ", tree_size(x))
  }
  if(min(index, na.rm=TRUE) < 1){
    stop("Invalid id, node ids must be greater than 0")
  }
  index
}


#' Get the names of all leafs and nodes in a tree
#'
#' @param x phylo object
#' @export
#' @name tree_names
tree_names <- function(x){
  UseMethod('tree_names', x)
}

#' @rdname tree_names
#' @export
tree_names.Strata <- function(x){
  .fmap(x, tree_names, fout=.identity)
}

#' @rdname tree_names
#' @export
tree_names.phylo <- function(x){
  c(x$tip.label, x$node.label)
}


#' Get leaf indices
#'
#' @param x phylo object
#' @param byname logical Return leaf labels, rather than indices
#' @param ... Arguments passed on to \code{leafs.phylo}
#' @export
#' @name leafs
leafs <- function(x, ...){
  UseMethod('leafs', x)
}

#' @rdname leafs
#' @export
leafs.Strata <- function(x, ...){
  .fmap(x, leafs, fout=.name, ...)
}

#' @rdname leafs
#' @export
leafs.phylo <- function(x, byname=FALSE, ...){
  ids <- 1:length(x$tip.label)
  if(byname){
    x$tip.label[ids]
  } else {
    ids
  }
}

#' Get node indices
#'
#' @param x phylo object
#' @param rebase logical Make indices base 1
#' @param byname logical Return node labels, rather than indices
#' @param ... Arguments passed on to \code{nodes.phylo}
#' @export
#' @name nodes
nodes <- function(x, ...){
  UseMethod('nodes', x)
}

#' @rdname nodes
#' @export
nodes.Strata <- function(x, ...){
  .fmap(x, nodes, fout=.identity, byname=TRUE, ...)
}

#' @rdname nodes
#' @export
nodes.phylo <- function(x, rebase=FALSE, byname=FALSE, ...){
  ids <- (length(x$tip.label)+1):(nleafs(x) + x$Nnode)
  ids_base1 <- ids - nleafs(x)
  if(byname){
    if(is.null(x$node.label))
      stop("Cannot get nodes by name since tree has no node labels")
    x$node.label[ids_base1]
  } else if(rebase){
    ids_base1
  } else {
    ids
  }
}

#' Get the number of leafs
#'
#' @param x phylo object
#' @export
#' @name nleafs
nleafs <- function(x){
  UseMethod('nleafs', x)
}

#' @rdname nleafs
#' @export
nleafs.Strata <- function(x){
  .fmap(x, nleafs, fout=.identity)
}

#' @rdname nleafs
#' @export
nleafs.phylo <- function(x){
  length(leafs(x))
}

#' Find the ancestors of a node
#'
#' @param x phylo object
#' @param id A single name or index
#' @param use_name Should names be used rather than internal ids
#' @param type id type: ['name', 'id', 'auto']
#' @param ... Arguments passed on to \code{lineage.phylo}
#' @examples
#' data(atree)
#' lineage(atree, 't1')
#' lineage(atree, 1)
#' @export
#' @name lineage
lineage <- function(x, ...){
  UseMethod('lineage', x)
}

#' @rdname lineage
#' @export
lineage.Strata <- function(x, id=x@focal_species, type='name', use_name=TRUE, ...){
  .fmap(x, lineage, fout=.identity, id=id, type=type, use_name=use_name)
}

#' @rdname lineage
#' @export
lineage.phylo <- function(x, id, use_name=FALSE, type='auto', ...){
  id <- clean_phyid(x, id, len=1, type=type)
  id <- if(is_root(x, id, type='index')){
    id
  } else {
    c(lineage(x, parent(x, id, type='index'), type='index'), id)
  }
  if(use_name){
    tree_names(x)[id]
  } else {
    id
  }
}

#' Get the number of species and ancestors in the tree
#'
#' @param x phylo or Strata object
#' @return integer
#' @export
#' @name tree_size
tree_size <- function(x){
  UseMethod('tree_size', x)
}

#' @rdname tree_size
#' @export
tree_size.Strata <- function(x){
  .fmap(x, tree_size, fout=.identity)
}

#' @rdname lineage
#' @export
tree_size.phylo <- function(x){
  max(nodes(x))
}

#' Vectorized node parent getter
#'
#' @param x phylo object
#' @param id vector of ids or names
#' @param type id type: ['name', 'id', 'auto']
#' @return vector of parent ids (NA if root)
#' @param ... Arguments passed on to \code{parent.phylo}
#' @export
#' @name parent
parent <- function(x, ...){
  UseMethod('parent', x)
}

#' @rdname parent
#' @export
parent.Strata <- function(x, id=id, type='name', ...){
  .fmap(x, parent, fout=.name, id=id, type=type)
}

#' @rdname parent
#' @export
parent.phylo <- function(x, id=1:tree_size(x), type='auto', ...){
  id <- clean_phyid(x, id, type=type)
  x$edge[match(id, x$edge[, 2]), 1]
}

#' Get the immediate children of a node
#'
#' @param x phylo object
#' @param id vector of ids or names
#' @param type id type: ['name', 'id', 'auto']
#' @param ... Arguments passed on to \code{children.phylo}
#' @return vector of children, (integer(0) if there are no children)
#' @export
#' @name children
children <- function(x, ...){
  UseMethod('children', x)
}

#' @rdname children
#' @export
children.Strata <- function(x, id, type='name', ...){
  .fmap(x, children, fout=.name, id=id, type=type)
}

children.phylo <- function(x, id, type='auto', ...){
  id <- clean_phyid(x, id, len=1, type=type)
  sort(x$edge[x$edge[,1] == id, 2])
}

#' Vectorized tree root finder
#'
#' @param x phylo object
#' @param id vector of ids or names
#' @param type id type: ['name', 'id', 'auto']
#' @param ... Arguments passed on to \code{is_root.phylo}
#' @return logical vector
#' @export
#' @name is_root
is_root <- function(x, ...){
  UseMethod('is_root', x)
}

#' @rdname is_root
#' @export
is_root.Strata <- function(x, ...){
  .fmap(x, is_root, fout=.identity, ...)
}

#' @rdname is_root
#' @export
is_root.phylo <- function(x, id=1:tree_size(x), type='auto', ...){
  id <- clean_phyid(x, id, type=type)
  is.na(parent(x, id, type='index'))
}

#' Get the root of a tree
#'
#' There should normally be only one, but if there are more, this function will
#' find all of them.
#'
#' @param x phylo object
#' @param ... Arguments passed to \code{is_root}
#' @return index of the root node (or nodes if there are more than one)
#' @export
#' @name get_root
get_root <- function(x){
  UseMethod('get_root', x)
}

#' @rdname get_root
#' @export
get_root.Strata <- function(x){
  .fmap(x, get_root, fout=.name)
}

#' @rdname get_root
#' @export
get_root.phylo <- function(x){
  which(is_root(x))
}

#' Get all nodes (leafs and ancestors) descending from, and including, a set of ids
#'
#' @param x phylo object
#' @param id vector of ids or names
#' @param type id type: ['name', 'id', 'auto']
#' @param ... Arguments passed on to \code{descendent_nodes.phylo}
#' @return integer vector of leaf and node ids
#' @export
#' @name descendent_nodes
descendent_nodes <- function(x, ...){
  UseMethod('descendent_nodes', x)
}

#' @rdname descendent_nodes
#' @export
descendent_nodes.Strata <- function(x, id, type='name', ...){
  .fmap(x, descendent_nodes, fout=.name, id=id, type=type, ...)
}

#' @rdname descendent_nodes
#' @export
descendent_nodes.phylo <- function(x, id, type='auto', ...){
  id <- clean_phyid(x, id, type=type)
  children <- x$edge[x$edge[,1] %in% id, 2]
  ids <-
    if(length(children) > 0)
      c(id, descendent_nodes(x, children))
    else
      id
  sort(ids)
}

#' Get all leafs descending from a set of ids
#'
#' This is a simple wrapper around \code{descendent_nodes}.
#'
#' @param x phylo object
#' @param id vector of ids or names
#' @param ... Arguments sent to \code{descendent_nodes}
#' @return integer vector of leaf ids
#' @export
#' @name descendents
descendents <- function(x, ...){
  UseMethod('descendents', x)
}

#' @rdname descendents
#' @export
descendents.Strata <- function(x, ...){
  .fmap(x, descendents, fout=.name, ...)
}

#' @rdname descendents
#' @export
descendents.phylo <- function(x, id, ...){
  descendent_nodes(x, id, ...) %>%
    intersect(leafs(x)) %>%
    sort
}

#' Select a subset of nodes from a tree
#'
#' @param x phylo object
#' @param id vector of ids or names
#' @param collapse logical collapse edges around nodes with a single descendent
#' @param descend logical include all descendents of each id
#' @param type id type: ['name', 'id', 'auto']
#' @param ... Arguments sent to \code{subtree.phylo}
#' @return phylo object
#' @examples
#' data(atree)
#' subtree(atree, c('n15', 'n19'))
#' subtree(atree, c('t7', 't4', 't1'))
#' @export
#' @name subtree
subtree <- function(x, id, ...){
  UseMethod('subtree', x)
}

#' @rdname subtree
#' @export
subtree.Strata <- function(x, id, type='name', ...){
  .fmap(x, subtree, fout=.strata, id=id, type=type, ...)
}

#' @rdname subtree
#' @export
subtree.phylo <- function(x, id, collapse=TRUE, descend=TRUE, type='auto', ...){
  id <- clean_phyid(x, id, type=type)
  if(descend){
    id <- descendent_nodes(x, id, type='index')
  }
  leaf_ids <- intersect(leafs(x), id)
  id <- lapply(leaf_ids, lineage, x=x, type='index') %>%
    do.call(what=c) %>% unique
  new_edge <- x$edge[(x$edge[, 1] %in% id) & (x$edge[, 2] %in% id), , drop=FALSE]
  new_Nnode <- length(unique(new_edge[, 1])) 
  new_size <- length(unique(c(new_edge[, 1], new_edge[, 2])))
  leaf_ids <- intersect(leafs(x), id)
  node_ids <- intersect(nodes(x), id)
  idmap <- seq_along(c(leaf_ids, node_ids))
  names(idmap) <- c(leaf_ids, node_ids)
  new_edge[, 1] <- idmap[as.character(new_edge[, 1])]
  new_edge[, 2] <- idmap[as.character(new_edge[, 2])]
  new_tip.label <- x$tip.label[leaf_ids]
  new_node.label <- x$node.label[node_ids - nleafs(x)]
  new_tree <- list(
    edge        = new_edge,
    tip.label   = new_tip.label,
    edge.length = x$edge.length[x$edge[, 1] %in% id],
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
    new_tree$edge <- matrix(c(2,1), ncol=2)
    new_tree$Nnode <- 1
  } else if(new_size > 1 && collapse){
    # collapse the edges of nodes that have only a single descendent 
    new_tree <- ape::collapse.singles(new_tree)
  }
  new_tree
}

#' Get the sisters of a node
#'
#' @param x phylo object
#' @param id vector of ids or names
#' @param type id type: ['name', 'id', 'auto']
#' @param ... Arguments sent to \code{sisters.phylo}
#' @return integer vector of indices
#' @export
#' @name sisters
sisters <- function(x, id, ...){
  UseMethod('sisters', x)
}

#' @rdname sisters
#' @export
sisters.Strata <- function(x, id, ...){
  .fmap(x, sisters, id=id, ...)
}

#' @rdname sisters
#' @export
sisters.phylo <- function(x, id, type='auto', ...){
  id <- clean_phyid(x, id, len=1, type=type)
  if(is_root(x, id, type='index')){
    integer(0)
  } else {
    setdiff(children(x, parent(x, id)), id)
  }
}

#' Get list of sister trees for a given node
#'
#' @param x phylo object
#' @param id vector of ids or names
#' @param type id type: ['name', 'id', 'auto']
#' @param ... Arguments passed to \code{sister_trees.phylo} 
#' @return list of phylo objects
#' @export
#' @name sister_trees
sister_trees <- function(x, id, ...){
  UseMethod('sister_trees', x)
}

#' @rdname sister_trees
#' @export
sister_trees.Strata <- function(x, id, ...){
  .fmap(x, sister_trees, id=id, ...)
}

#' @rdname sister_trees
#' @export
sister_trees.phylo <- function(x, id, type='auto', ...){
  sister_ids <- sisters(x, id, type=type)
  sister_trees <- lapply(sister_ids, subtree, x=x, type='index')
  if(length(sister_trees) > 0){
    names(sister_trees) <- sister_ids
  }
  sister_trees
}

#' Remove the specified indices, and their descendents, from the tree
#'
#' @param x phylo object
#' @param id vector of indices or names 
#' @param type id type: ['name', 'id', 'auto']
#' @param ... Arguments sent to \code{prune.phylo}
#' @export
#' @name prune
prune <- function(x, id, ...){
  UseMethod('prune', x)
}

#' @rdname prune
#' @export
prune.Strata <- function(x, id, ...){
  .fmap(x, prune, id=id, ...)
}

#' @rdname prune
#' @export
prune.phylo <- function(x, id, type='auto', ...){
  id <- clean_phyid(x, id, type=type)
  id <- descendents(x, id, type='index')
  subtree(x, setdiff(leafs(x), id), type='index')
}

#' Merge fully named subtrees according to a reference tree
#'
#' @param x phylo object
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
#' @name merge_phylo
merge_phylo <- function(x, subtrees){
  UseMethod('merge_phylo', x)
}

#' @rdname merge_phylo
#' @export
merge_phylo.Strata <- function(x, subtrees){
  stop("Not Implemented")
}

#' @rdname merge_phylo
#' @export
merge_phylo.phylo <- function(x, subtrees){
  lapply(subtrees, tree_names) %>%
    do.call(what='c') %>%
    unique %>%
    subtree(x=x, type='name')
}

#' Organize the tips with the focal_id on tip
#'
#' @param x phylo object
#' @param focal_id The name of the taxon to be ordered relative to
#' @name make_tree_relative_to
make_tree_relative_to <- function(x, focal_id){
  UseMethod('make_tree_relative_to', x)
}

#' @rdname make_tree_relative_to
#' @export
make_tree_relative_to.Strata <- function(x, focal_id){
  .fmap(x, make_tree_relative_to, focal_id)
}

#' @rdname make_tree_relative_to
#' @export
make_tree_relative_to.phylo <- function(x, focal_id){
  if(!any(focal_id %in% x$tip.label)){
    stop("'focal_id' is not one of the tips of 'tree'")
  }
  tip_vector <- lapply(lineage(x, focal_id, type='name'), function(i){
    lapply(sister_trees(x, i, type='index'), function(x) x$tip.label) %>%
      unlist %>% unname
  }) %>% unlist
  tip_vector <- c(tip_vector, focal_id)
  if(!setequal(tip_vector, x$tip.label)){
    stop("Unexpected error in getting tip_vector, probably the focal_id argument bad")
  }
  ape::rotateConstr(x, tip_vector)
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
        clean_phyid(a, a$tip.label, type='name'),
        clean_phyid(b, a$tip.label, type='name')
      ), ncol=2))
  # iteratively find the parents of nodes, until root is reached
  while(TRUE){
    last_a <- tail(idmaps, n=1)[[1]][, 1]
    last_b <- tail(idmaps, n=1)[[1]][, 2]
    if(!all(is_root(a, last_a))){
      idmaps <- append(
        idmaps,
        list(unique(matrix(c(
            parent(a, last_a, type='index'),
            parent(b, last_b, type='index')
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
