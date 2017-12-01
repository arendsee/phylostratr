#' Make a tree from the union of lineages
#'
#' @param lineages list of lineages, where each lineage is a data.frame with
#' the columns id and name
#' @param remove_subspecies Should leafs with subtrees be removed. Currently if
#' one lineage uses a specific species as a leaf and a second lineage uses a
#' child of the species, a malformed tree is created where one node is both a
#' tip and a node.
#' @return phylo object with node names and no branch lengths
lineages_to_phylo <- function(lineages, remove_subspecies=FALSE){

  if(remove_subspecies)
    lineages <- Filter(function(x){
      # either the lineage terminates in a species
      (tail(x$rank, 1) == 'species') ||
      # or it is a subspecies, where no other member of the lineage is
      # species-terminal
      (
        sum(
          x[x$rank == 'species', 'id'] == 
          (lapply(lineages, function(y) y$id) %>% unlist %>% unname)
        ) == 1
      )
    }, lineages)

  # get tip names (usually species)
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

  ape::collapse.singles(tree)
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
    edges      = matrix(c(idmap[as.character(from)], idmap[as.character(to)]), ncol=2),
    tip.label  = unname(tip.label),
    node.label = unname(node.label),
    Nnode      = length(node.label)
  )
  class(tree) <- 'phylo'

  tree
}


