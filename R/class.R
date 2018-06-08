setOldClass("phylo")
setClass(
  "Strata",
  representation(
    tree = "phylo",
    data = "list",
    focal_species = "character"
  )
)

#' Initialize a Strata object
#'
#' @param focal_species The name of the focal species. This name should match the tip label from the tree.
#' @param tree 'phylo' object
#' @param data List of data. This is where filenames for downloaded proteins,
#' among other things, will be stored.
#' @return Strata object
#' @export
Strata <- function(focal_species, tree, data=list()){
  new('Strata',
    tree = tree,
    data = data,
    focal_species = as.character(focal_species)
  )
}


setClass(
  "CountMatrix",
  representation(
    x = "matrix",
    xlab = "character",
    ylab = "character"
  )
)

#' Initialize a CountMatrix object
#'
#' @param x The matrix
#' @param ylab A label for the y axis (used in plotting, for example)
#' @param xlab A label for the x axis
#' @return CountMatrix object
#' @export
CountMatrix <- function(x, ylab="from", xlab="to"){
  m <- new("CountMatrix") 
  m@x = x
  m@xlab = xlab
  m@ylab = ylab
  m
}
