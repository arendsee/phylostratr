#' Strata print generic function
#'
#' @param x Strata object
#' @param ... Additional arguments (unused)
#' @export
print.Strata <- function(x, ...){
  print(glue::glue(
    "
    Strata object
    focal species: {x@focal_species}
    data:          {paste(names(x@data), collapse=', ')}
    #species:      {length(x@tree$tip.label)}
    S4 Slots:
      x@tree
      x@data
      x@focal_species
    "
  ))
}
setMethod("show", "Strata",
  function(object) print(object)
)

#' Plot a Strata object
#'
#' @param x Strata object
#' @param ... Additional arguments sent to ggtree
#' @export
plot.Strata <- function(x, ...){
  if(class(x) == "Strata"){
    x <- x@tree
  }
  ggtree::ggtree(x, layout='slanted', ladderize=FALSE) +
    ggtree::geom_tiplab(size=2) +
    ggtree::geom_nodelab(size=2)
}
