#' Strata print generic function
#'
#' @param x Strata object
#' @param ... Additional arguments (unused)
#' @export
print.Strata <- function(x, ...){
  print(glue::glue(
    "
    Strata object
    focal species: {x@focal_name}
    focal id:      {x@focal_id}
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
