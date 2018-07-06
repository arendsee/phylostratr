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
  plot(x, show.node.label=TRUE, ...)
}

#' Transpose a CountMatrix
#'
#' Take the traspose of the matrix and swap the labels
#'
#' @param x CountMatrix object
#' @export
t.CountMatrix <- function(x){
  CountMatrix(x=t(x@x), ylab=x@xlab, xlab=x@ylab)
}

#' Plot a CountMatrix object
#'
#' @param x CountMatrix object
#' @param y Not used
#' @param value_trans A function for transforming the values
#' @param normalize A function for normalizing the matrix
#' @param scheme A color scheme
#' @param ... Additional arguments sent to plot 
#' @export
plot.CountMatrix <- function(x, y=NULL,
  value_trans = identity,
  normalize   = normalize_matrix_by_row,
  scheme      = ggplot2::scale_fill_gradient(low = "grey", high = "red", labels=scales::percent),
  ...
){
  cnt <- x@x
  m <- normalize(cnt)

  cnt <- reshape2::melt(cnt) %>% dplyr::filter(value > 0)
  m <- reshape2::melt(m)

  names(m)[1:3] <- c("a", "b", "value")
  names(cnt)[1:3] <- c("a", "b", "n")

  m$value <- value_trans(m$value)
  cnt$value <- value_trans(cnt$value)

  ggplot2::ggplot() +
    ggplot2::geom_tile(data=m, ggplot2::aes(b,a, fill=value)) +
    ggplot2::xlab(x@xlab) +
    ggplot2::ylab(x@ylab) +
    scheme +
    ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle=270, hjust=0, vjust=1),
        legend.title = ggplot2::element_blank()
    ) +
    ggplot2::geom_text(data=cnt, mapping=ggplot2::aes(x=b, y=a, label=n), size=2)

}
