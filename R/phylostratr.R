#' @importFrom magrittr "%>%"
#' @importFrom rlang .data
#' @importFrom utils head tail
#' @importFrom grDevices dev.off pdf
utils::globalVariables(c("%>%", "%$%", "."))
NULL

#' phylostratr: Phylostratigraphy execution and analysis
#'
#' This is a stub
#'
#' @section prokaryote_sample:
#'
#' A selection of bacterial and archeael species which have UniProt reference
#' proteoms. These species were selected by randomly choosing one
#' representative from each class (using the \code{uniprot_sample_prokaryotes}
#' function with random seed 42).
#'
#' This sample is used as the root stratum in both the arabidopsis and yeast
#' vignettes.
#'
#' @docType package
#' @name phylostratr
NULL
