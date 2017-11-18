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

#' Get a species tree tracing back to origin for the focal species
#'
#' @param x The focal species taxon ID
#' @param FUN A function of x and ... that builds a species tree
#' @param ... Additional arguments passed to FUN
#' @return A tree with species leaves and optional node names
species_to_species_tree <- function(x, FUN, ...){ }

#' Acquire the protein sequence of each species
#'
#' @param x A species tree
#' @param FUN A function of the species tree and ... that retrieves the raw
#' data (usually protein sequences) needed for the phylostratigraphy analysis
#' @param ... Additional arguments passed to FUN
#' @return A tree of raw data
species_tree_to_data_tree <- function(x, FUN, ...){ }

#' Convert sequence tree to result tree
#'
#' @param x Tree of protein FASTA sequences
#' @param FUN A function of the focal species raw data and a single subject
#' species raw data that produces a raw comparison result
#' @param ... Additional arguments passed to FUN
#' @return A tree of raw results
data_tree_to_result_tree <- function(x, FUN, ...){ }

#' Convert tree of raw results to a p-value tree
#'
#' @param x Tree of raw results for pairwise comparisons
#' @param FUN A function of one branch of a tree that assigns a p-value to each leaf
#' @param ... Additional arguments passed to FUN
#' @return A p-value tree
result_tree_to_pvalue_tree <- function(x, FUN, ...){ }

#' Convert tree of p-values to a homolog tree
#'
#' @param x Tree of raw results for pairwise proteome comparison
#' @param FUN A function of one branch of a tree that assigns a binary value to each leaf
#' @param ... Additional arguments passed to FUN
#' @return A tree with binary leaves
result_tree_to_homolog_tree <- function(x, FUN, ...){ }

#' Convert a homolog tree to a phylostrata table
#'
#' @param x Tree of binary homology results
#' @param FUN A function of a a binary homology tree and ... that produces a
#' phylostrata table
#' @param ... Additional arguments passed to FUN
#' @return A phylostrata table
homolog_tree_to_phylostrata <- function(x, FUN, ...){ }
