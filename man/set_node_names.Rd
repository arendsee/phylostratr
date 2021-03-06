% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tree_manipulation.R
\name{set_node_names}
\alias{set_node_names}
\alias{set_node_names.Strata}
\alias{set_node_names.phylo}
\title{If nodes are not named, give them default names}
\usage{
set_node_names(x, ...)

\method{set_node_names}{Strata}(x, ...)

\method{set_node_names}{phylo}(x, default_names = paste0("n", nodes(x)), ...)
}
\arguments{
\item{x}{phylo or Strata object}

\item{...}{Arguments passed to \code{set_node_names.phylo}}

\item{default_names}{node names to use if they are missing. Must be a vector
with a number of elements equal to the number of nodes in the tree.}
}
\value{
phylo or Strata object
}
\description{
In phylostratr, all nodes in a tree need to be named, since the names are
used to merge trees.
}
