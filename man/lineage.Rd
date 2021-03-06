% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tree_manipulation.R
\name{lineage}
\alias{lineage}
\alias{lineage.Strata}
\alias{lineage.phylo}
\alias{tree_size.phylo}
\title{Find the ancestors of a node}
\usage{
lineage(x, ...)

\method{lineage}{Strata}(x, id = x@focal_species, type = "name", use_name = TRUE, ...)

\method{lineage}{phylo}(x, id, use_name = FALSE, type = "auto", ...)

\method{tree_size}{phylo}(x)
}
\arguments{
\item{x}{phylo object}

\item{...}{Arguments passed on to \code{lineage.phylo}}

\item{id}{A single name or index}

\item{type}{id type: ['name', 'index', 'auto']}

\item{use_name}{Should names be used rather than internal ids}
}
\description{
Find the ancestors of a node
}
\examples{
data(atree)
lineage(atree, 't1')
lineage(atree, 1)
}
