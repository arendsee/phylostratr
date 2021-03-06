% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generics.R
\name{plot.CountMatrix}
\alias{plot.CountMatrix}
\title{Plot a CountMatrix object}
\usage{
\method{plot}{CountMatrix}(
  x,
  y = NULL,
  value_trans = identity,
  normalize = normalize_matrix_by_row,
  scheme = ggplot2::scale_fill_gradient(low = "grey", high = "red", labels =
    scales::percent),
  ...
)
}
\arguments{
\item{x}{CountMatrix object}

\item{y}{Not used}

\item{value_trans}{A function for transforming the values}

\item{normalize}{A function for normalizing the matrix}

\item{scheme}{A color scheme}

\item{...}{Additional arguments sent to plot}
}
\description{
Plot a CountMatrix object
}
