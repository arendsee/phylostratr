% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/uniprot-access.R
\name{uniprot_fill_strata}
\alias{uniprot_fill_strata}
\title{Download sequence data for each species in a UniProt-based strata}
\usage{
uniprot_fill_strata(strata, ...)
}
\arguments{
\item{strata}{Strata object where all species are represented in UniProt.}

\item{...}{Additional arguments for \code{uniprot_retrieve_proteome}}
}
\value{
Strata object where 'data' slot is filled with protein FASTA files
}
\description{
Download sequence data for each species in a UniProt-based strata
}
