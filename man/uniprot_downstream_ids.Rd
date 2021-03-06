% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/uniprot-access.R
\name{uniprot_downstream_ids}
\alias{uniprot_downstream_ids}
\title{Get the uniprot ids downstream of a node}
\usage{
uniprot_downstream_ids(taxid, reference_only = FALSE, ...)
}
\arguments{
\item{taxid}{Ancestral node of clade of interest}

\item{reference_only}{Should only reference proteomes be considered?}

\item{...}{Additional args (the logical 'delay', currently)}
}
\value{
A numeric vector of NCBI taxon ids listing all species in this clade
for which Uniprot has a complete proteome.
}
\description{
Get the uniprot ids downstream of a node
}
