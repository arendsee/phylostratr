% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/util.R
\name{get_max_hit}
\alias{get_max_hit}
\title{Get the maximum hit}
\usage{
get_max_hit(d)
}
\arguments{
\item{d}{A dataframe with the columns 'qseqid', 'staxid', and 'score'. Other
columns may be present, if so, they are left unchanged.}
}
\description{
The resulting data.frame will be complete, with every (qseqid, staxid) pair represented
}
