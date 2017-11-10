#' Arabidopsis thaliana vignette data
#'
#' A data frame with the following columns:
#' 
#' \enumerate{
#'    \item qseqid unique gene identifier for the focal species query
#'    \item staxid subject NCBI taxon ID
#'    \item evalue BLAST e-value for the best hit of qseqid against staxid
#'    \item score raw score for the best hit
#'    \item mrca Most Recent Common Ancestor of the query and subject species
#'    \item ps phylostratum level (where 1 is root)
#' }
#'
#' The Arabidopsis thaliana gene models are from the Araport11 annotation. Only
#' 1/100 of the proteins are used (sampled via the regular expression 'AT.G.99'
#' on the TAIR model IDs). The TAIR model ids have the format 'AT<chromosome
#' number>G<position>.<model_id>', for example, AT3G30270.2 is the second gene
#' model for the locus AT3G30270 and is located on the 3rd chromosome.
#'
#' @format data.frame
"arabidopsis"
