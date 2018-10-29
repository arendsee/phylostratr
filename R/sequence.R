#' Select ids from a FASTA file, creating a new file
#'
#' @param f filename of FASTA protein sequence file
#' @param ids a vector for subseting (either logical with one element for each
#' gene, or integers with 1 <= i <= length(seqs).
#' @param prefix a file basename prefix for the output FASTA file
#' @return filename
#' @export
subset_fasta <- function(f, ids, prefix="zzz_"){
  seq <- Biostrings::readAAStringSet(f)
  seq <- seq[ids]
  newFilename <- file.path(dirname(f), paste0(prefix, basename(f)))
  Biostrings::writeXStringSet(seq, newFilename)
  newFilename
}

#' Select n entries from a FASTA file, creating a new file
#'
#' Selects a total of n entries from a FASTA file by selecting every entry where
#'
#' i % n == n / 2
#'
#' Where i is the element index (1,2, ... N), N is the total number of
#' sequences in the genome, '%' is the modulo operator, and '/' is integer
#' division.
#'
#' The final term 'n / 2' ensures the first and last selected entries are
#' equidistant (for even values of n) from the beginning and end of the
#' sequences.
#'
#' @param f filename of FASTA protein sequence file
#' @param n integer number of entries to extract
#' @param prefix a file basename prefix for the output FASTA file
#' @return filename
#' @export
thin_fasta <- function(f, n, ...){
  seq <- Biostrings::readAAStringSet(f)
  ids <- (seq_along(seq) %% n) == (n %/% 2)
  subset_fasta(f, ids, ...)
}
