#' Add proteome stats for each proteome
#'
#' @param strata A Strata object with an 'faa' field
#' @param overwrite logical. If TRUE, then the 'proteome_stats' field will be
#' overwritten if it already exists.
#' @return A Strata object with a new 'proteome_stats' data field
#' @export
add_proteome_stats <- function(strata, overwrite=FALSE){
  is_valid_strata(strata, required='faa')

  if(overwrite || !('proteome_stats' %in% names(strata@data))){
    strata@data$proteome_stats <- lapply(strata@data$faa, function(x){
      seq <- Biostrings::readAAStringSet(x)
      list(n = length(seq), lengths=Biostrings::width(seq))
    })
  }

  strata
}


#' Build a table of proteome stratistics 
#'
#' @param strata Strata object
#' @return data.frame with the fields 'species', 'protein_length', 'mrca', 'ps', and 'index'.
#' @export
proteome_stats_table <- function(strata){
  is_valid_strata(strata)
  if(!('proteome_stats' %in% names(strata@data))){
    strata <- add_proteome_stats(strata)
  }

  strata@data$proteome_stats %>%
    tuplify %>%
    lapply(function(x){
      tibble::data_frame(
        species=rep_len(x$name, length(x$value$length)),
        protein_length=x$value$length
      ) %>%
      dplyr::arrange(.data$protein_length) %>%
      dplyr::mutate(index = seq_along(.data$protein_length))
    }) %>%
    do.call(what=rbind) %>%
    merge(get_phylostrata_map(strata), by='species')
}
