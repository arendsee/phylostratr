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
