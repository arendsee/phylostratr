add_proteome_stats <- function(strata){
  is_valid_strata(strata, required='faa')

  strata@data$proteome_stats <- lapply(strata@data$faa, function(x){
    seq <- Biostrings::readAAStringSet(x)
    list(n = length(seq), lengths=Biostrings::width(seq))
  })

  strata
}
