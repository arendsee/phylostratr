proteome_report_table <- function(strata){
  d <- proteome_stats_table(strata)

  d %>%
    dplyr::group_by(.data$species, .data$mrca) %>%
    dplyr::summarize(
      N      = length(.data$ps),
      min    = min(.data$protein_length),
      q25    = quantile(.data$protein_length, probs=0.25),
      mean   = mean(.data$protein_length),
      median = median(.data$protein_length),
      q75    = quantile(.data$protein_length, probs=0.75),
      max    = max(.data$protein_length),
    )
}
