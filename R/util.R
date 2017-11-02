get_max_hit <- function(d){
  dplyr::group_by(d, .data$qseqid, .data$staxid) %>%  
    dplyr::filter(.data$score == max(.data$score))
}
