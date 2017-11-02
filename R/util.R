get_max_hit <- function(d){
  dplyr::group_by(d, qseqid, staxid) %>%  
    dplyr::filter(score == max(score))
}
