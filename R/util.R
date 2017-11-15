#' Get the maximum hit 
#'
#' The resulting data.frame will be complete, with every (qseqid, staxid) pair represented
#'
#' @param d A dataframe with the columns 'qseqid', 'staxid', and 'score'
#' @export
get_max_hit <- function(d){
  d %>%
    dplyr::group_by(.data$qseqid, .data$staxid) %>%
    dplyr::filter(.data$score == max(.data$score)) %>%
    # the filter step can lead to multiple hits with equal score, I want just
    # one, so have to pass through distinct
    dplyr::distinct(qseqid, staxid, .keep_all=TRUE)
}

maybe_message <- function(msg, verbose=TRUE, ...){
  if(verbose)
    message(sprintf(msg, ...))
}
