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
    dplyr::ungroup() %>%
    dplyr::distinct(.data$qseqid, .data$staxid, .keep_all=TRUE)
}

maybe_message <- function(msg, verbose=TRUE, ...){
  if(verbose)
    message(sprintf(msg, ...))
}

.check_type <- function(
  m,
  type,
  test   = function(x) { setequal(class(x), type) },
  nframe = sys.nframe()-1,
  place  = if(nframe > 0) { deparse(sys.calls()[[nframe]]) } else { 'global' }
){
  if(!test(m)){
    varname <- deparse(substitute(m)) # NOTE: this has to be outside of glue
    stop(glue::glue(
      "In 'phylostratr::{place}', expected '{name}' to be of class {exp_type} but got '{obs_type}'",
      obs_type = class(m),
      name     = varname,
      place    = place,
      exp_type = type
    ))
  }
}
