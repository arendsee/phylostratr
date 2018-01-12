check_hit_table <- function(d, has_mrca=FALSE, has_ps=FALSE){

  if(!is.data.frame(d)){
    msg <- "The table of hits '%s', should be a data.frame, but found a '%s'"
    stop(sprintf(msg, deparse(substitute(d)), class(d)))
  }

  if(!all(c('qseqid', 'evalue', 'score', 'staxid') %in% names(d))){
    msg <- "The table of hits '%s', should have the columns [qseqid, evalue, score, staxid], but found a [%s]"
    stop(sprintf(msg, deparse(substitute(d)), paste(names(d), collapse=", ")))
  }

  if(has_mrca && !('mrca' %in% names(d))){
    msg <- "This table of hits '%s', needs an 'mrca' field, but this is missing"
    stop(sprintf(msg, deparse(substitute(d))))
  }

  if(has_ps && !('ps' %in% names(d))){
    msg <- "This table of hits '%s', needs an 'ps' field, but this is missing"
    stop(sprintf(msg, deparse(substitute(d))))
  }

  if(nrow(d) == 0){
    warning("The hit table is empty")
  }

}

check_noise <- function(x) {
  # This is a stub
}
