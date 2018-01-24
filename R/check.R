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

check_focal_species <- function(x) {
  # TODO: assert that it is a valid NCBI taxonomy ID

  # assert that it has not children, warn if it does
  nkids <- nrow(taxizedb::children(x)[[1]])
  if(nkids > 0){
    warning("The focal species is not a leaf in the tree. This may be fine. But
            if any of the children of the focal species (e.g. subspecies) are
            in the analysis, the program will crash later, since the focal
            species is assumed to be a terminal node. The solution, in this
            case, is to replace the focal species to the appropriate
            subspecies. For example, use 'Apis mellifera mellifera', instead
            'Apis mellifera'")
  }
}
