load_blast_data <- function(filename=system.file('extdata', 'araport11', 'araport11_subset9.tab', package='phylostratr')){
  readr::read_tsv(
    filename,
    col_names=TRUE,
    col_types='cddc'
  ) %>% {
    # TODO: Resolve the multiple taxid nonsense
    #       Currently I just throw them all out
    .[!grepl(';', .$staxid, fixed=TRUE), ]
  }
}
