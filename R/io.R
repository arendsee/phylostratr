load_blast_data <- function(filename=system.file('araport11_subset9.tab', package='phylostratr')){
  readr::read_tsv(
    filename,
    col_names=c('qseqid', 'evalue', 'bitscore', 'staxid'),
    col_types='cdic'
  ) %>% {
    # TODO: Resolve the multiple taxid nonsense
    #       Currently I just throw them all out
    .[!grepl(';', .$staxid, fixed=TRUE), ] %>% as.integer
  }
}
