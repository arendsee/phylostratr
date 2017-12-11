#' Find domains in a protein file
#'
#' @param fastafile The name of a protein FASTA file
#' @param db The name of a HMM database (e.g. PFAM)
#' @param nthreads Number of cores to use (leave two or more threads free)
#' @return A list containing two dataframe: tblout and domtblout (see HMMer manual)
run_hmmscan <- function(
  fastafile,
  db,
  nthreads=1
){
  tblout_file    <- tempfile()
  domtblout_file <- tempfile()
  system2(
    'hmmscan',
    args = c(
      '--tblout', tblout_file,
      '--domtblout', domtblout_file,
      '--cpu', nthreads,
      '--noali',
      db, fastafile
    )
  )
  domtblout <- hmmer_parse_domtblout(domtblout_file)
  tblout <- hmmer_parse_tblout(tblout_file)
  unlink(tblout_file)
  unlink(domtblout_file)
  list(domtblout=domtblout, tblout=tblout)
}

hmmer_parse_tblout <- function(file){
  parse_hmmer_output(file, 'tblout')
}

hmmer_parse_domtblout <- function(file){
  parse_hmmer_output(file, 'domtblout')
}

parse_hmmer_output <- function(file, type){
  column_names <- if(type == 'tblout'){
    c(
      'domain_name',
      'domain_accession',
      'query_name',
      'query_accession',
      'sequecne_evalue',
      'sequence_score',
      'sequence_bias',
      'best_domain_evalue',
      'best_domain_score',
      'best_domain_bis',
      'domain_number_exp',
      'domain_number_reg',
      'domain_number_clu',
      'domain_number_ov',
      'domain_number_env',
      'domain_number_dom',
      'domain_number_rep',
      'domain_number_inc',
      'description'
    )
  } else if(type == 'domtblout') {
    c(
      'domain_name',
      'domain_accession',
      'domain_len',
      'query_name',
      'query_accession',
      'qlen',
      'sequence_evalue',
      'sequence_score',
      'sequence_bias',
      'domain_N',
      'domain_of',
      'domain_cevalue',
      'domain_ievalue',
      'domain_score',
      'domain_bias',
      'hmm_from',
      'hmm_to',
      'ali_from',
      'ali_to',
      'env_from',
      'env_to',
      'acc',
      'description'
    )
  }

  N <- length(column_names)

  readr::read_lines(file) %>%
    Filter(f=function(x) grepl('^[^#]', x)) %>%
    sub(
      pattern = sprintf("(%s) *(.*)", paste0(rep('\\S+', N-1), collapse=" +")),
      replacement = '\\1\t\\2',
      perl = TRUE
    ) %>%
    paste0(collapse="\n") %>%
    readr::read_tsv(col_names=c('X', 'description')) %>%
    tidyr::separate(X, head(column_names, -1), sep=' +')
}
