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
  domtblout <- rhmmer::read_domtblout(domtblout_file)
  tblout <- rhmmer::read_tblout(tblout_file)
  unlink(tblout_file)
  unlink(domtblout_file)
  list(domtblout=domtblout, tblout=tblout)
}

make_hmmscan_filter <- function(by='domain_ievalue', k=1e-5, op='<'){
  f <- get(op)
  function(x){
    x$domtblout[f(x$domtblout[[by]], k), ] %>%
      dplyr::select(
        domain = 'domain_accession',
        qseqid = 'query_name',
        dlen   = 'domain_len',
        by,
        qlen   = 'qlen',
        dfrom  = 'hmm_from',
        dto    = 'hmm_to',
        qfrom  = 'ali_from',
        qto    = 'ali_to'
      )
  }
}

strata_hmmscan <- function(
  strata,
  finder  = run_hmmscan,
  filter  = make_hmmscan_filter(by='domain_ievalue', k=1e-5),
  ...
){
  strata@data$hmmscan <- lapply(
    names(strata@data$faa),
    function(species) filter(finder(strata@data$faa[[species]], ...))
  )
  strata
}

strata_load_hmmscan <- function(
  strata,
  dir = 'hmmscan-results',
  filter  = make_hmmscan_filter(by='domain_ievalue', k=1e-5),
  ...
){
  if(!dir.exists(dir)){
    stop("Expected to find hmmscan resuls in folder 'hmmscan-results', but this folder does not exist")
  }
  strata@data$hmmscan <- lapply(strata@tree$tip.label, function(taxid){
    list(
      domtblout = rhmmer::read_domtblout(file.path(dir, paste0(taxid, '.domtblout.tab'))),
      tblout = rhmmer::read_tblout(file.path(dir, paste0(taxid, '.tblout.tab')))
    )
  })
  strata
}

retrievable_pfam_taxids <- function(){
  pfam_url <- "ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam31.0/proteomes/"
  RCurl::getURL(pfam_url, verbose=TRUE, dirlistonly=TRUE) %>%
    readr::read_lines() %>%
    sub(pattern='\\.tsv.gz', replacement='') %>%
    { .[grepl('^[0-9]+$', .)] } %>%
    as.integer
}

strata_retrieve_PFAM_domains <- function(strata, dir='pfam-domains', version='31.0'){
  col_types <- c(
    seq_id          = readr::col_character(),
    alignment_start = readr::col_integer(),
    alignment_end   = readr::col_integer(),
    envelope_start  = readr::col_integer(),
    envelope_end    = readr::col_integer(),
    hmm_acc         = readr::col_character(),
    hmm_name        = readr::col_character(),
    type            = readr::col_character(),
    hmm_start       = readr::col_integer(),
    hmm_end         = readr::col_integer(),
    hmm_length      = readr::col_integer(),
    bit_score       = readr::col_double(),
    E_value         = readr::col_double(),
    clan            = readr::col_character()
  )

  pfam_taxids <- retrievable_pfam_taxids()

  strata@data$pfam <- rep(NA, length(strata@tree$tip.label))
  names(strata@data$pfam) <- strata@tree$tip.label

  dir.create(dir, showWarnings=FALSE)

  strata@data$pfam[names(strata@data$pfam) %in% pfam_taxids] <- 
    lapply(strata@tree$tip.label, function(taxid){
      url=glue::glue('ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam{version}/proteomes/{taxid}.tsv.gz')
      filename <- file.path(dir, paste0(taxid, ".tsv"))
      if(!file.exists(filename)){
        readr::read_tsv(url, comment="#", col_names=col_names) %>%
          readr::write_tsv(path=filename)
      }
      filename
    })

  strata
}
