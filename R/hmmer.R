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

#' Retrieve uniprot to pfam maps for all taxa
#'
#' @param strata Strata object
#' @param dir director to which maps should be written
#' @return Strata object with uniprot2pfam field in data, which stores the
#'         filenames of the retrieved maps.
#' @export
strata_uniprot_pfam_map <- function(strata, dir='pfam'){
  dir.create(dir, recursive=TRUE, showWarnings=FALSE)
  taxa <- strata@tree$tip.label
  strata@data$uniprot2pfam <- lapply(taxa, function(taxid){
    filename <- file.path(dir, paste0(taxid, ".tab"))
    if(!file.exists(filename))
      readr::write_tsv(uniprot_map2pfam(taxid), path=filename)
    filename
  })
  names(strata@data$uniprot2pfam) <- taxa
  strata
}

stratify_by_pfam_domain <- function(strata){
  if(!('uniprot2pfam' %in% names(strata@data))){
    stop("You need to run `strata_uniprot_pfam_map` before this function")
  }

  reps <- strata@tree %>%
    lineage(strata@focal_species, type='name') %>%
    lapply(sister_trees, tree=strata@tree, type='index') %>%
    lapply(function(x){ lapply(x, leafs, byname=TRUE) %>% unlist %>% unique })

  names(reps) <- tree_names(strata@tree)[lineage(strata@tree, strata@focal_species)]
  reps <- reps[-1]

  domains <- lapply(reps, function(taxa){
    lapply(taxa, function(taxid){
      readr::read_tsv(strata@data$uniprot2pfam[[as.character(taxid)]])$pfamID
    }) %>% unlist %>% unique
  })

  focal_domains <- readr::read_tsv(strata@data$uniprot2pfam[[strata@focal_species]])

  domstrat <- lapply(domains, function(x) {
    focal_domains$pfamID %in% x
  }) %>%
    do.call(what=cbind) %>%
    apply(1, which) %>%
    lapply(head, 1) %>% unlist

  # TODO: the above just gets that stratification of specific domains, I want
  # to stratify genes BY domain. Anyway, tomorrow, continue from here.
    
}

merge_pfam <- funtion(strata){
  pfam <- lapply(seq_along(strata@data$uniprot2pfam), function(i){
    d <- readr::read_tsv(strata@data$uniprot2pfam[[i]], col_types='cc')
    d$staxid <- names(strata@data$uniprot2pfam)[i]
    d
  }) %>%
    do.call(what=rbind) %>%
    dplyr::filter(!is.na(pfamID))
}
