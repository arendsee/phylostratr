#' Get the uniprot ids downstream of a node
#'
#' @param taxid Ancestral node of clade of interest
#' @param reference_only Should only reference proteomes be considered?
#' @param delay Sleep for 0.3 seconds before retrieving (polite in loops)
#' @return A numeric vector of NCBI taxon ids listing all species in this clade
#' for which Uniprot has a complete proteome.
#' @export
uniprot_downstream_ids <- function(taxid, reference_only=FALSE, delay=FALSE){
  ref_str <- if(reference_only){
    'reference:yes'
  } else {
    '*'
  }
  url_str <- glue::glue(
    'http://www.uniprot.org/taxonomy/?query=ancestor:{taxid}&format=list&fil=proteome:({ref_str})'
  )
  if(delay)
    Sys.sleep(0.3)
  con <- curl::curl(url_str)
  children <- readLines(con) %>% as.integer
  close(con)
  children
}

#' Download a UniProt proteome
#'
#' Thre proteome is written to a file with the name '<taxid>.faa', for example,
#' Arabidopsis thaliana, which has the id 3702, will be written the '3702.faa'.
#'
#' @param taxid An NCBI taxonomy id
#' @param keep_isoforms Should all isoforms of each protein be included?
#' @param dir Directory in which to write all FASTA files
#' @param dryrun If TRUE, do not download genomes, just print the URLs and destination files
#' @param verbose If TRUE, print progress messages
#' @export
#' @examples
#' \dontrun{
#' # uniprot_retrieve_proteome(3702)
#' }
uniprot_retrieve_proteome <- function(
  taxid,
  keep_isoforms = FALSE,
  dir           = 'uniprot-seqs',
  dryrun        = FALSE,
  verbose       = FALSE
){
  if(!dir.exists(dir) && !dryrun){
    dir.create(dir, recursive=TRUE)
  }
  inc_str <- if(keep_isoforms){ 'include=yes' } else { 'include=no' }
  for_str <- 'format=fasta'
  fastafile <- file.path(dir, paste0(taxid, '.faa'))
  if(file.exists(fastafile)){
    maybe_message("Skipping %s - already retrieved", verbose, taxid)
  } else {
    url_str <- glue::glue(
      "http://www.uniprot.org/uniprot/?query=organism:{taxid}&{for_str}&{inc_str}"
    )
    if(dryrun){
      message(sprintf("Checking %s ...", taxid))
      message(sprintf("  url: %s", url_str))
      message(sprintf("  destination: %s", fastafile))
    } else {
      maybe_message("Retrieving %s ...", verbose, taxid)
      # FIXME: if the ID is not in UniProt, no warning is emitted and an empty
      # file is created
      curl::curl_download(url_str, fastafile)
      Sys.sleep(1) # for good manner
    }
  }
  fastafile
}

#' Download sequence data for each species in a UniProt-based strata
#'
#' @param strata data.tree where leaves are named by NCBI taxon ID. These IDs
#' are expected to be in UniProt.
#' @param ... Additional arguments for \code{uniprot_retrieve_proteome}
#' @return A data.tree with the additional field 'protein_file'
#' @export
uniprot_fill_strata <- function(strata, ...){
  data.tree::Do(
    data.tree::Traverse(strata, filterFun=isLeaf),
    function(node){
      node$protein_file <- uniprot_retrieve_proteome(node$name, ...)
    }
  )
  # Since data.tree uses reference semantics, this return is not strictly
  # necessary. But by returning here we can at pretend we aren't sinning.
  strata
}

#' Given a focal taxid, find the uniprot descendents of a taxon's uncles
#'
#' @param taxid The focal species NCBI taxon id
#' @param ... Additional arguments sent to \code{uniprot_downstream_ids}
#' @return data.tree with strata from the NCBI lineage. Each stratum has one or
#' more representative clades, which are the immediate, outgroup children of
#' the stratum's mose recent common ancestor. For each of the 'uncles', all
#' descendent species that are represented in UniProt are added included. Note
#' this is a flat tree, the topology between uncle and descendent is lost.
#' @export
uniprot_cousins <- function(taxid, filter){
  tree <- ncbi_uncles(taxid)
  data.tree::Do(
    data.tree::Traverse(tree, filterFun=data.tree::isLeaf),
    function(node){
      children <- uniprot_downstream_ids(taxid)[-1]
      for(child in children){
        node$AddChild(child)
      }
      node$type = 'uncle'
  })

  is_stratum <- function(node) {
    child_types <- unlist(sapply(node$children, function(n) n$type))
    all(!is.null(child_types) & (child_types == "uncle"))
  }

  data.tree::Do(data.tree::Traverse(tree, filterFun=is_stratum), filter)

  tree
}

#' Retrive the sequences from the uniprot cousins
#'
#' cousin_sets A named list of vectors of species taxon ids. Each vector
#' in the list bears the name of an uncle, all at the same level in the tree.
#' There may be multiple uncles in multifurcating nodes of the tree (this is
#' very common in the NCBI common tree).
#'
#' @param cousins The output of \code{uniprot_cousins}
#' @param dir Directory in which to write sequences
#' @param prefix The phylostratum-level prefix
#' @param dryrun If TRUE, do not download files or create directories
#' @param verbose Print progress messages
#' @param ... Additional arguments passed to \code{uniprot_retrieve_proteome}
#' @return Nothing, this function is run for its effects
#' @examples
#' \dontrun{
#' uniprot_cousins(3702) %>%
#'   lapply(take_first) %>%
#'   uniprot_cousin_proteomes
#'
#' cfilter <- make_do_if_over(3, take_first)
#' uniprot_cousins(3702) %>%
#'   lapply(cfilter) %>%
#'   uniprot_cousin_proteomes
#' }
uniprot_cousin_proteomes <- function(
  cousins,
  dir     = 'strata',
  prefix  = 'ps_',
  dryrun  = FALSE,
  verbose = TRUE,
  ...
){
  strata <- names(cousins)
  for(ps in seq_along(cousins)){
    stratum_str <- strata[ps]
    for(node_id in names(cousins[[ps]])){
      psdir <- file.path(dir, paste0(prefix, ps), node_id)
      if(dryrun){
        maybe_message("Retrieving descendents of '%s' ...", verbose, node_id)
      } else {
        maybe_message("Checking descendents of '%s' ...", verbose, node_id)
        if(!dir.exists(psdir) && length(cousins[[ps]][[node_id]] > 0)){
          dir.create(psdir, recursive=TRUE)
        }
      }
      for(taxid in cousins[[ps]][[node_id]]){
        maybe_message(node_id, verbose)
        uniprot_retrieve_proteome(node_id, dir=psdir, dryrun=dryrun, verbose=verbose, ...)
      }
    }
  }
}


uniprot_sample_prokaryotes <- function(downto='class'){
  # Get all bacterial and Archael classes (class is one level below phylum)
  prokaryote_classes <- taxize::downstream(c('eubacteria', 'Archaea'), downto=downto, db='ncbi')

  # Get all uniprot reference genomes for each bacterial class
  bacteria_class_reps <- lapply(
      prokaryote_classes$eubacteria$childtaxa_id,
      uniprot_downstream_ids, reference_only=TRUE, delay=TRUE
    )
  names(bacteria_class_reps) <- prokaryote_classes$eubacteria$childtaxa_id

  # Get all uniprot reference genomes for each bacterial class
  archaea_class_reps <- lapply(
      prokaryote_classes$Archaea$childtaxa_id,
      uniprot_downstream_ids, reference_only=TRUE, delay=TRUE
    )
  names(archaea_class_reps) <- prokaryote_classes$Archaea$childtaxa_id

  # From each class, randomly select a single uniprot reference genome
  sample_taxids <- function(x, ...){
      # This is required, because R is evil. If x is of length 1 and is numeric,
      # then `sample` treats it as the upperbound of a discrete distribution
      # between 1 and x. Otherwise it is treated as a set to be sampled from.
      sample(as.character(x), ...) %>% as.integer
  }
  bacteria_taxids <- bacteria_class_reps %>%
      Filter(f = function(x){ length(x) > 0}) %>%
      lapply(sample_taxids, size=1) %>% unlist
  archaea_taxids <- archaea_class_reps %>%
      Filter(f = function(x){ length(x) > 0}) %>%
      lapply(as.character) %>%
      lapply(sample_taxids, size=1) %>% unlist

  list(`2`=bacteria_taxids, `2157`=archaea_taxids)
}

#' Use the prokaryotic species that I use
#'
#' This set of species is currently not all that great. It is a sample of
#' bacterial and archaeal species from each class in each domain.
#'
#' @param x A named list of named lists, where all names are NCBI taxon IDs and
#' the leaves are vectors of species-level IDs.
#' @return x with the basal stratum (ID=131567) replaced
#' @export
use_recommended_prokaryotes <- function(x){
  prokaryote_sample <- readRDS(system.file(
    'extdata',
    'prokaryote_sample.rda',
    package='phylostratr'
  ))

  # If the basal stratum is already cellular organisms, replace it
  if(as.integer(names(x)[1]) == 131567L){
    x[[1]] <- prokaryote_sample
  # if not, prepend it
  } else {
    x <- append(list(`131567`=prokaryote_sample), x)
  }
  x
}
