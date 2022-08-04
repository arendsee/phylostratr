#' Get the uniprot ids downstream of a node
#'
#' @param taxid Ancestral node of clade of interest
#' @param reference_only Should only reference proteomes be considered?
#' @param ... Additional args (the logical 'delay', currently)
#' @return A numeric vector of NCBI taxon ids listing all species in this clade
#' for which Uniprot has a complete proteome.
#' @export
uniprot_downstream_ids <- function(taxid, reference_only=FALSE, ...){
  ref_str <- if(reference_only){
    'reference:yes'
  } else {
    '*'
  }
  query <- glue::glue('ancestor:{taxid}&fil=proteome:({ref_str})')
  wrap_uniprot_id_retrieval(db='taxonomy', query=query, cast=as.integer, ...)
}

#' Get the uniprot ids for organelle proteins of a given taxid
#'
#' @param taxid Taxon ID for species of interest 
#' @param organelle string, one of 'mitochondrion', 'chloroplast'
#' @param ... Additional args (the logical 'delay', currently)
#' @return A numeric vector of NCBI taxon ids
#' @export
uniprot_organelle_ids <- function(taxid, organelle='Mitochondrion', ...){
  query <- glue::glue('organelle:{organelle}+organism_id:{taxid}')
  wrap_uniprot_id_retrieval(db='uniprotkb', query=query, ...)
}

list_content <- function(con){
  stringi::stri_split_lines(stringi::stri_trim_both(httr::content(con, encoding="UTF-8")))[[1]]
}

paginate <- function(url, parser, reducer){
  con <- httr::GET(url)
  values <- list()
  values[[1]] <- parser(con)
  while(TRUE){
    link_str <- con$all_headers[[1]]$headers$link
    if(is.null(link_str)){
      break
    } else {
      link_url <- sub(x=link_str, pattern="<(.*)>.*", replacement="\\1", perl=TRUE)
      con <- httr::GET(link_url)
      new_values <- list()
      new_values[[1]] <- parser(con)
      values <- append(values, new_values)
    }
  }
  Reduce(reducer, values)
}

#' Internal function for wrapping ID retrieval from UniProt
#'
#' @param db UniProt database to search
#' @param query UniProt query string
#' @param delay Sleep for 0.3 seconds before retrieving (polite in loops)
#' @param cast A function for type casting the resulting IDs (e.g. as_integer)
#' @param date Retrieve only IDs last modified before this date (e.g. "20180601")
#' @return vector of IDs or other literals (nothing structured)
wrap_uniprot_id_retrieval <- function(db, query, date=NULL, delay=FALSE, cast=identity){
  date <- if(is.null(date)){
    ""
  } else {
    glue::glue("created%3A%5B19860101+TO+{date}%5D+AND+")
  }
  url <- glue::glue('https://rest.uniprot.org/{db}/search?query={date}{query}&format=list&size=500')
  message(glue::glue("Accessing uniprot: {url}"))
  if(delay)
    Sys.sleep(0.3)
  cast(paginate(url, parser=list_content, reducer=c))
}

#' Parse a UniProt ID from a FASTA file
#'
#' UniProt FASTA files have the initial pattern:
#'   <db>|<UniProt ID>|<alternative ID> <desc>'
#'
#' Here we parse out the UniProt ID. If the FASTA header does not match the
#' pattern, we return NULL
#'
#' @param fastafile A protein FASTA file, assumed to be UniProt format
#' @return character vector of UniProt IDs IF this is a UniProt FASTA, NULL otherwise
#' @export
extract_uniprot_id_from_fasta <- function(fastafile){
  headers <- names(Biostrings::readAAStringSet(fastafile))
  if(all(grepl('^..\\|', headers))){
    sub('^...([^\\|]+).*', '\\1', headers) 
  } else {
    # This is not a UniProt FASTA file
    NULL
  }
}

#' Make reference-species preferring weight vector for diverse_subtree
#'
#' @param weight How much to prefer the reference species. The default, 1.05,
#' will weakly prefer them, acting mostly as a tie-breaker. Higher weights
#' could lead to reduced diversity.
#' @param clade Weight taxa descending from this clade. The default is 2759
#' (Eukaryota).
#' @export
uniprot_weight_by_ref <- function(weight=1.05, clade=2759){
  refs <- uniprot_downstream_ids(clade, reference_only=TRUE)
  weights <- rep(1.1, length(refs))
  names(weights) <- refs
  weights 
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

  if(! keep_isoforms){
    warning("You asked to remove isoforms from the proteome, but I'm afraid I don't know how to do that under Uniprot's new API.")
  }

  for_str <- 'format=fasta'
  fastafile <- file.path(dir, paste0(taxid, '.faa'))
  if(file.exists(fastafile)){
    maybe_message("Skipping %s - already retrieved", verbose, taxid)
  } else {

    url_str <- glue::glue(
      "https://rest.uniprot.org/proteomes/search?query=organism_id:{taxid}+proteome_type:1"
    )

    # This object contains a lot of potentially useful data
    proteome_data <- httr::content(httr::GET(url_str, httr::accept_json()))
    # After downloading, we can check that the proper number of items have been retrieved
    proteome_size <- proteome_data[[1]][[1]]$proteinCount
    # The proteome id, for example "UP000006548"
    proteome_id <- proteome_data[[1]][[1]]$id

    url_str <- glue::glue("https://rest.uniprot.org/uniprotkb/search?format=fasta&query=proteome:{proteome_id}")

    message(glue::glue("Retrieving proteome from uniprot with: {url_str}"))
    if(dryrun){
      message(sprintf("Checking %s ...", taxid))
      message(sprintf("  url: %s", url_str))
      message(sprintf("  destination: %s", fastafile))
    } else {
      maybe_message("Retrieving %s ...", verbose, taxid)
      url_str <- glue::glue("{url_str}&size=500")

      if(verbose){
        progress <- txtProgressBar(min = 0, max = ceiling(proteome_size / 500), style=3)
        progress_index <- 0
      }

      parse_fasta <- function(con){
        fasta <- httr::content(con, encoding="UTF-8")
        if(verbose){
          progress_index <<- progress_index + 1
          setTxtProgressBar(pb=progress, value=progress_index)
        }
        fasta
      }

      fasta_entries <- paginate(url=url_str, parser=parse_fasta, reducer=c)

      if(verbose){
        close(progress)
      }
      cat(fasta_entries, file=fastafile, sep = "")
    }
  }
  fastafile
}

#' Download sequence data for each species in a UniProt-based strata
#'
#' @param strata Strata object where all species are represented in UniProt.
#' @param ... Additional arguments for \code{uniprot_retrieve_proteome}
#' @return Strata object where 'data' slot is filled with protein FASTA files
#' @export
uniprot_fill_strata <- function(strata, ...){
  species <- strata@tree$tip.label
  strata@data$faa <- lapply(strata@tree$tip.label, uniprot_retrieve_proteome, ...)
  names(strata@data$faa) <- strata@tree$tip.label
  strata
}

#' Given a focal taxid, build a UniProt-based Strata object
#'
#' @param taxid The focal species NCBI taxon id
#' @param from stratum to begin from, where 1 is 'cellular organisms'
#' @return Strata object
#' @export
uniprot_strata <- function(taxid, from=2){

  # ensure the focal gene is included, even if not in uniprot (fix #10)
  add_focal <- function(xs){
    if(!(taxid %in% xs)){
      message("The focal species is not present in UniProt. ",
              "You may add it after retrieving uniprot sequences ",
              "(i.e. with 'uniprot_fill_strata') with a command such as: ",
              "strata_obj@data$faa[[focal_taxid]] <- '/path/to/your/focal-species.faa'")
      xs <- c(taxid, xs)
    }
    xs
  }

  taxizedb::classification(taxid)[[1]]$id[from] %>%
    uniprot_downstream_ids %>%
    add_focal %>%
    taxizedb::classification() %>%
    Filter(f=is.data.frame) %>%
    lineages_to_phylo(clean=TRUE) %>%
    Strata(
      focal_species = taxid,
      tree       = .,
      data       = list()
    )
}

#' Map UniProt IDs for an organism to PFAM IDs
#'
#' @param taxid Species NCBI taxonomy ID
#' @return A data.frame with columns 'uniprotID' and 'pfamID'
#' @export
uniprot_map2pfam <- function(taxid){
  base='https://rest.uniprot.org/uniprotkb'
  format='format=tsv'
  columns='fields=accession,xref_pfam' # see https://www.uniprot.org/help/return_fields for new field list
  url <- glue::glue('{base}/search?query=organism_id:{taxid}&{format}&{columns}&size=500')
  message(glue::glue("uniprot_map2pfam url: {url}"))

  csv_content <- function(x) {
    readr::read_tsv(
      I(httr::content(x, encoding="UTF-8")),
      col_names=c("uniprotID", "pfamID"),
      col_types="cc"
    ) %>%
    dplyr::mutate(pfamID = sub(';$', '', .data$pfamID)) %>%
    tidyr::separate_rows(.data$pfamID, sep=';')
  }

  paginate(url, parser=csv_content, reducer=rbind)
}


uniprot_sample_prokaryotes <- function(downto='class', remake=FALSE){

  old.file <- system.file('extdata', 'prokaryote_sample.rda', package='phylostratr')
  if(!remake && file.exists(old.file)){
    return(readRDS(old.file))
  }

  # Get all bacterial and Archael classes (class is one level below phylum)
  prokaryote_classes <- taxizedb::downstream(c('eubacteria', 'Archaea'), downto=downto, db='ncbi')

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

  c(bacteria_taxids, archaea_taxids) %>%
    taxizedb::classification() %>%
    lineages_to_phylo
}

#' Use the prokaryotic species that I use
#'
#' This set of species is currently not all that great. It is a sample of
#' bacterial and archaeal species from each class in each domain.
#'
#' @param x Strata object
#' @return x Strata object with the basal stratum (ID=131567) replaced
#' @export
use_recommended_prokaryotes <- function(x){
  prokaryote_sample <- readRDS(system.file(
    'extdata',
    'prokaryote_sample.rda',
    package='phylostratr'
  ))
  # extract the Eukaryota branch
  x@tree <- subtree(x@tree, '2759', type='name')
  # get 'cellular_organism -> Eukaryota' tree (just these two nodes)
  root <- taxizedb::classification(2759)[[1]]$id %>% lineage_to_ancestor_tree
  # bind the prokaryotes to 'cellular_organism'
  y <- ape::bind.tree(root, prokaryote_sample)
  # bind the eukaryote input tree to 'Eukaryota'
  y <- ape::bind.tree(y, x@tree, where=which(tree_names(y) == '2759'))
  x@tree <- y
  x
}

#' Extract an ID map from uniprot fasta files
#'
#' This assumes the fasta files respect the UniProt fasta header format.
#'
#' @param strata Strata object 
#' @return Strata object with 'idmap' data field
#' @export
uniprot_add_idmap <- function(strata){
  is_valid_strata(strata, required='faa')
  strata@data$idmap <- lapply(strata@data$faa, function(x){
     Biostrings::readAAStringSet(x) %>%
       names %>%
       strsplit('\\|') %>%
       purrr::transpose() %>%
       lapply(unlist) %>%
       {tibble::data_frame(
         db = .[[1]],
         uniprot_id = .[[2]],
         other_id = sub(" .*", "", .[[3]])
       )}
  })
  names(strata@data$idmap) <- names(strata@data$faa)
  strata
}
