#' Execute a SPARQL query
#'
#' @param template filename A SPARQL query file that may contain variables that will be replaced with values from `macros`
#' @param macros named vector where names are macros that occur in the template SPARQL query and values are what will replace the macro
#' @return table containing the rows of data returned from the query
query_sparql <- function(template, macros){
  query <- readLines(template) %>%
    stringr::str_replace_all(macros) %>%
    paste0(collapse="\n") %>%
    URLencode %>%
  {
    httr::content(httr::POST(
      url = 'https://sparql.uniprot.org/sparql',
      httr::accept('text/csv'),
      httr::content_type("application/x-www-form-urlencoded"),
      body = glue::glue("query={.}")
    ), show_col_types = FALSE)
  }
}

#' Get the uniprot ids downstream of a node
#'
#' @param taxid Ancestral node of clade of interest
#' @return A numeric vector of NCBI taxon ids listing all species in this clade
#' for which Uniprot has a complete proteome.
#' @export
uniprot_downstream_ids <- function(taxid){

  template <- system.file("sparql", "get-representative-proteome-taxids.rq", package="phylostratr")

  rows <- query_sparql(template, macros=c("TAXID" = as.character(taxid)))

  # The taxon ids are stored as UIDs, e.g., "http://purl.uniprot.org/taxonomy/3711"
  # So we need to strip out everything but the number and convert to int
  sub("http.*/", "", rows$taxon) %>% as.integer
}

#' Get the uniprot ids for organelle proteins of a given taxid
#'
#' @param taxid Taxon ID for species of interest 
#' @param organelle string, one of 'mitochondrion', 'chloroplast'
#' @return A numeric vector of Uniprot ids
#' @export
uniprot_organelle_ids <- function(taxid, organelle='Mitochondrion'){
  template <- system.file("sparql", "get-all-mitochondrial-uniprot-ids.rq", package="phylostratr")
  macros <- c(TAXID=taxid, ORGANELLE=organelle)
  rows <- query_sparql(template, macros)
  sub("http.*/", "", rows$protein)
}


list_content <- function(con){
  stringi::stri_split_lines(stringi::stri_trim_both(httr::content(con, encoding="UTF-8")))[[1]]
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
  refs <- uniprot_downstream_ids(clade)
  weights <- rep(1.1, length(refs))
  names(weights) <- refs
  weights 
}

#' Retrieve UniProt proteome
#'
#' The proteome is written to a file with the name '<taxid>.faa', for example,
#' Arabidopsis thaliana, which has the id 3702, will be written the '3702.faa'.
#'
#' @param taxid An NCBI taxonomy id
#' @param dir Directory in which to write all FASTA files
#' @param verbose Be super chatty
#' @export
#' @examples
#' \dontrun{
#' # uniprot_retrieve_proteome(3702)
#' }
uniprot_retrieve_proteome <- function(taxid, dir = 'uniprot-seqs', verbose = TRUE){
  if(!dir.exists(dir)){
    dir.create(dir, recursive=TRUE)
  }
  fastafile <- file.path(dir, paste0(taxid, '.faa'))
  if(file.exists(fastafile)){
    maybe_message("Skipping %s - already retrieved", verbose, taxid)
  } else {
    uniprot_retrieve_proteome_table(taxid) %>%
    {
      paste0(">", .$uniprot_uid, "\n", .$aa_sequence, collapse="\n")
    } %>%
      writeLines(fastafile)
  }
  fastafile
}

#' Retrieve UniProt proteome table
#'
#' Tries to retrieve all proteins in a "representative" proteome for the given
#' species. If the species does not have a representative proteome, retrieve
#' all proteins for the species.
#'
#' @param taxid An NCBI taxonomy id
#' @export
#' @examples
#' \dontrun{
#' x <- uniprot_retrieve_proteome_table(3702)
#' }
uniprot_retrieve_proteome_table <- function(taxid){
  template <- system.file("sparql", "get-representative-proteome-sequences.rq", package="phylostratr")
  macros <- c(TAXID=as.character(taxid))
  rows <- query_sparql(template, macros)

  if(nrow(rows) == 0){
    warning(glue::glue("Could not find representative proteome for {taxid}, retrieving all records, may be incomplete or contain redundant entries"))
    template <- system.file("sparql", "get-all-taxid-sequences.rq", package="phylostratr")
    macros <- c(TAXID=as.character(taxid))
    rows <- query_sparql(template, macros)
  }

  rows$uniprot_uid <- sub("http.*/", "", rows$uniprot_uid)
  rows
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
uniprot_strata <- function(taxid, from=2, 
                           drop.names=NULL){ # Jan 3 2025, LTC: add drop.names, a character vector containing taxa ids to completely exclude

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

  # edited by LTC Dec 9 2024 to make this function work when from = 1
  if (from==1) {
    # using the "celluar organisms" taxon id (131567) doesn't work in the sparql query for sublcasses
    # instead, just manually search the downstream ids of each of the 3 sub-classes of cellular organisms:
     c(uniprot_downstream_ids("2"), # bacteria
      uniprot_downstream_ids("2157"), # archaea
      uniprot_downstream_ids("2759") # eukarya
      ) %>%
      add_focal %>%
      taxizedb::classification() %>%
      Filter(f=is.data.frame) %>%
      lineages_to_phylo(clean=TRUE) %>%
      ape::drop.tip.phylo(drop.names) %>% # completely remove some taxa
      Strata(
        focal_species = taxid,
        tree       = .,
        data       = list()
      )
  } else {
  taxizedb::classification(taxid)[[1]]$id[from] %>%
    uniprot_downstream_ids %>%
    add_focal %>%
    taxizedb::classification() %>%
    Filter(f=is.data.frame) %>%
    lineages_to_phylo(clean=TRUE) %>%
    ape::drop.tip.phylo(drop.names) %>% # completely remove some taxa
    Strata(
      focal_species = taxid,
      tree       = .,
      data       = list()
    )
  }
}

#' Map UniProt IDs for an organism to PFAM IDs
#'
#' @param taxid Species NCBI taxonomy ID
#' @return A data.frame with columns 'uniprotID' and 'pfamID'
#' @export
uniprot_map2pfam <- function(taxid){
  template <- system.file("sparql", "get-representative-proteome-pfam-map.rq", package="phylostratr")
  query_sparql(template, macros=c("TAXID" = as.character(taxid))) %>%
    dplyr::mutate(
      uniprotID = sub("http.*/", "", uniprotID),
      pfamID = sub("http.*/", "", pfamID)
    )
}

#' Randomly sample prokaryotic representatives
#'
#' This function was used to select the prokaryotes used in the paper. I do not
#' recommend using it now though. It can select unclassified taxa.
#'
#' @param downto the lowest phylogenetic rank that shoul be sampled
#' @param remake whether to replace an existing file
#' @return phylo object containing the prokaryptic sample tree
# add weights to uniprot_sample_prokaryotes: Jan 6 2025, LTC
# add ability to drop taxa entirely (e.g., with weights==0 by drop.names=names(my.weights)[my.weights==0]) Jan 29 2025 LTC
uniprot_sample_prokaryotes <- function(downto='class', weights=NULL, drop.names=NULL){
  
  # Get all bacterial and Archael classes (class is one level below phylum)
  prokaryote_classes <- taxizedb::downstream(c('eubacteria', 'Archaea'), downto=downto, db='ncbi')
  
  # Get all uniprot reference genomes for each bacterial class
  bacteria_class_reps <- lapply(
    prokaryote_classes$eubacteria$childtaxa_id,
    uniprot_downstream_ids
  )
  names(bacteria_class_reps) <- prokaryote_classes$eubacteria$childtaxa_id
  
  # Get all uniprot reference genomes for each bacterial class
  archaea_class_reps <- lapply(
    prokaryote_classes$Archaea$childtaxa_id,
    uniprot_downstream_ids
  )
  names(archaea_class_reps) <- prokaryote_classes$Archaea$childtaxa_id
  
  clean_reps <- function(taxa){
    # taxa should be a list of clades (class level by default) that each have
    # a list of NCBI taxonomy ids that have representative uniprot proteomes
    taxa %>%
      # remove any clades that have no representative
      Filter(f = function(x) length(x) > 0) %>%
      # get the lineage for each representative
      lapply(taxizedb::classification) %>%
      # remove any cases where lineage was missing
      lapply(Filter, f = function(x) is.data.frame(x)) %>%
      # remove all unclassified entries
      lapply(function(x) Filter(f = function(lineage){ ! any(grepl("unclassified", lineage$name)) }, x)) %>%
      # # remove any clades where all representatives are unclassified
      Filter(f = function(x) length(x) > 0) %>%
      # convert back from classification table to species/strain taxonomy id
      lapply(function(xs) sapply(xs, function(x) x$id[nrow(x)]) %>% unname)
  }
  
  # From each class, randomly select a single uniprot reference genome
  sample_taxids <- function(x, ...){
    # This is required, because R is evil. If x is of length 1 and is numeric,
    # then `sample` treats it as the upperbound of a discrete distribution
    # between 1 and x. Otherwise it is treated as a set to be sampled from.
    sample(as.character(x), ...) %>% as.integer
  }

  # From each class, use weights to select a single uniprot reference genome
sample_taxids_weights <- function(taxid_list, weights) {
  # Use this to choose from each clade using weights
  new.taxid_list <- vector("list", length=length(taxid_list))
  for (i in 1:length(taxid_list)) {
    # pull weights for this clade
    current.taxid <- taxid_list[[i]]
    current.weights <- weights[current.taxid]
    
    # if weight not available, default = 1
    if(sum(is.na(current.weights))>0) {
      na.ids <- is.na(current.weights)
      current.weights[na.ids] <- 1
      names(current.weights)[na.ids] <- current.taxid[na.ids]
    }
    
    stopifnot(names(current.weights)==current.taxid)
    
    # and use them to "sample"
    new.taxid_list[[i]] <- taxid_list[[i]][which.max(current.weights)]
  }
  return(new.taxid_list) 
}
   # clean taxa and remove names to drop:
  bacteria_taxids <- bacteria_class_reps %>%
    clean_reps %>%
    lapply(function(x) {x <- x[!(x %in% drop.names)]}) %>% # drop excluded taxa
    Filter(f = function(x) length(x) > 0) # remove any clades where all representatives are excluded
 
  archaea_taxids <- archaea_class_reps %>%
    clean_reps %>%
    lapply(function(x) {x <- x[!(x %in% drop.names)]}) %>% # drop excluded taxa
    Filter(f = function(x) length(x) > 0) # remove any clades where all representatives are excluded
                                 
     # Use this to choose from each clade using provided weights
  if(is.null(weights)) {
    # if no weights provided, just choose at random
    bacteria_taxids <- bacteria_taxids %>%
      lapply(sample_taxids, size=1) %>% unlist
    
    archaea_taxids <- archaea_taxids %>%
      lapply(sample_taxids, size=1) %>% unlist
    
  } else {
    # otherwise, incorporate the provided weights
    bacteria_taxids <- bacteria_taxids %>%
      sample_taxids_weights(weights) %>% unlist
    
    archaea_taxids <- archaea_taxids %>%
      sample_taxids_weights(weights) %>% unlist
  }
  
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
