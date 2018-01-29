#' Add proteome stats for each proteome
#'
#' @param strata A Strata object with an 'faa' field
#' @param overwrite logical. If TRUE, then the 'proteome_stats' field will be
#' overwritten if it already exists.
#' @return A Strata object with a new 'proteome_stats' data field
#' @export
add_proteome_stats <- function(strata, overwrite=FALSE){
  is_valid_strata(strata, required='faa')

  if(overwrite || !('proteome_stats' %in% names(strata@data))){
    strata@data$proteome_stats <- lapply(strata@data$faa, function(x){
      seq <- Biostrings::readAAStringSet(x)
      list(n = length(seq), lengths=Biostrings::width(seq))
    })
  }

  strata
}


#' Build a table of proteome statistics 
#'
#' @param strata Strata object
#' @return data.frame with the fields 'species', 'protein_length', 'mrca', 'ps', and 'index'.
#' @export
proteome_stats_table <- function(strata){
  is_valid_strata(strata)
  if(!('proteome_stats' %in% names(strata@data))){
    strata <- add_proteome_stats(strata)
  }

  strata@data$proteome_stats %>%
    tuplify %>%
    lapply(function(x){
      tibble::data_frame(
        species=rep_len(x$name, length(x$value$length)),
        protein_length=x$value$length
      ) %>%
      dplyr::arrange(.data$protein_length) %>%
      dplyr::mutate(index = seq_along(.data$protein_length))
    }) %>%
    do.call(what=rbind) %>%
    merge(get_phylostrata_map(strata), by='species')
}

#' Add organelle ID lists for each proteome
#'
#' @param strata A Strata object where leafs are NCBI taxa represented in UniProt
#' @param overwrite logical. If TRUE, then the 'proteome_stats' field will be
#' overwritten if it already exists.
#' @return A Strata object with a new 'organelle' data field
#' @export
add_organelle_proteins <- function(strata, overwrite=FALSE){
  is_valid_strata(strata)

  get_organelle_ids <- function(id, ...){
    intersect(
      uniprot_organelle_ids(id, delay=TRUE, ...),
      extract_uniprot_id_from_fasta(strata@data$faa[[id]])
    )
  }

  if(overwrite || !('organelle' %in% names(strata@data))){
    ids <- leafs(strata)
    # TODO: add test is_taxid, I should implement this in taxizedb
    strata@data$organelle <- lapply(ids, function(id){
      list(
        mitochondrion = get_organelle_ids(id, organelle='Mitochondrion'),
        chloroplast   = get_organelle_ids(id, organelle='Chloroplast')
      )
    })
    names(strata@data$organelle) <- ids
  }

  strata
}

#' Build a table of organelle count statistics 
#'
#' @param strata Strata object
#' @param ... Additional arguments passed to \code{add_organelle_proteins}
#' @return data.frame with the fields 'species', 'n_mitochondrial', 'n_chloroplast', 'mrca', and 'ps'.
#' @export
organelle_table <- function(strata, ...){
  strata <- add_organelle_proteins(strata, ...)

  tuplify(strata@data$organelle) %>% {
    data.frame(
      species = sapply(., function(x) x$name),
      n_mitochondrial = sapply(., function(x) length(x$value$mitochondrion)),
      n_chloroplast   = sapply(., function(x) length(x$value$chloroplast))
    )
  } %>%
  merge(get_phylostrata_map(strata), by='species')
}
