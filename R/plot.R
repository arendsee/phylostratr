#' Plots for the hit distribution of all genes
#'
#' @param hits A table of best hits
#' @param filename Name of the output PDF file
#' @param tree An optional phylo object, if it is not given, a tree will be
#' infrerred from the taxon IDs in the 'staxid' column. This will work only if
#' all these IDs are valid NCBI taxonomy IDs.
#' @param n The number of genes to display per page
#' @param focal_id The focal taxonomy ID (if given, this will be used to order
#' the tree relative with the focal species on top)
#' @param to_name If TRUE, then the tip labels will be converted from taxonomy
#' IDs to scientific names
#' @param scheme Color scheme
#' @export
#' @examples
#' \dontrun{
#' plot_heatmaps(hits, "heatmaps.pdf", tree=strata@tree) 
#'
#' # You can change the colorscheme and cutoffs:
#' funky_scheme <- list(
#'   cutoff = c(1e-20, 1e-5, 1e-3, 1), 
#'   color  = c('blue', 'green', 'yellow', 'darkorange1', 'darkred')
#' )
#' plot_heatmaps(hits, "heatmaps.pdf", tree=strata@tree, scheme=funky_scheme)
#' }
plot_heatmaps <- function(hits, filename, tree=NULL, n=50, focal_id=NULL, to_name=TRUE, scheme=scheme2){
  # This is a wrapper for making the whole PDF
  pdf(filename)
  print(plot_phyloheatmap(
    hits     = hits,
    tree     = tree,
    n        = n,
    focal_id = focal_id,
    to_name  = to_name,
    scheme   = scheme
  ))
  dev.off()
}

# Do not use directly, wrapped by plot_heatmaps
# Preps data for building individual heat maps, partitions data to pages
plot_phyloheatmap <- function(hits, tree=NULL, n=50, focal_id=NULL, to_name=TRUE, ...){
  dat <- base::split(hits, f=factor(hits$qseqid))

  if(is.null(tree)){
    tree <- lineages_to_phylo(taxizedb::classification(dat[[1]]$staxid))
  }
  if(!is.null(focal_id)){
    tree <- make_tree_relative_to(tree, focal_id)
  }
  if(to_name){
    tree$tip.label <- partial_id_to_name(tree$tip.label)
  }
  N <- length(dat)

  lapply(seq(1, N, by=n), function(i){
    indices <- i:min(i+n-1, N)
    stat <- do.call(what=rbind, dat[indices])
    if(to_name){
      stat$staxid <- partial_id_to_name(stat$staxid)
    }
    plot_one_phyloheatmap(tree=tree, stat=stat, ...)
  })
}

# @param scheme - set the color scheme
# @ggtree_args - named list of arguments sent to ggtree::ggtree
# @ggtree_heatmap_args - named list of arguments sent to ggtree::heatmap
# @ggtree_tiplab - a ggtree::geom_tiplab object
plot_one_phyloheatmap <- function(
  tree,
  stat,
  phylostrata = NULL,
  scheme              = scheme2,
  ggtree_args         = list(),
  ggtree_heatmap_args = list(),
  ggtree_tiplab       = ggtree::geom_tiplab(size=1, color="black")
){
  set <- function(x, ...){
    defaults <- list(...)
    for(n in names(x)){
      defaults[[n]] <- x[[n]]
    }
    defaults
  }
  # set ggtree arguments with defaults
  ggtree_args = set(ggtree_args, layout='slanted', ladderize=FALSE)
  # set heatmpa arguments
  ggtree_heatmap_args = set(ggtree_heatmap_args, offset=14, width=6,
                            colnames=TRUE, colnames_angle=-90, hjust=0,
                            font.size=1)
  stat$evalue_bin <- factor(eval_bins(stat$evalue, scheme))
  # if(!is.null(phylostrata)){
  #   id_levels <- dplyr::select(phylostrata, .data$qseqid, .data$ps)
  #   d <- d %>%
  #     merge(id_levels) %>%
  #     dplyr::arrange(-.data$ps, .data$qseqid) %>%
  #     dplyr::select(-.data$ps)
  # }
  name_order <- unique(stat$qseqid)
  stat <- reshape2::dcast(stat, staxid ~ qseqid, value.var='evalue_bin', fill=5)
  rownames(stat) <- stat[, 1]
  stat <- stat[, -1]  
  stat <- stat[, name_order]
  g <- do.call(ggtree::ggtree, append(list(tr=tree), ggtree_args)) + ggtree_tiplab
  g <- do.call(ggtree::gheatmap, append(list(p=g, data=as.matrix(stat)), ggtree_heatmap_args)) +
    ggplot2::scale_fill_manual(
      values = scheme$color,
      labels = c(scheme$cutoff, paste0(tail(scheme$cutoff,1), '+'))
    )
  g
}

scheme1 <- list(
  cutoff = c(1e-20, 1e-5, 1e-3, 1), 
  color  = c('blue', 'green', 'yellow', 'darkorange1', 'darkred')
)

scheme2 = list(
  cutoff = c(1e-100, 1e-20, 1e-5, 1e-1),
  color = c('#0000FF', '#14A9FF', '#7AFDFF', '#CC9500', '#FF1E00')
)

eval_bins <- function(evalue, scheme=scheme2){
  stopifnot((length(scheme$cutoff)+1) == length(scheme$color))
  # treat NA as maximally insignificant
  evalue[is.na(evalue)] <- Inf
  .bincode(evalue, c(-Inf, scheme$cutoff, Inf))
}

#' Plot the ordered lengths of all proteins in all proteomes
#'
#' @param strata Strata object
#' @param normalize Normalize protein index, this loses proteome length, but
#' makes protein distribution more comparable
#' @return ggplot object
#' @export
plot_proteome_lengths <- function(strata, normalize=FALSE){

  d <- proteome_stats_table(strata)

  if(normalize){
    d <- d %>%
      dplyr::group_by(.data$species) %>%
      dplyr::mutate(index = .data$index / length(.data$index)) %>%
      dplyr::ungroup()
  }

  ggplot2::ggplot(d) +
    ggplot2::geom_path(ggplot2::aes_string(x='index', y='protein_length', group='species', color='mrca')) +
    ggplot2::scale_y_continuous(
        trans='log2',
        breaks=scales::trans_breaks('log2', function(y) round(2^y))
    ) +
    ggplot2::xlab("Protein index, ordered by length") +
    ggplot2::ylab("Protein length") +
    ggplot2::ggtitle(sprintf("Lengths of proteomes used in %s phylostratigraph", strata@focal_species))

}

#' Plot summary statistics for all proteomes 
#'
#' @param strata Strata object
#' @export
plot_proteome_stats <- function(strata){
  d <- proteome_report_table(strata) %>%
    reshape2::melt()
  ggplot2::ggplot(d) +
    ggplot2::geom_point(ggplot2::aes_string(x='species', y='value', color='mrca')) +
    ggplot2::facet_wrap("variable", scales='free') +
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank()
    ) +
    ggplot2::xlab("Species") +
    ggplot2::ylab("Protein length") +
    ggplot2::ggtitle(sprintf("Lengths of proteomes used in %s phylostratigraph", strata@focal_species))
}


###### I need to sort through this old code and figure why it was duplicated ...
# #' Plot the max hit scores against each species
# #'
# #' @param d data.frame of the best hits of focal genes against subject species.
# #' It must have the columns [staxid, qseqid, evalue, score, mrca, ps]
# #' @param ... Additional arguments sent to .arrange_hits
# #' @return A named list of ggplot2 objects, with names being the qseqid
# #' @export
# plot_phyloheatmap <- function(d, ...){
#
#   # A lower bound at which to truncate scores
#   lower.bound=0
#   # An upper bound at which to truncate scores
#   upper.bound=Inf
#   if(!is.null(upper.bound))
#     d$score <- ifelse(d$score > upper.bound, upper.bound, d$score)
#   if(!is.null(lower.bound))
#     d$score <- ifelse(d$score < lower.bound, lower.bound, d$score)
#
#   taxidmap <- d %>%
#     dplyr::select(.data$staxid, .data$ps) %>%
#     dplyr::distinct() %>%
#     dplyr::arrange(-.data$ps)
#
#   d <- dplyr::arrange(d, -.data$ps, .data$staxid)
#   # It is essential to turn the taxon IDs into ordered factors, otherwise they
#   # will be interpreted as numbers, and plotted accordingly.
#   d$staxid <- factor(as.character(d$staxid), levels=as.character(taxidmap$staxid))
#
#   stat <- d %>%
#     dplyr::select(.data$qseqid, .data$score, .data$evalue, .data$staxid) %>%
#     merge(taxidmap, by='staxid') %>%
#     dplyr::group_by(.data$qseqid) %>%
#     dplyr::mutate(index=seq_along(.data$qseqid)) %>%
#     dplyr::ungroup() %>%
#     dplyr::mutate(eval_bins=eval_bins(.data$evalue)) %>%
#     {
#       base::split(., f=factor(.$qseqid))
#     } %>%
#     lapply(dplyr::arrange, -.data$ps, .data$staxid)
#
#   nspecies <- length(unique(d$staxid))
#
#   xlines <- taxidmap %>%
#     dplyr::group_by(.data$ps) %>%
#     dplyr::summarize(n = length(.data$ps)) %>%
#     dplyr::arrange(-.data$ps) %>%
#     { Reduce(sum, .$n, accumulate=TRUE) } %>%
#     { .[-length(.)] } %>%
#     { . + 0.5 } %>%
#     as.data.frame %>%
#     magrittr::set_names('xlines')
#
#   lapply(
#     names(stat),
#     function(qseqid) {
#       plot_one_phyloheatmap(stat[[qseqid]], qseqid, xlines)
#     }
#   ) %>% magrittr::set_names(stat$qseqid)
# }
#
# #' Make PDf with multiple phyloheatmap plots per page
# #'
# #' @param d Besthits table
# #' @param file Filename for output PDF
# #' @param width Integer number of plots per row
# #' @param height Integer number of plots per column
# #' @param ... Additional arguments for \code{plot_phyloheatmap}
# #' @export
# make_phyloheatmap_pdf <- function(d, file='phyloheatmap.pdf', width=1, height=5, ...){
#   plots <- plot_phyloheatmap(d, ...)
#   nloci <- length(unique(d$qseqid))
#   pdf(file)
#   for(page.num in 0:((nloci - 1) %/% (width * height))){
#     i <- width * height * page.num + 1
#     j <- min(width * height * page.num + (width*height), length(plots))
#     gridExtra::grid.arrange(grobs=plots[i:j], ncol=width)
#   }
#   dev.off()
# }
#
# #' Plot the max hit scores against each species for one gene
# #'
# #' This is an internal method, a user would normally use \code{plot_phyloheatmap} or
# #' \code{make_phyloheatmap_pdf}.
# #'
# #' @param stats A data.frame with the columns [index, score, eval_bins]
# #' @param qseqid The query gene ID
# #' @param xlines A table specifying where the phylostrata lines should be drawn
# #' @return ggplot2 object
# plot_one_phyloheatmap <- function(stats, qseqid, xlines){
#   stats$index = 1:nrow(stats)
#   ggplot2::ggplot(stats) +
#     ggplot2::geom_point(ggplot2::aes_string(x='index', y='score', color='eval_bins'), size=0.5) +
#     ggplot2::scale_color_manual(values=c('darkred', 'red', 'orange', 'green', 'blue')) +
#     ggplot2::geom_vline(data=xlines, ggplot2::aes_string(xintercept='xlines'), alpha=0.4) +
#     ggplot2::xlim(1, nrow(stats)) +
#     ggplot2::theme_minimal() +
#     ggplot2::ggtitle(qseqid) +
#     ggplot2::theme(
#       legend.position = 'none',
#       axis.title      = ggplot2::element_blank(),
#       axis.ticks      = ggplot2::element_blank(),
#       axis.text       = ggplot2::element_blank(),
#       panel.grid      = ggplot2::element_blank()
#     )
# }

