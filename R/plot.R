scheme1 <- list(
  cutoff = c(1e-20, 1e-5, 1e-3, 1), 
  color  = c('blue', 'green', 'yellow', 'darkorange1', 'darkred')
)

scheme2 = list(
  cutoff = c(1e-100, 1e-20, 1e-5, 1e-1),
  color = c('#0000FF', '#14A9FF', '#7AFDFF', '#CC9500', '#FF1E00')
)

eval_bins <- function(d, scheme=scheme2){
  stopifnot((length(scheme$cutoff)+1) == length(scheme$color))
  # treat NA as maximally insignificant
  evalue <- d$evalue
  evalue[is.na(evalue)] <- Inf
  .bincode(evalue, c(-Inf, scheme$cutoff, Inf))
}

plot_one_obo_tree <- function(
  tree,
  stat,
  phylostrata = NULL,
  color_scheme=scheme2
){
  d <- stat
  d$evalue_bin <- factor(eval_bins(stat, color_scheme))
  if(!is.null(phylostrata)){
    id_levels <- dplyr::select(.data$phylostrata, .data$qseqid, level=.data$ps)
    d <- d %>%
      merge(id_levels) %>%
      dplyr::arrange(-.data$level, .data$qseqid) %>%
      dplyr::select(-.data$level)
  }
  name_order <- unique(d$qseqid) 
  d <- reshape2::dcast(d, staxid ~ qseqid, value.var='evalue_bin', fill=5)
  rownames(d) <- d[, 1]
  d <- d[, -1]  
  d <- d[, name_order]
  g <- ggtree::ggtree(tree, layout='slanted', ladderize=FALSE) +
    ggtree::geom_tiplab(size=1, color="black")
  g <- ggtree::gheatmap(g, d, offset=14, width=6, colnames=TRUE,
                   colnames_angle=-90, hjust=0,
                   font.size=1) +
    ggplot2::scale_fill_manual(
      values = color_scheme$color,
      labels = c(color_scheme$cutoff, paste0(tail(color_scheme$cutoff,1), '+'))
    )
  g
}

#' Plots for the hit distribution of all genes
#'
#' @param hits A table of best hits
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
#' scheme = list(
#'   cutoff = c(1e-100, 1e-20, 1e-5, 1e-1),
#'   color = c('#0000FF', '#14A9FF', '#7AFDFF', '#CC9500', '#FF1E00')
#' )
#' plot_one_obo_tree(tree, stat, scheme)
#' }
plot_obo_trees <- function(hits, tree=NULL, n=50, focal_id=NULL, to_name=TRUE, scheme=scheme2){
  dat <- base::split(hits, f=factor(hits$qseqid))

  if(is.null(tree)){
    tree <- lineages_to_phylo(taxizedb::classification(dat[[1]]$staxid))
  }
  if(!is.null(focal_id)){
    tree <- make_tree_relative_to(tree, focal_id)
  }
  if(to_name){
    tree$tip.label <- taxizedb::taxid2name(tree$tip.label)
  }
  N <- length(dat)

  lapply(seq(0, N, by=n), function(i){
    indices <- i:min(i+n-1, N)
    stat <- do.call(what=rbind, dat[indices])
    if(to_name){
      stat$staxid <- taxizedb::taxid2name(stat$staxid)
    }
    plot_one_obo_tree(tree, stat, scheme)
  })
}

#' Plot the max hit scores against each species for one gene
#'
#' This is an internal method, a user would normally use \code{plot_obo} or
#' \code{make_obo_pdf}.
#'
#' @param stats A data.frame with the columns [index, score, eval_bins]
#' @param qseqid The query gene ID
#' @param xlines A table specifying where the phylostrata lines should be drawn
#' @return ggplot2 object
plot_one_obo <- function(stats, qseqid, xlines){
  stats$index = 1:nrow(stats)
  ggplot2::ggplot(stats) +
    ggplot2::geom_point(ggplot2::aes_string(x='index', y='score', color='eval_bins'), size=0.5) +
    ggplot2::scale_color_manual(values=c('darkred', 'red', 'orange', 'green', 'blue')) +
    ggplot2::geom_vline(data=xlines, ggplot2::aes_string(xintercept='xlines'), alpha=0.4) +
    ggplot2::xlim(1, nrow(stats)) +
    ggplot2::theme_minimal() +
    ggplot2::ggtitle(qseqid) +
    ggplot2::theme(
      legend.position = 'none',
      axis.title      = ggplot2::element_blank(),
      axis.ticks      = ggplot2::element_blank(),
      axis.text       = ggplot2::element_blank(),
      panel.grid      = ggplot2::element_blank()
    )
}


# @param lower.bound A lower bound at which to truncate scores
# @param upper.bound An upper bound at which to truncate scores
.common <- function(d, lower.bound=0, upper.bound=Inf){
  if(!is.null(upper.bound))
    d$score <- ifelse(d$score > upper.bound, upper.bound, d$score)
  if(!is.null(lower.bound))
    d$score <- ifelse(d$score < lower.bound, lower.bound, d$score)

  taxidmap <- d %>%
    dplyr::select(.data$staxid, .data$ps) %>%
    dplyr::distinct() %>%
    dplyr::arrange(-.data$ps)

  d <- dplyr::arrange(d, -.data$ps, .data$staxid)
  # It is essential to turn the taxon IDs into ordered factors, otherwise they
  # will be interpreted as numbers, and plotted accordingly.
  d$staxid <- factor(as.character(d$staxid), levels=as.character(taxidmap$staxid))

  stat <- d %>%
    dplyr::select(.data$qseqid, .data$score, .data$evalue, .data$staxid) %>%
    merge(taxidmap, by='staxid') %>%
    dplyr::group_by(.data$qseqid) %>%
    dplyr::mutate(index=seq_along(.data$qseqid)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(eval_bins=eval_bins(.data$evalue)) %>%
    {
      base::split(., f=factor(.$qseqid))
    } %>%
    lapply(dplyr::arrange, -.data$ps, .data$staxid)

  list(d=d, stat=stat, taxidmap=taxidmap)
}


#' Plot the max hit scores against each species
#'
#' @param d data.frame of the best hits of focal genes against subject species.
#' It must have the columns [staxid, qseqid, evalue, score, mrca, ps]
#' @param ... Additional arguments sent to .arrange_hits
#' @return A named list of ggplot2 objects, with names being the qseqid
#' @export
plot_obo <- function(d, ...){

  dat <- .common(d)

  nspecies <- length(unique(dat$d$staxid))

  xlines <- dat$taxidmap %>%
    dplyr::group_by(.data$ps) %>%
    dplyr::summarize(n = length(.data$ps)) %>%
    dplyr::arrange(-.data$ps) %>%
    { Reduce(sum, .$n, accumulate=TRUE) } %>%
    { .[-length(.)] } %>%
    { . + 0.5 } %>%
    as.data.frame %>%
    magrittr::set_names('xlines')

  lapply(
    names(dat$stats), 
    function(qseqid) {
      plot_one_obo(dat$stats[[qseqid]], qseqid, xlines)
    }
  ) %>% magrittr::set_names(dat$stats$qseqid)
}

#' Make PDf with multiple obo plots per page
#'
#' @param d Besthits table
#' @param file Filename for output PDF
#' @param width Integer number of plots per row
#' @param height Integer number of plots per column
#' @param ... Additional arguments for \code{plot_obo}
#' @export
make_obo_pdf <- function(d, file='obo.pdf', width=1, height=5, ...){
  plots <- plot_obo(d, ...)
  nloci <- length(unique(d$qseqid))
  pdf(file)
  for(page.num in 0:((nloci - 1) %/% (width * height))){
    i <- width * height * page.num + 1
    j <- min(width * height * page.num + (width*height), length(plots))
    gridExtra::grid.arrange(grobs=plots[i:j], ncol=width)
  }
  dev.off()
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

plot_revenant <- function(
  d,
  cutoff      = 40,
  title       = 'revenant.pdf',
  strange     = TRUE,
  lower.bound = NULL,
  upper.bound = NULL
){
  # FIXME: STUB
}

# plot_constriction <- function(d, max_score=100, min_score=0, ...){
#   m <- constrict(d, ...)
#   m$mrca_name <- get_mrca_names(m)[m$ps] %>% unname
#   m$score <- ifelse(m$score > max_score, max_score, m$score)
#   m$score <- ifelse(m$score > min_score, min_score, m$score)
#   m$species <-
#   ggplot(m) +
#     geom_line(aes(x=species, y=score, group=qseqid, color=mrca))
# }

plot_noise <- function(){
  # FIXME: STUB
}
