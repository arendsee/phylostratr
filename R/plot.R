eval_bins <- function(evalue){
  bins <- c("e >= 1", "e < 1", "e < 1e-3", "e < 1e-5", "e < 1e-20")
  bin <- ifelse(evalue >= 1e-20, bins[4], bins[5])
  bin <- ifelse(evalue >= 1e-5,  bins[3], bin)
  bin <- ifelse(evalue >= 1e-3,  bins[2], bin)
  bin <- ifelse(evalue >= 1,     bins[1], bin)
  bin <- factor(bin, levels=bins)
  bin
}

plot_one_obo_tree <- function(tree, stat, cutoff=1e-20){
  d <- stat
  d$evalue <- ifelse(is.na(d$evalue), 10, d$evalue)
  d$nlogE <- ifelse(d$evalue < cutoff, cutoff, d$evalue) %>% {-1 * log10(.)}
  d <- reshape2::dcast(d, staxid ~ qseqid, value.var='nlogE')
  rownames(d) <- d[, 1]
  d <- d[, -1]
  d <- sapply(d, function(x) ifelse(is.na(x), -1, x)) %>%
    magrittr::set_rownames(rownames(d))
  g <- ggtree::ggtree(tree, layout='slanted') +
    ggtree::geom_tiplab(size=2, color="black")
  ggtree::gheatmap(g, d, offset=4, width=0.8, colnames=TRUE, colnames_angle=-45, hjust=0)
}

plot_obo_trees <- function(hits, tree=NULL, n=30, focal_id=NULL, to_name=FALSE){
  dat <- .common(hits)
  if(is.null(tree)){
    tree <- lineages_to_phylo(taxizedb::classification(unique(hits$staxid)))
  }
  if(!is.null(focal_id)){
    if(!any(focal_id %in% tree$tip.label)){
      stop("'focal_id' is not one of the tips of 'tree'")
    }
    tip_vector <- lapply(lineage(tree, focal_id, type='name'), function(i){
      lapply(sister_trees(tree, i, type='index'), function(x) x$tip.label) %>% unlist %>% unname
    }) %>% unlist
    tip_vector <- c(tip_vector, focal_id)
    if(!setequal(tip_vector, tree$tip.label)){
      stop("Unexpected error in getting tip_vector, probably the focal_id argument bad")
    }
    tree <- ape::rotateConstr(tree, tip_vector)
  }
  if(to_name)
    tree$tip.label <- taxid2name(tree$tip.label) %>% substr(1, 30) 
  N <- length(dat$stat)

  lapply(seq(0, N, by=n), function(i){
    indices <- i:min(i+n-1, N)
    plot_one_obo_tree(tree, do.call(what=rbind, dat$stat[indices]))
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
