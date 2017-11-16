eval_bins <- function(evalue){
  bins <- c("e >= 1", "e < 1", "e < 1e-3", "e < 1e-5", "e < 1e-20")
  bin <- ifelse(evalue >= 1e-20, bins[4], bins[5])
  bin <- ifelse(evalue >= 1e-5,  bins[3], bin)
  bin <- ifelse(evalue >= 1e-3,  bins[2], bin)
  bin <- ifelse(evalue >= 1,     bins[1], bin)
  bin <- factor(bin, levels=bins)
  bin
}

plot_one_obo <- function(stats, qseqid, xlines){
  stats$index = 1:nrow(stats)
  ggplot2::ggplot(stats) +
    ggplot2::geom_point(ggplot2::aes(x=.data$index, y=.data$score, color=eval_bins), size=0.5) +
    ggplot2::scale_color_manual(values=c('darkred', 'red', 'orange', 'green', 'blue')) +
    ggplot2::geom_vline(data=xlines, ggplot2::aes(xintercept=xlines), alpha=0.4) +
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

plot_obo <- function(d, lower.bound=0, upper.bound=100){

  d <- dplyr::arrange(d, -.data$ps, .data$staxid)

  if(!is.null(upper.bound))
    d$score <- ifelse(d$score > upper.bound, upper.bound, d$score)
  if(!is.null(lower.bound))
    d$score <- ifelse(d$score < lower.bound, lower.bound, d$score)

  nspecies <- length(unique(d$staxid))

  taxidmap <- d %>%
    dplyr::select(.data$staxid, .data$ps) %>%
    dplyr::distinct() %>%
    dplyr::arrange(-.data$ps)

  # It is essential to turn the taxon IDs into ordered factors, otherwise they
  # will be interpreted as numbers, and plotted accordingly.
  d$staxid <- factor(as.character(d$staxid), levels=as.character(taxidmap$staxid))

  xlines <- taxidmap %>%
    dplyr::group_by(.data$ps) %>%
    dplyr::summarize(n = dplyr::n()) %>%
    dplyr::arrange(-.data$ps) %>%
    { Reduce(sum, .$n, accumulate=TRUE) } %>%
    { .[-length(.)] } %>%
    { . + 0.5 } %>%
    as.data.frame

  stats <- d %>%
    dplyr::select(.data$qseqid, .data$score, .data$evalue, .data$staxid) %>%
    merge(taxidmap, by='staxid') %>%
    dplyr::group_by(.data$qseqid) %>%
    dplyr::mutate(index=1:dplyr::n()) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(eval_bins=eval_bins(.data$evalue)) %>%
    {
      base::split(., f=factor(.$qseqid))
    } %>%
    lapply(dplyr::arrange, -.data$ps, .data$staxid)

  lapply(
    names(stats), 
    function(qseqid) {
      plot_one_obo(stats[[qseqid]], qseqid, xlines)
    }
  )
}

make_obo_pdf <- function(d, title='obo.pdf', width=1, height=5, ...){
  plots <- plot_obo(d, ...)
  nloci <- length(unique(d$qseqid))
  pdf(title)
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

plot_constriction <- function(d){
  # FIXME: STUB
}

plot_noise <- function(){
  # FIXME: STUB
}
