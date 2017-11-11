plot_obo <- function(d, title='obo.pdf', lower.bound=NULL, upper.bound=NULL){

  d <- dplyr::arrange(d, -ps, staxid)

  if(!is.null(upper.bound))
    d$score <- ifelse(d$score > upper.bound, upper.bound, d$score)
  if(!is.null(lower.bound))
    d$score <- ifelse(d$score < lower.bound, lower.bound, d$score)

  nloci <- length(unique(d$qseqid))

  nspecies <- length(unique(d$staxid))

  # blue      |  e < 1e-20  |  safe
  # green     |  e < 1e-5   |  under the default threshold
  # orange    |  e < 0.001  |  maybe something?
  # red       |  e < 1      |  nothing
  # dark red  |  e > 1      |  really nothing
  choose.colors <- function(evalue){
    colors <- ifelse(evalue >= 1e-20, "green",    "blue")
    colors <- ifelse(evalue >= 1e-5,  "orange",   colors)
    colors <- ifelse(evalue >= 1e-3,  "red",      colors)
    colors <- ifelse(evalue >= 1,     "dark red", colors)
    colors
  }

  taxidmap <- d %>%
    dplyr::select(staxid, ps) %>%
    dplyr::distinct() %>%
    dplyr::arrange(-ps)

  xlines <- taxidmap %>%
    dplyr::group_by(ps) %>%
    dplyr::summarize(n = n()) %>%
    dplyr::arrange(-ps) %>%
    { Reduce(sum, .$n, accumulate=TRUE) } %>%
    { .[-length(.)] }

  nw = 1
  nh = 5

  qseqids <- unique(d$qseqid)

  pdf(title_)
  for(page.num in 0:((nloci - 1) %/% (nw * nh))){
    par(mfrow=c(nw,nh), mar=c(0,0,1,0))
    for(j in 1:(nw*nh)){
      index <- nw * nh * page.num + j
      if(index > nloci) next
      qseqid <- qseqids[index]

      staxid_scores <- d[which(d$qseqid == qseqid), ] %>%
        dplyr::select(score, evalue, staxid) %>%
        merge(taxidmap, by='staxid')

      plot(
        1:nspecies,
        scores,
        col  = choose.colors(scores),
        tck  = 0,
        xaxt = "n",
        yaxt = "n",
        main = qseqid
      )
      # Phylostrata delimiting lines
      abline(v=(xlines + 0.5), col=rgb(0, 0, 0, 75, maxColorValue=255))

    }
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
