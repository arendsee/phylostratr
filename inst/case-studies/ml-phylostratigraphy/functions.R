setClass(
  "CountMatrix",
  representation(
    x = "matrix",
    xlab = "character",
    ylab = "character"
  )
)

countMatrix <- function(x, ylab="from", xlab="to"){
  m <- new("CountMatrix") 
  m@x = x
  m@xlab = xlab
  m@ylab = ylab
  m
}

plot.CountMatrix <- function(x, y=NULL,
  labels      = names(m)[1:2],
  value_trans = identity,
  normalize   = normalize_matrix_by_row,
  scheme      = ggplot2::scale_fill_gradient(low = "grey", high = "red"),
  ...
){

  cnt <- x@x
  m <- normalize(cnt)

  cnt <- reshape2::melt(cnt) %>% dplyr::filter(value > 0)
  m <- reshape2::melt(m)

  names(m)[1:3] <- c("a", "b", "value")
  names(cnt)[1:3] <- c("a", "b", "n")

  m$value <- value_trans(m$value)
  cnt$value <- value_trans(cnt$value)

  ggplot2::ggplot() +
    ggplot2::geom_bin2d(data=m, ggplot2::aes(x=b, y=a, fill=value)) +
    ggplot2::xlab(x@xlab) +
    ggplot2::ylab(x@ylab) +
    ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle=270, hjust=0, vjust=1),
        legend.title = ggplot2::element_blank()
    ) +
    scheme +
    ggplot2::geom_text(data=cnt, mapping=ggplot2::aes(x=b, y=a, label=n), size=2)
}


run_comparison <- function(results){
  list(
    holm05 = stratify(results, classify_assuming_iid( 5e-2, method='holm' ))[,       c('qseqid', 'mrca_name')],
    bon05  = stratify(results, classify_assuming_iid( 5e-2, method='bonferroni' ))[, c('qseqid', 'mrca_name')],
    bon5   = stratify(results, classify_assuming_iid( 1e-5, method='bonferroni' ))[, c('qseqid', 'mrca_name')]
  )
}

access_cache <- function(cache, FUN, ...){
  if(file.exists(cache)){
    d <- readRDS(cache)
  } else {
    d <- run_comparison(results) 
    saveRDS(d, cache)
  }
  d
}

make_matrix_from_two_strata <- function(d){
  # for now, only use the first two columns
  d <- d[, c(1,2)]
  labels <- names(d)
  names(d) <- c("a", "b")

  mat <- dplyr::group_by(d, a, b) %>%
    dplyr::count() %>%
    reshape2::acast(a ~ b, fill=0)

  countMatrix(x=mat, ylab=labels[1], xlab=labels[2])
}

# Build a matrix relating 2014 study to 2018, with 2014 classifications on the
# y axis and 2018 on the x axis.
make_TIPS_comparison_matrix <- function(strata, results, classifier){
  # Get the strata from the current analysis
  a2018 <- stratify(results)
  # Get and clean the (Arendsee et al., 2014) data (stored in phylostratr package)
  a2014 <- readr::read_tsv(system.file('extdata', 'arendsee2014_strata.tab', package='phylostratr'))
  # Add NCBI taxonomy ID to tips data
  a2014$mrca <- taxizedb::name2taxid(a2014$name)
  # match column naming conventions
  a2014 <- dplyr::select(a2014, qseqid, mrca, ps, mrca_name=name)
  # underscores to spaces 
  a2014$mrca_name <- sub('_', ' ', a2014$mrca_name)
  # get a lineage map for Arabidopsis
  strata_map <- taxizedb::classification('3702')[[1]]
  # factor the backbone
  a2014$mrca_name <- droplevels(factor(a2014$mrca_name, levels=strata_map$name))

  ps <- standardize_strata(list(
      a2014 = a2014,
      a2018 = a2018
  ))
  ps <- lapply(ps, function(x) x[c('qseqid', 'mrca_name')])

  species2018 <- strata %>% strata_fold(function(s) taxizedb::taxid2name(s@tree$tip.label)) %>%
    { names(.) <- taxizedb::taxid2name(names(.)); . } %>%
    reshape2::melt() %>%
    dplyr::select(phylostratum = L1, species = value)

  species2014 <- readr::read_tsv(system.file('case-studies', 'ml-phylostratigraphy', 'tips-ages.tab', package='phylostratr'))

  d <- merge(ps[[1]], ps[[2]], by="qseqid", suffixes=paste0(".", names(ps)))

  make_matrix_from_two_strata(d[, -1])
}

make_significance_matrices <- function(results, cache='significance_list.Rda'){
  d <- access_cache(cache, run_comparison, results)

  qseqids <- lapply(d, function(x) x$qseqid) %>%
    {Reduce(f=intersect, .[-1], .[[1]])}
  d <- lapply(d, function(x) x[x$qseqid %in% qseqids, ]) %>%
    lapply(dplyr::arrange, qseqid)
  m <- as.data.frame(lapply(d, function(x) x$mrca_name))

  gs <- list()
  for(i in 1:(ncol(m) - 1)){
    gss <- list()
    for(j in (i+1):ncol(m)){
      gss[[names(m)[j]]] <- make_matrix_from_two_strata(m[, c(i, j)])
    }
    gs[[names(m)[i]]] <- gss
  }

  gs
}


# Create matrix counting the distances between factor levels across groups
make_jump_matrix <- function(d, labels=NULL){
  # expect two columns
  stopifnot(ncol(d) == 2)
  names(d) <- c("group", "index")
  stopifnot(is.numeric(d$index))

  mat <- d %>%
    dplyr::select(group, index) %>%
    dplyr::distinct() %>%
    dplyr::arrange(index) %>%
    dplyr::group_by(group) %>%
    dplyr::transmute(
      from = c(NA, head(index, -1)),
      to   = c(NA, tail(index, -1))
    ) %>%
    dplyr::ungroup() %>%
    dplyr::filter(!is.na(from) & !is.na(to)) %>%
    dplyr::select(from, to) %>%
    dplyr::group_by(from, to) %>%
    dplyr::count() %>%
    dplyr::ungroup() %>%
    reshape2::acast(from ~ to, fill=0) %>% t

  if(!is.null(labels)){
    dimnames(mat) <- list(labels[as.integer(dimnames(mat)[[1]])], labels[as.integer(dimnames(mat)[[2]])])
  }

  countMatrix(x=mat)
}

# dataframe_significance_comparison <- function(results, cache='significance_list.Rda'){
#   d <- access_cache(cache, run_comparison, results)
#   m <- lapply(d, function(x) summary(x$mrca_name)) %>%
#     as.data.frame()
# }

# plot_skipped_distance <- function(results, is_present, strata_levels, ...){
#   results[is_present, c("qseqid", "ps")] %>%
#     make_jump_matrix(labels=strata_levels)
#
#     reshape2::melt() %>%
#     dplyr::select(from=Var1, to=Var2, n=value) %>%
#     dplyr::mutate(
#       from = factor(strata_levels[from], levels=strata_levels),
#       to   = factor(strata_levels[to],   levels=strata_levels)
#     ) %>%
#     dplyr::select(to, from, n)
#   m <- cnt %>%
#     dplyr::group_by(to) %>%
#     dplyr::mutate(n = n / sum(n))
#   cnt <- dplyr::filter(cnt, n != 0)
#   plot_heatmap_with_counts(m, cnt, ...)
# }
#
# plot_significance_comparison <- function(results, cache='significance_list.Rda', ...){
#   d <- access_cache(cache, run_comparison, results)
#   qseqids <- lapply(d, function(x) x$qseqid) %>%
#     {Reduce(f=intersect, .[-1], .[[1]])}
#   d <- lapply(d, function(x) x[x$qseqid %in% qseqids, ]) %>%
#     lapply(dplyr::arrange, qseqid)
#   m <- as.data.frame(lapply(d, function(x) x$mrca_name))
#   gs <- list()
#   for(i in 1:(ncol(m) - 1)){
#     for(j in (i+1):ncol(m)){
#       gs[[paste0(names(m)[i], "-", names(m)[j])]] <- plot_compare_two_strata(m[, c(i, j)], ...)
#     }
#   }
#   gs
# }
#
#
# plot_significance_lines <- function(results, cache='significance_list.Rda'){
#   d <- access_cache(cache, run_comparison, results)
#   m <- apply(m, 2, function(x) log2(x / m[, 1])) %>% as.data.frame
#   m$phylostratum <- factor(rownames(m), levels=rownames(m))
#
#   ggplot(melt(m)) +
#     geom_path(aes(x = phylostratum, y=value, color=variable, group=variable)) +
#     scale_color_brewer(palette="Paired") +
#     theme(
#         axis.text.x = element_text(angle=270, hjust=0, vjust=1)
#     )
# }
#
# plot_skipped_counts <- function(revenants, strata_table, ...){
#   strata_levels <- levels(strata_table$mrca_name)
#   m <- revenants[, c('basal_ps', 'ps')]
#   m$ps       <- factor(strata_levels[m$ps],       levels=strata_levels)
#   m$basal_ps <- factor(strata_levels[m$basal_ps], levels=strata_levels)
#   names(m) <- c("skipped_ps", "assigned_ps")
#   normalize_with = function(d){
#     t(apply(d, 1, function(x) x / summary(strata_table$mrca_name)))
#   }
#   plot_compare_two_strata(m, normalize_with=normalize_with, ...)
# }
#
#
# # Find the number of hits of genes assigned to each phylostratum against each
# # phylostratum
# plot_hits_by_phylostratum <- function(results, is_present, strata_table, ...){
#   m <- results[, c('qseqid', 'ps')]
#   m$present = is_present
#
#   strata_levels <- levels(strata_table$mrca_name)
#
#   ### make a table with the columns: qseqid | ps | present | ps_class
#   ### Where ps_class is the oldest stratum with an inferred homolog
#   focal_species = length(strata_levels)
#   mm <- m %>%
#     dplyr::distinct() %>%
#     dplyr::group_by(qseqid, ps) %>%
#     dplyr::summarize(present = any(present)) %>%
#     dplyr::group_by(qseqid) %>%
#     dplyr::mutate(ps_class = ifelse(any(present), min(ps[present]), focal_species))
#   mm$ps       <- factor(strata_levels[mm$ps],       levels=strata_levels)
#   mm$ps_class <- factor(strata_levels[mm$ps_class], levels=strata_levels)
#
#   cnt <- dplyr::select(mm, ps_class, ps, present) %>%
#     dplyr::group_by(ps_class, ps) %>%
#     dplyr::summarize(n = sum(present))
#
#   strata_counts <- as.vector(table(strata_table$mrca_name))
#
#   normalize_with = function(d){
#     apply(d, 2, function(x) x / strata_counts)
#   }
#
#   cnt %>%
#     reshape2::acast(ps_class ~ ps, fill=0) %>%
#     normalize_with  %>%
#     reshape2::melt(value.name='value') %>%
#     plot_heatmap_with_counts(cnt, labels=names(cnt)[1:2], ...) +
#       ggplot2::labs(
#         x = "Phylostratrum",
#         y = "Assigned phylostratum",
#         title = "Number of focal genes with inferred homologs"
#       )
# }
#
# # like above but genewise
# plot_hits_by_gene <- function(results, is_present, strata_table, ...){
#
#   strata_levels <- levels(strata_table$mrca_name)
#   focal_species <- length(strata_levels)
#
#   mm <- results[, c('qseqid', 'staxid', 'ps')] %>%
#     dplyr::mutate(present = is_present) %>%
#     dplyr::distinct() %>%
#     dplyr::group_by(qseqid) %>%
#     dplyr::mutate(ps_class = ifelse(any(present), min(ps[present]), focal_species)) %>%
#     dplyr::arrange(ps, staxid)
#
#   taxid_map  <- taxizedb::taxid2name(unique(mm$staxid))
#   names(taxid_map) <- unique(mm$staxid)
#   mm$species <- factor(taxid_map[as.character(mm$staxid)], levels=taxid_map)
#   mm$ps_class <- factor(strata_levels[mm$ps_class], levels=strata_levels)
#
#   cnt <- dplyr::select(mm, ps_class, species, present) %>%
#     dplyr::group_by(ps_class, species) %>%
#     dplyr::summarize(n = sum(present))
#
#   strata_counts <- as.vector(table(strata_table$mrca_name))
#
#   normalize_with = function(d){
#     apply(d, 2, function(x) x / strata_counts)
#   }
#
#   cnt %>%
#     reshape2::acast(ps_class ~ species, fill=0) %>%
#     normalize_with  %>%
#     reshape2::melt(value.name='value') %>%
#     plot_heatmap_with_counts(cnt, labels=names(cnt)[1:2], ...) +
#       ggplot2::labs(
#         x = "Species",
#         y = "Assigned phylostratum",
#         title = "Number of focal genes with inferred homologs"
#       )
# }


# plot_compare_to_TIPS <- function(strata, results, classifier, ...){
#
#   plot_compare_two_strata(d[, -1], ...) +
#     ggplot2::labs(
#       x = "Manually curated proteome set (2014 TIPS)",
#       y = "Automatically selected set (phylostratr)",
#       title = "Comparison between manual and automatic data selection"
#     )
# }
