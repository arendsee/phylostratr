require(ggplot2, quietly=TRUE)
require(reshape2, quietly=TRUE)
require(grid, quietly=TRUE)
require(scales, quietly=TRUE)

source('phylogenetic_functions.R')

# =========
# CONSTANTS
# =========

SCORE_SCHEMES <- function(){ 
    c('mpscore', 'mmscore', 'upscore', 'umscore')
}

ID_VARS <- function(){
    c('qgb', 'hspecies', 'mrca', 'ps')
}

# ==================
# PLOTTING FUNCTIONS
# ==================

# All plotting functions require a dataframe of identical format as 
# their first argument.
# It must have the columns:
# qgb|qlocus|species|phylostratum|mrca|sciname|qlen|palen|malen|salen|pscore|mscore|sscore 

# Plot one-by-one
# plotobo <- function(d, title_='obo.pdf', lower.bound=NULL, upper.bound=NULL){
#     d <- Truncate(d, lower.bound, upper.bound)
#     d.loci <- split(d, d$qgb)
#     xlines <- as.numeric(summary(d.loci[[1]]$mrca))
#     xlines <- Reduce(sum, xlines, accumulate=TRUE)
#     xlines <- xlines[1:(length(xlines)-1)]
#     N <- length(d.loci)
#     nw = 5
#     nh = 5
#     choose.colors <- function(x){
#         colors <- ifelse(x <= 100, "green", "blue")
#         colors <- ifelse(x <= 50, "orange", colors)
#         colors <- ifelse(x <= 40, "red", colors)
#         colors <- ifelse(x <= 30, "dark red", colors)
#     }
#     pdf(title_)
#     for(page.num in 0:((N - 1) %/% (nw * nh))){
#         par(mfrow=c(nw,nh), mar=c(0,0,1,0))
#         for(j in 1:(nw*nh)){
#             index <- nw * nh * page.num + j
#             if(index > N) next
#             d.locus <- d.loci[[index]]
#             uscores <- d.locus$upscore
#             mscores <- d.locus$mpscore
#             title_ <- d.locus$qgb[1]
#             # Unmasked socres
#             plot(1:nrow(d.locus), uscores, col=choose.colors(uscores), pch=".",
#                  tck=0, xaxt="n", yaxt="n", main=title_)
#             points(1:nrow(d.locus), uscores, col=alpha("blue", 0.1), type="l")
#             # Masked scores
#             points(1:nrow(d.locus), mscores, col=choose.colors(mscores), pch=".")
#             points(1:nrow(d.locus), mscores, col=alpha("green", 0.1), type="l")
#             # Phylostrata delimiting lines
#             abline(v=(xlines - 0.5), col=rgb(0, 0, 0, 75, maxColorValue=255))
#         }
#     }
#     dev.off()
# }

plotobo <- function(d, title_='obo.pdf', lower.bound=NULL, upper.bound=NULL){
    d <- Truncate(d, lower.bound, upper.bound)
    nloci <- length(unique(d$qgb))
    nspecies <- length(unique(d$hspecies))
    # a.masked <- acast(d, qgb ~ hspecies, value.var='mmscore')
    # a.unmasked <- acast(d, qgb ~ hspecies, value.var='umscore')
    # The above used to work (maybe my current thing broke it)
    a.masked <- acast(d, qgb ~ hspecies ~ variable)[,,'mmscore']
    a.unmasked <- acast(d, qgb ~ hspecies ~ variable)[,,'umscore']
    
    choose.colors <- function(x){
        colors <- ifelse(x <= 100, "green", "blue")
        colors <- ifelse(x <= 50, "orange", colors)
        colors <- ifelse(x <= 40, "red", colors)
        colors <- ifelse(x <= 30, "dark red", colors)
    }
    xlines <- as.numeric(summary(d[1:length(unique(d$hspecies)),]$mrca))
    xlines <- Reduce(sum, xlines, accumulate=TRUE)
    xlines <- xlines[1:(length(xlines)-1)]
    nw = 5
    nh = 5
    pdf(title_)
    for(page.num in 0:((nloci - 1) %/% (nw * nh))){
        par(mfrow=c(nw,nh), mar=c(0,0,1,0))
        for(j in 1:(nw*nh)){
            index <- nw * nh * page.num + j
            if(index > nloci) next
            uscores <- a.unmasked[index, ]
            mscores <- a.masked[index, ]
            gb <- rownames(a.unmasked)[index]
            # Unmasked socres
            plot(1:nspecies, uscores, col=choose.colors(uscores), pch=".",
                 tck=0, xaxt="n", yaxt="n", main=gb)
            points(1:nspecies, uscores, col=alpha("blue", 0.1), type="l")
            # Masked scores
            points(1:nspecies, mscores, col=choose.colors(mscores), pch=".")
            points(1:nspecies, mscores, col=alpha("green", 0.1), type="l")
            # Phylostrata delimiting lines
            abline(v=(xlines + 0.5), col=rgb(0, 0, 0, 75, maxColorValue=255))
        }
    }
    dev.off()
}

SpeciesScoreBoxplot <- function(d, drop.zeros=FALSE, schemes=NULL, ...){
    d <- melt(d, id.vars=ID_VARS())
    title_ <- 'Score Comparison by Species'
    if(drop.zeros == TRUE){
        d <- d[which(d$value != 0), ]
        title_ <- paste(title_, '(zeros ignored)')
    }
    else{
        d$value <- ifelse(d$value == 0, 20, d$value)
        title_ <- paste(title_, '(zeros included)')
    }
    if(! is.null(schemes)){
        d <- d[which(d$variable %in% schemes), ]
    }
    N <- length(unique(d$mrca))
    mycol <- rep(c("#CC6666", "#9999CC"), N)[1:N]
    g <- ggplot(d, aes(x=hspecies, y=log10(value), group=hspecies)) + 
        coord_trans(ytrans='log10') +
        geom_boxplot(aes(color=mrca, outlier.color=mrca), 
                     outlier.size=1, outlier.shape=46) +
        scale_color_manual(values=mycol) +
        ylab('bitscore (log10)') +
        ggtitle(title_) + 
        theme(
            axis.text.x = element_text(angle=330, hjust=0),
            axis.title.x = element_blank(),
            legend.position="none",
            plot.margin=unit(c(0.5,2,0.5,0.5), "cm")
        ) +
        geom_hline(aes(yintercept=log10(20), alpha=0.1)) +
        facet_grid(variable ~ ., scales="free_y")
    return(g)
}


ScoreComp <- function(d, title_="None", ...){
    melted.strat <- ToStratcount(d, ...)
    g <- ggplot(melted.strat, aes(stratum, count, group=cutoff, colour=cutoff)) + 
        geom_line() +
        labs(title=title_) +
        theme(
            plot.margin=unit(c(0.5,2,0.5,0.5), "cm"), 
            legend.position='none',
            axis.text.x = element_text(angle=330, hjust=0)
             )
    return(g)
}

MaskingComp <- function(d, title_="MaskingComp", ...){
    cutoffs=c(40, 60, 80, 100)
    masked.strat <- ToStratcount(d, cutoffs=cutoffs, column='mpscore')
    unmasked.strat <- ToStratcount(d, cutoffs=cutoffs, column='upscore')
    N <- nrow(masked.strat)
    stratdat <- rbind(masked.strat, unmasked.strat)
    stratdat$mask <- c(rep('Seg', N), rep('None', N))
    g <- ggplot(stratdat, aes(x=stratum, y=count, group=mask, color=mask)) +
        geom_line() +
        ggtitle(title_) +
        theme(axis.text.x = element_text(angle=330, hjust=0)) +
        facet_grid(. ~ cutoff)
    return(g)
}

# Plots all bestscores for all input sequences
plotall <- function(d, title_='Untitled', lower.bound=NULL, upper.bound=NULL,
                    is.molten=FALSE){
    if(is.molten == FALSE){
        d <- Truncate(d, lower.bound, upper.bound)
        d <- melt(d, id.vars=ID_VARS())
    }
    N <- length(unique(d$mrca))
    mycol <- rep(c("#CC6666", "#9999CC"), N)[1:N]
    g <- ggplot(d) +
         geom_line(aes(x=hspecies, y=value, group=qgb, colour=mrca),
                   alpha=I(0.2)) +
         scale_color_manual(values=mycol) +
         theme(
            axis.text.x = element_text(angle=330, hjust=0),
            axis.title.x = element_blank(),
            legend.position="none",
            plot.margin=unit(c(0.5,2,0.5,0.5), "cm")
            ) +
         facet_grid(variable ~ .) +
         labs(y='bitscore', title=title_)
    return(g)
}
# source('io.R')
# masked <- '~/ohome/data/blast_csvs/old/masked_at.csv'
# unmasked <- '~/ohome/data/blast_csvs/old/unmasked_at.csv'
# d <- PrepareDataframe(masked=masked, unmasked=unmasked)
# dr <- BuildConstricted(d, strat=c('malvids'))
# dr <- dr[,-c(4,5)]
# dr <- dr[-which(dr$hspecies %in% unique(dr$hspecies)[30:36]), ]
# colnames(dr) <- c('ref', 'species', 'stratum', 'score')
# drt <- dr
# drt$score <- ifelse(drt$score > 100, 100, drt$score)
# pic <- function(g., basename., dir.='~/ohome/presentations/bcblab-ggplot2/', ext='jpeg', ...){
#     f <- paste0(dir., basename., '.', ext)
#     ggsave(f, g., ...)
# }
# g1 <- ggplot() +
#      geom_line(
#         aes(x=species, y=score, group=ref),
#         data=dr,
#         alpha=I(0.2))
# 
# g2 <- ggplot() +
#      geom_line(
#         aes(x=species, y=score, group=ref, color=stratum),
#         data=drt,
#         alpha=I(0.2))
#     scale_y_continuous(
#         limits=c(0,100)
#     )
# 
# g3 <- g2 +
#     theme(
#         axis.text.x = element_text(angle=330, hjust=0),
#         axis.title.x = element_blank(),
#         legend.position="none",
#         plot.margin=unit(c(0.5,2.5,0.5,0.5), "cm")
#     )
# 
# g4 <- g2 +
#     theme(
#         axis.text.x = element_text(angle=330, hjust=0),
#         axis.title.x = element_blank(),
#         legend.position="none",
#         panel.grid=element_blank(),
#         panel.background=element_rect(linetype=3, fill='blue',color='orange', size=3),
#         plot.margin=unit(c(0.5,2.5,0.5,0.5), "cm"),
#         plot.background=element_rect(linetype=2, color='red', fill='pink', size=5)
#     )
# pic(g1, 'e3p1')
# pic(g2, 'e3p2')
# pic(g3, 'e3p3')
# pic(g4, 'e3p4')


# TODO reduce memory requirements
RevenantPlot <- function(d, strata, cutoff=40, title_=NULL, strange=TRUE,
                             lower.bound=NULL, upper.bound=NULL){
    # Prepare title
    nseqs <- length(unique(d$qgb))
    if(is.null(title_)){
        cutoff.str <- paste0("bitscore cutoff=", cutoff)
        nseqs.str <- paste0("total sequences=", nseqs)
        bounds.str <- ifelse(is.null(lower.bound) & is.null(upper.bound),
                             '',
                             paste0('bounds=[', lower.bound, ',', upper.bound, ']'))
        title_ <- paste0('Constriction Plot (',
                         paste(cutoff.str, nseqs.str, bounds.str, sep='; '),
                         ') ',
                         paste(strata, collapse=', '))
    }

    # Prepare ps label colors
    species <- unique(d$hspecies)
    nspecies <- length(species)
    mycol <- rep('black', nspecies)
    stratum.species <- unique(d[which(d$mrca %in% strata), 'hspecies'])
    mycol[which(species %in% stratum.species)] = 'red'

    # Get data as list
    l <- Constrict(d=d, strata=strata, cutoff=cutoff,
                   lower.bound=lower.bound, upper.bound=upper.bound,
                   strange=strange) 

    # # Choose dataset:
    # # strange.dat has only seqs constricted at strata AND having significant
    # # matches in lower strata
    if(strange) {
        dat <- l$strange
    } else {
        dat <- l$cons
    }
    dat <- Truncate(dat, is.molten=TRUE,
                    lower.bound=lower.bound, upper.bound=upper.bound)

    # Set labels to upper right corner
    l$lab$x <- rep(nspecies, nrow(l$lab))
    l$lab$y <- rep(max(dat$value), nrow(l$lab))

    # Modify plotall output graph
    eval(substitute(
        g <- plotall(dat, title_=title_, is.molten=TRUE,
                     upper.bound=upper.bound, lower.bound=lower.bound) +
            geom_hline(aes(yintercept=cutoff, alpha=0.3)) +
            theme(
                  axis.text.x=element_text(color=mycol),
                  plot.title=element_text(hjust=0)
                 ) +
            geom_text(aes(x, y, label=labs, vjust=1, hjust=1, alpha=0.3),
                      data=l$lab),
    env=list(cutoff=cutoff, mycol=mycol)))
    return(g)
}

# masked <- '~/ohome/projects/ghosts/csvs/by_gb/masked_hs_100-bit-ghosts.csv'
# unmasked <- '~/ohome/projects/ghosts/csvs/by_gb/unmasked_hs_100-bit-ghosts.csv'
# d <- PrepareDataframe(masked=masked, unmasked=unmasked)
# g <- RevenantPlot(d, 'Amniota', cutoff=40, strange=TRUE, upper.bound=150, lower.bound=30)

NoisePlot <- function(d, zero.strata, noise.strata, cutoffs=c(40, 60, 80, 100),
                      schemes=c('umscore', 'mmscore'), title_=NULL, cons=NULL,
                      binwidth=NULL){
    if(is.null(title_)){
        title_ <- paste0('Noise Plot (Bound on ',
                         paste(zero.strata, collapse=','),
                         ') ',
                         paste(noise.strata, collapse=','))
    }
    if(is.null(cons)){
        cons <- MakeNoise(d, zero.strata, noise.strata, cutoffs)
    }
    # cons$value <- log2(cons$value)
    cons <- cons[which(cons$variable %in% schemes), ]
    if(is.null(binwidth)){
        binwidth <- log2(max(cons$value) - min(cons$value)) / 240
    }
    g <- ggplot(cons, aes(x=value)) +
        geom_histogram(binwidth=1/8) +
        scale_x_continuous(
            trans='log2',
            breaks=trans_breaks('log2', function(x) round(2^x))
        ) +
        ggtitle(title_) +
        # theme(
        #     axis.text.y=element_blank(),
        #     axis.ticks.y=element_blank(),
        #     axis.title.y=element_blank()
        # ) +
        xlab('bitscore') +
        facet_grid(. ~ cutoff ~ variable ~ hspecies)
    return(g)
}

# =================
# UTILITY FUNCTIONS
# =================

MakeStrangeLabels <- function(cons, strange){
    s <- SCORE_SCHEMES()
    k <- length(s) 
    labs <- rep(0, k)
    for (i in 1:k){
        # Number of total constricted sequences
        nc <- length(unique(cons[which(cons$variable == s[i]), ]$qgb))
        # Number of constricted sequences that have later significant matches
        ns <- length(unique(strange[which(strange$variable == s[i]), ]$qgb))
        labs[i] <- paste0(ns, '|', nc)
    }
    baddat <- data.frame(x=rep(0, k),
                         y=rep(0, k),
                         variable=s,
                         labs=labs)
    return(baddat)
}

# Truncates all scores
Truncate <- function(d, lower.bound=NULL, upper.bound=NULL, is.molten=FALSE){
    s <- SCORE_SCHEMES()
    if(! is.null(lower.bound)){
        if(is.molten){
            d$value <- ifelse(d$value < lower.bound, lower.bound, d$value)
        }
        else{
            d[, s] <- apply(d[, s], c(1,2), max, lower.bound)
        }
    }
    if(! is.null(upper.bound)){
        if(is.molten){
            d$value <- ifelse(d$value > upper.bound, upper.bound, d$value)
        }
        else{
            d[, s] <- apply(d[, s], c(1,2), min, upper.bound)
        }
    }
    return(d)
}

# masked <- '~/ohome/data/blast_csvs/masked_at.csv'
# unmasked <- '~/ohome/data/blast_csvs/unmasked_at.csv'
# d <- PrepareDataframe(masked, unmasked) 
# g <- ScoreComp(d, title_="Arabisopsis thaliana phylostratigraph",
#                cutoffs=100)
# ggsave('100bit_phylostratigraph.jpeg')
# 
