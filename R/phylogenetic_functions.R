Constrict <- function(d, strange=TRUE, labs=TRUE, ...){
    if(labs & ! strange){
        strange = TRUE
    }
    out <- list()
    out$cons <- BuildConstricted(d, ...)
    if(strange){
        out$strange <- BuildStrange(out$cons, ...)
    }
    if(labs){
        out$lab <- MakeStrangeLabels(out$cons, out$strange)
    }
    return(out)
}

BuildConstricted <- function(d, strata, cutoff=40, schemes='mmscore'){
    cons <- EmptyDataframe(is.molten=TRUE)
    for(column in schemes){
        scores <- d[which(d$mrca %in% strata), c('qgb', column)]
        has.homo <- unique(scores$qgb[which(scores[,2] > cutoff)])
        no.homo <- setdiff(unique(d$qgb), has.homo)
        dd <- d[which(d$qgb %in% no.homo), c(ID_VARS(), column)]
        dd.melt <- melt(dd, id.vars=ID_VARS())
        cons <- rbind(cons, dd.melt)
    }
    return(cons)
}

BuildStrange <- function(cons, strata, cutoff=40, schemes='mmscore'){
    strange <- EmptyDataframe(is.molten=TRUE)
    for (column in schemes){
        ps.level <- min(cons[which(cons$mrca %in% strata), 'ps']) 
        s <- cons[which(cons$variable == column),]
        ss <- s[which(s$ps < ps.level), ]
        strange.qgb <- unique(ss[which(ss$value > cutoff), 'qgb'])
        s <- s[which(s$qgb %in% strange.qgb), ]
        strange <- rbind(strange, s)
    }
    return(strange)
}

# Return all sequences which resurrect after a double constriction
# OUTPUT
#       vector of qgb revenants
FindRevenants <- function(d, strata=NULL, ...){
    if(is.null(strata)){
        strata <- unique(d$mrca)
        strata <- strata[-c(1, length(strata))]
    }
    revenants <- c()
    for(i in 1:(length(strata)-1)){
        for(j in (i+1):length(strata)){
            cons <- Constrict(d, strata=c(strata[i], strata[j]), strange=TRUE, labs=FALSE, ...)
            if(nrow(cons$strange) == 0) { next }
            revenants <- union(revenants, unique(cons$strange$qgb))
        }
    }
    return(revenants)
}

# INPUT: molten dataframe
#       hspecies|mrca|$score|<<other columns (ignored)>>
# OUTPUT: molten dataframe
#       cutoff|stratum|count
ToStratcount <- function(d, cutoffs=c(5 * (8:20)), score='umscore'){
    if(score == 'mmeval' | score == 'umeval'){
        strata <- ToStratCountEvalue(d, cutoffs, score)
    } else {
        comp <- function(x, cutoff) {max(1, which(x > cutoff))} 
        strata <- ToStratCountBase(d, cutoffs, comp, score)
    }
    return(strata)
}

ToStratCountEvalue <- function(d, cutoffs=c(1e-4, 1e-5, 1e-6), fill='mmeval'){
    comp <- function(x, cutoff) {max(1, which(x < cutoff))} 
    strata <- ToStratCountBase(d, cutoffs, comp, fill)
    return(strata)
}

ToStratCountBase <- function(d, cutoffs, comp, fill){
    strat.names <- d[1:length(unique(d$hspecies)), c('hspecies', 'mrca')]
    out <- data.frame(stratum=character, cutoff=numeric, count=numeric)
    a <- acast(d, qgb ~ hspecies ~ variable)[,,fill]
    for(cutoff in cutoffs){
        last <- apply(a, 1, function(x) comp(x, cutoff))
        strat <- strat.names[last, 2] 
        count <- as.numeric(summary(strat))
        strat.dat <- data.frame(stratum=levels(strat),
                                cutoff=rep(cutoff, length(count)),
                                count=count)
        out <- rbind(out, strat.dat)
    }
    out$stratum <- factor(out$stratum, levels=levels(strat.names$mrca))
    return(out)
}

GetLSG.by.evalue <- function(...){
    comp <- function(x, cutoff) {max(1, which(x < cutoff))} 
    return(GetLSG(comp=comp, fill='mevalue', ...))
}
GetLSG.by.score <- function(...){
    comp <- function(x, cutoff) {max(1, which(x > cutoff))} 
    return(GetLSG(comp=comp, ...))
}
# Get lineage specific genes
# INPUT: melted data
GetLSG <- function(d., cutoff, comp, fill, sciname=FALSE){
    strat.names <- d.[1:length(unique(d.$hspecies)), c('hspecies', 'mrca', 'mrca_sciname')]
    a <- acast(d., qgb ~ hspecies ~ variable)[unique(d.$qgb),strat.names$hspecies,fill]
    # Index of most distant match
    last <- apply(a, 1, function(x) comp(x, cutoff))
    mrca.col = if(sciname) 3 else 2
    strat <- data.frame(gb=attributes(last), stratum=strat.names[last, mrca.col]) 
    return(strat)
}

MakeNoise <- function(d, zero.strata, noise.strata, cutoffs=c(40, 60, 80, 100)){
    cons <- EmptyDataframe(is.molten=TRUE)
    ps.levels <- d[which(d$mrca %in% noise.strata), 'ps'] 
    for (cutoff in cutoffs){
        cc <- BuildConstricted(d, zero.strata, cutoff)
        cc <- cc[which(cc$ps %in% ps.levels), ]
        cc <- cc[which(cc$value != 0), ]
        cc$cutoff <- rep(cutoff, nrow(cc))
        cons <- rbind(cons, cc)
    }
    return(cons)
}
