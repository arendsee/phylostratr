ghost_analysis_2013.9.26 <- function(d, species='unnamed',
                                     outfile=paste0(species, '_100-bit-ghosts.pdf')){
    # Plot window width and height
    w=20
    h=10.3
    
    # Writes to a potentially massive, multi-page pdf depicting the scores for
    # each locus
    plotobo(d, title_=paste0(species, "_obo_100-bit-ghosts.pdf"),
            lower.bound=30, upper.bound=150)

    pdf(file=outfile, width=w, height=h)
        print(plotall(d, title_=paste(species, 'Gene Scores'),
                lower.bound=30, upper.bound=150))
        print(ScoreComp(d, title_=paste(species, 'Bitscore Comparison')))
        print(MaskingComp(d, title_=paste(species, "Masking Comparison")))
        print(SpeciesScoreBoxplot(d, drop.zeros=TRUE))
        print(SpeciesScoreBoxplot(d, drop.zeros=FALSE))

        MassConstrict(d, cutoff=40,
                      lower.bound=30, upper.bound=150,
                      strange=TRUE)
    dev.off()
}


genome_analysis_2013.9.27 <- function(d, species='unnamed',
                                     outfile=paste0(species, '_full_genome.pdf')){
    # Plot window width and height
    w=20
    h=10.3
    
    # Writes to a potentiall massive, multi-page pdf depicting the scores for
    # each locus
    plotobo(d, title_=paste0(species, "_obo.pdf"),
            lower.bound=30, upper.bound=150)
    # Without a 150-bistscore upper bound
    plotobo(d, title_=paste0(species, "_obo.pdf"),
            lower.bound=30, upper.bound=NULL)


    pdf(file=outfile, width=w, height=h)
        print(SpeciesScoreBoxplot(d, drop.zeros=TRUE))
        print(SpeciesScoreBoxplot(d, drop.zeros=FALSE))
        MassConstrict(d, cutoff=40,
                      lower.bound=30, upper.bound=150,
                      strange=TRUE)
    dev.off()
}

human_genome_analysis_2013.9.30 <- function(d, species='Human',
                                     outfile=paste0(species, '_full_genome.pdf')){
    # Plot window width and height
    w=20
    h=10.3
    
    # Writes to a potentiall massive, multi-page pdf depicting the scores for
    # each locus
    plotobo(d, title_=paste0(species, "_obo.pdf"),
            lower.bound=30, upper.bound=150)
    # Without a 150-bistscore upper bound
    plotobo(d, title_=paste0(species, "_obo.pdf"),
            lower.bound=30, upper.bound=NULL)


    pdf(file=outfile, width=w, height=h)
        print(SpeciesScoreBoxplot(d, drop.zeros=TRUE, schemes=c('mmscore', 'mpscore')))
        print(SpeciesScoreBoxplot(d, drop.zeros=FALSE, schemes=c('mmscore', 'mpscore')))
        MassConstrict(d, cutoff=40,
                      lower.bound=30, upper.bound=150,
                      strange=TRUE)
    dev.off()
}

At_noise_analysis_2013.9.30 <- function(d, outfile='noise.pdf'){
    co <- c(40, 60, 80, 100, 40, 100, 250, Inf)
    s <- c('mmscore', 'umscore', 'mpscore', 'upscore')
    l1 <- 'malvids'
    l2 <- 'rosids'
    l3 <- 'core eudicotyledons'
    l4 <- 'Magnoliophyta'
    l5 <- 'Eukaryota'
    l6 <- 'cellular organisms'

    w=20
    h=10.3
    plotset <- function(d, zs, ns, co, s, cplot=FALSE){
        if(cplot){
            print(RevenantPlot(d, strata=zs, strange=TRUE, cutoff=40,
                                   lower.bound=30, upper.bound=150))
        }
        co1 <- co[1:4]
        cons <- MakeNoise(d, zs, ns, co1)
        print(NoisePlot(d, zs, ns, co1, s[1:2], cons=cons))
        print(NoisePlot(d, zs, ns, co1, s[3:4], cons=cons))
        co2 <- co[5:8]
        cons <- MakeNoise(d, zs, ns, co2)
        print(NoisePlot(d, zs, ns, co2, s[1:2], cons=cons))
        print(NoisePlot(d, zs, ns, co2, s[3:4], cons=cons))
    }
    pdf(file=outfile, width=w, height=h)
        plotset(d=d, zs=l1, ns=l2, co=co, s=s, TRUE)
        plotset(d=d, zs=l2, ns=l3, co=co, s=s, TRUE)
        plotset(d=d, zs=l3, ns=l4, co=co, s=s, TRUE)
        plotset(d=d, zs=l4, ns=l5, co=co, s=s, TRUE)
        plotset(d=d, zs=l5, ns=l6, co=co, s=s, TRUE)

        plotset(d=d, zs=c(l1, l2), ns=l3, co=co, s=s, TRUE)
        plotset(d=d, zs=c(l1, l2), ns=l4, co=co, s=s)
        plotset(d=d, zs=c(l1, l2), ns=l5, co=co, s=s)
        plotset(d=d, zs=c(l1, l2), ns=l6, co=co, s=s)

        plotset(d=d, zs=c(l2, l3), ns=l4, co=co, s=s, TRUE)
        plotset(d=d, zs=c(l2, l3), ns=l5, co=co, s=s)
        plotset(d=d, zs=c(l2, l3), ns=l6, co=co, s=s)
        
        plotset(d=d, zs=c(l3, l4), ns=l5, co=co, s=s, TRUE)
        plotset(d=d, zs=c(l3, l4), ns=l6, co=co, s=s)
    dev.off()
}


Hs_noise_analysis_2013.9.29 <- function(d, outfile='Human_noise.pdf'){
    co <- c(40, 60, 80, 100, 40, 100, 250, Inf)
    s <- c('mmscore', 'mpscore')
    l1 <- 'Homininae'
    l2 <- 'Catarrhini'
    l3 <- 'Simiiformes'
    l4 <- 'Euarchontoglires'
    l5 <- 'Eutheria'
    l6 <- 'Theria'
    l7 <- 'Amniota'
    l8 <- 'Euteleostomi'
    l9 <- 'Chordata'
    l10 <- 'Bilateria'
    l11 <- 'Opisthokonta'
    l12 <- 'Eukaryota'
    l13 <- 'cellular organisms'

    plotset <- function(d, zs, ns, co, s, cplot=FALSE){
        if(cplot){
            print(RevenantPlot(d, strata=zs, strange=TRUE, cutoff=40,
                               lower.bound=30, upper.bound=150))
        }
        cons <- MakeNoise(d, zs, ns, co[1:4])
        print(NoisePlot(d, zs, ns, co1, s, cons=cons))
        cons <- MakeNoise(d, zs, ns, co[5:8])
        print(NoisePlot(d, zs, ns, co2, s, cons=cons))
    }
    pdf(file=outfile, width=20, height=10.3)
        plotset(d=d, zs=l1, ns=l2, co=co, s=s, TRUE)
        plotset(d=d, zs=l2, ns=l3, co=co, s=s, TRUE)
        plotset(d=d, zs=l3, ns=l4, co=co, s=s, TRUE)
        plotset(d=d, zs=l4, ns=l5, co=co, s=s, TRUE)
        plotset(d=d, zs=l5, ns=l6, co=co, s=s, TRUE)
        plotset(d=d, zs=l6, ns=l7, co=co, s=s, TRUE)
        plotset(d=d, zs=l7, ns=l8, co=co, s=s, TRUE)
        plotset(d=d, zs=l8, ns=l9, co=co, s=s, TRUE)
        plotset(d=d, zs=l9, ns=l10, co=co, s=s, TRUE)
        plotset(d=d, zs=l10, ns=l11, co=co, s=s, TRUE)
        plotset(d=d, zs=l11, ns=l12, co=co, s=s, TRUE)
        plotset(d=d, zs=l12, ns=l13, co=co, s=s, TRUE)
        
        plotset(d=d, zs=c(l4, l5), ns=l6, co=co, s=s, TRUE)
        plotset(d=d, zs=c(l4, l5), ns=l7, co=co, s=s)
        plotset(d=d, zs=c(l4, l5), ns=l8, co=co, s=s)
        plotset(d=d, zs=c(l4, l5), ns=l9, co=co, s=s)
        plotset(d=d, zs=c(l4, l5), ns=l10, co=co, s=s)
        plotset(d=d, zs=c(l4, l5), ns=l11, co=co, s=s)
        plotset(d=d, zs=c(l4, l5), ns=l12, co=co, s=s)
        plotset(d=d, zs=c(l4, l5), ns=l13, co=co, s=s)

        plotset(d=d, zs=c(l6, l7), ns=l8, co=co, s=s, TRUE)
        plotset(d=d, zs=c(l6, l7), ns=l9, co=co, s=s)
        plotset(d=d, zs=c(l6, l7), ns=l10, co=co, s=s)
        plotset(d=d, zs=c(l6, l7), ns=l11, co=co, s=s)
        plotset(d=d, zs=c(l6, l7), ns=l12, co=co, s=s)
        plotset(d=d, zs=c(l6, l7), ns=l13, co=co, s=s)

        s1 <- c('Tetrapoda', 'Sarcopterygii')
        plotset(d=d, zs=s1, ns=l8, co=co, s=s, TRUE)
        plotset(d=d, zs=s1, ns=l9, co=co, s=s)
        plotset(d=d, zs=s1, ns=l10, co=co, s=s)
        plotset(d=d, zs=s1, ns=l11, co=co, s=s)
        plotset(d=d, zs=s1, ns=l12, co=co, s=s)
        plotset(d=d, zs=s1, ns=l13, co=co, s=s)
    dev.off()
}

At_Revenant_2013.10.1 <- function(d, outfile='At_revenant.pdf'){
    co <- c(35, 40, 45, 50, 55, 60, 70, 80)
    l1 <- 'Arabidopsis'
    l2 <- 'Camelineae'
    l3 <- 'Brassicaceae'
    l4 <- 'malvids'
    l5 <- 'rosids'
    l6 <- 'core eudicotyledons'
    l7 <- 'Magnoliophyta'
    l8 <- 'Spermatophyta'
    l9 <- 'Tracheophyta'
    l10 <- 'Embryophyta'
    l11 <- 'Viridiplantae'
    l12 <- 'Eukaryota'
    l13 <- 'cellular organisms'

    w=20
    h=10.3
    plotset <- function(d, s, co){
        for (cutoff in co){
            print(RevenantPlot(d, strata=s, strange=TRUE, cutoff=cutoff,
                                   lower.bound=15, upper.bound=125))
        }
    }
    pdf(file=outfile, width=w, height=h)

        plotset(d, c(l4, l6), co) 
        plotset(d, c(l4, l7), co) 
        plotset(d, c(l5, l6), co) 
        plotset(d, c(l5, l7), co) 
        plotset(d, c(l6, l7), co) 

        plotset(d, c(l1, l4), co) 
        plotset(d, c(l2, l5), co) 
        plotset(d, c(l3, l6), co) 
        plotset(d, c(l4, l7), co) 
        plotset(d, c(l5, l8), co) 
        plotset(d, c(l6, l9), co) 
        plotset(d, c(l7, l10), co) 
        plotset(d, c(l8, l11), co) 

        MassConstrict(d, cutoffs=co, strange=TRUE, lower.bound=15, upper.bound=125)

    dev.off()
}

# =================
# UTILITY FUNCTIONS
# =================

MassConstrict <- function(d, cutoffs=c(40), ...){
    strata <- unique(d$mrca)
    strata <- strata[-length(strata)]
    for (stratum in strata){
        for (cutoff in cutoffs){
            try(print(RevenantPlot(d, strata=stratum, cutoff=cutoff, ...)), silent=FALSE)
        }
    }
}
