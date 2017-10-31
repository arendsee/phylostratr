# Reads a csv file containing headers:
# qgb|qlocus|species|phylostratum|mrca|sciname|qlen|palen|malen|salen|pscore|mscore|sscore 
# Orders csv by locus, ascending phylostratum, and species 
PrepareDataframe <- function(masked, unmasked){
    prepare_ <- function(f){
        d <- read.csv(f, header=TRUE, row.names=NULL, stringsAsFactors=TRUE)
        d <- d[order(d$qgb, -d$phylostratum, d$hspecies), ] 
        return(d)
    }
    md <- prepare_(masked)
    ud <- prepare_(unmasked)
    # dataframe necessary for plot functions
    d <- data.frame (
            qgb=md$qgb,
            hspecies=factor(md$hspecies, levels=as.character(unique(md$hspecies))),
            mrca=factor(md$mrca_sciname, levels=as.character(unique(md$mrca_sciname))),
            ps=md$phylostratum,
            mmeval=md$mevalue, # Masked, max hsp evalue
            umeval=ud$mevalue, # Unmasked, max hsp evalue
            mpscore=md$pscore, # Masked, bestpath scores
            mmscore=md$mscore, # Masked, max hsp scores
            upscore=ud$pscore, # Unmasked, bestpath scores
            umscore=ud$mscore  # Unmasked, max hsp scores
    )
    return(d)
}

OrderedFullDataframe <- function(f){
    d <- read.csv(f, header=TRUE, row.names=NULL, stringsAsFactors=TRUE)
    d <- d[order(d$qgb, -d$phylostratum, d$hspecies), ] 
    return(d)
}

EmptyDataframe <- function(is.molten=TRUE){
    if(is.molten){
        return(data.frame(
            qgb=character,
            hspecies=character,
            mrca=character,
            ps=numeric,
            variable=character,
            value=numeric
        ))
    }
    else{
        return(data.frame(
            qgb=character,
            hspecies=character,
            mrca=character,
            ps=numeric,
            mmeval=numeric,
            umeval=numeric,
            mpscore=numeric,
            mmscore=numeric,
            upscore=numeric,
            umscore=numeric
        ))
    }
}
