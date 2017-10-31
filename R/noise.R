home <- '~/ohome/lib/R/'
source(paste0(home, 'io.R'), chdir=TRUE)

require(reshape2)
u <- OrderedFullDataframe('~/ohome/data/blast_csvs/unmasked_at.csv')
u$npscore <- pmax(0, u$pscore - log2(u$qlen))
u$nmscore <- pmax(0, u$mscore - log2(u$qlen))
u$nsscore <- pmax(0, u$sscore - log2(u$qlen))
m <- melt(u, id.vars=c('qgb', 'qlocus', 'hspecies', 'phylostratum', 'mrca', 'mrca_sciname'))
a <- acast(m, qgb ~ hspecies ~ variable, value.var='value')
a <- a[, as.character(unique(u$hspecies)), ]

z <- apply(a[,,'mscore'], 2, function(x) sum(x == 0) / sum(x < 40))
plot(1:length(z), z)

plotdens <- function(a., score='mscore'){
    par(mfrow=c(7,7), mar=c(0,0,0,0))
    for(i in 1:dim(a)[2]){
        dens <- density(log(a.[,i,score]))
        plot(dens$x, dens$y, type='l', tck=0, xaxt="n", yaxt="n")
    }
}

jpeg('Hs_score-kernel-density-plot.jpeg')
plotdens(a)
dev.off()

Score <- function(dens., u., s., a., b.){
    x <- dens$x
    dens.star <- dnorm(x, u., s.) + dgamma(x, a., b.)  
    rts <- mean((dens.star - dens$y)^2)
    return(rts)
}

plotboth <- function(dens., scale., u., s., a., b.){
    plot(dens.$x, dens.$y, type='l') 
    dens.star <- scale. * dnorm(dens.$x, u., s.) + (1-scale.) * dgamma(dens.$x, a., b.)  
    lines(dens.$x, dens.star, type='l', col='red')
}
