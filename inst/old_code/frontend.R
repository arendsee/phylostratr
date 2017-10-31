#! /usr/bin/Rscript

home <- '~/ohome/lib/R/'
source(paste0(home, 'visual_functions.R'), chdir=TRUE)
source(paste0(home, 'io.R'), chdir=TRUE)
source(paste0(home, 'analyses.R'), chdir=TRUE)
source(paste0(home, 'phylogenetic_functions.R'), chdir=TRUE)

args <- commandArgs(trailingOnly=TRUE)

masked <- '~/ohome/data/blast_csvs/masked_at.csv'
unmasked <- '~/ohome/data/blast_csvs/unmasked_at.csv'
base.name <- 'at'
d <- PrepareDataframe(masked=masked, unmasked=unmasked)

lin <- c('Arabidopsis','Camelineae','Brassicaceae','malvids','rosids','core eudicotyledons',
         'Magnoliophyta','Spermatophyta','Tracheophyta','Embryophyta','Viridiplantae','Eukaryota')
revenants <- FindRevenants(d, strata=lin, cutoff=40)
m <- read.csv(masked)
m <- m[which(m$qgb %in% revenants), ]
u <- read.csv(masked)
u <- u[which(u$qgb %in% revenants), ]
write.csv(m, 'masked_revenants.csv')
write.csv(u, 'unmasked_revenants.csv')
r <- d[which(d$qgb %in% revenants), ]
plotobo(r, title_=paste0(base.name, '_obo_revenants.pdf')) 

lin <- c('Viridiplantae','Eukaryota')
revenants <- FindRevenants(d, strata=lin, cutoff=40)
write.csv(m, 'masked_ve_revenants.csv')
write.csv(u, 'unmasked_ve_revenants.csv')
r <- d[which(d$qgb %in% revenants), ]
plotobo(r, title_=paste0(base.name, '_obo_ve_revenants.pdf')) 
