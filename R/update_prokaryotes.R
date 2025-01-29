# New functions by Laura Tibbs-Cortes Dec 2024
# Rationale: as of Dec 2024, 30 out of the 85 suggested prokaryote species
# have warnings (see https://github.com/arendsee/phylostratr/issues/34)
# To fix this, I am making a function to choose a diverse set of prokaryotes
# of desired size.
# The original phylostratr publication (doi: 10.1093/bioinformatics/btz171) states:
# "For both case studies, we selected a diverse set of 85 prokaryotic species, 
# comprised of one species from each class of Bacteria and Archaea, 
# sampled at random from the UniProt reference proteomes."
# However, there are currently 121 distinct classes in Bacteria alone according to UniProt.
# Therefore, I cannot select 85 species and simultaneously represent every class, 
# so I am making a function to build a new tree of desired size 
# and with the ability to provide weights as desired.
generate_prokaryote_tree <- function(domain.n=50, # approximate number of species desired (may not give exactly this tree size due to rounding)
                                     weights=NULL, 
                                     focal_taxid="3702", # focal taxa just needed to convert to strata object; not used otherwise
                                     current.domain) {  # use domain names from UniProt: '2' = bacteria, '2157' = archaea
  
  if(current.domain %in% c("2", "2157")) {errorCondition("Domain must be set to '2' (Bacteria) or '2157' (Archaea) ")}

# get prokaryote lineage for current domain;
# based on lineage_to_phylo
pro.lineage <- uniprot_downstream_ids(current.domain)%>%
  taxizedb::classification() %>%
  Filter(f=is.data.frame) 

to_edge <- function(xs){   # function from original package 
  # LTC comment: this makes a data table with all the unique connections between e.g., genus (col 1) and sp (col 2), order (col 1) to family (col 2) etc
  # build edge list based on taxonomy IDs
  lapply(xs, function(x){
    matrix(c(head(x$id, -1), tail(x$id, -1)), ncol=2)
  }) %>%
    do.call(what=rbind) %>%
    unique
}

pro.edges <- to_edge(pro.lineage)

# from edgelist_to_phylo code:
from <- pro.edges[,1]
to <- pro.edges[,2]
ids <- unique(c(pro.edges[,1], pro.edges[,2]))

tip.label <- setdiff(ids, from)
node.label <- unique(from)

# make a map from taxonomy ID to internal 1:n ids
idmap <- 1:(length(tip.label) + length(node.label))
names(idmap) <- c(tip.label, node.label)

 # build phylo tree 
tree <- list(
  edges      = matrix(c(idmap[as.character(from)], idmap[as.character(to)]), ncol=2),
  tip.label  = unname(tip.label),
  node.label = unname(node.label),
  Nnode      = length(node.label)
)
class(tree) <- 'phylo'

# completely remove specified taxa
tree <- tree %>%
  drop.tip.phylo(drop.names) 

  # convert to Strata object
strata.tree <- tree %>%
  Strata(
    focal_species = focal_taxid,
    tree       = .,
    data       = list()
  )

# get weights 
if(is.null(weights)){
  current.weights <- rep(1, nleafs(strata.tree))
} else {
  current.weights <- weights[strata.tree@tree$tip.label]
  current.weights[is.na(current.weights)] <- 1
}

# based on diverse_subtree code, but here using to make 
# a diverse tree for current domain of the desired size
current.subtree <- diverse_subtree(tree=strata.tree, n=domain.n,
                                   weights=current.weights,collapse=F) 

return(current.subtree)
}

# Use this function to add the custom tree to the analysis
# based on use_recommended_prokaryotes but enables using custom trees
use_custom_prokaryote_tree <- function(x, # x is a eukaryotic tree
                                       prokaryote.phylo) { # prokaryote.phyl is a phylo object of the desired prokaryotic tree. Generate customizable trees using generate_prokaryote_tree
  # extract the Eukaryota branch
  x@tree <- subtree(x@tree, '2759', #2759 = eukarya
                    type='name')
  # get 'cellular_organism -> Eukaryota' tree (just these two nodes)
  root <- taxizedb::classification(2759)[[1]]$id %>% lineage_to_ancestor_tree
  # bind the prokaryotes to 'cellular_organism'
  prokaryote.phylo <- ape::bind.tree(root, pro.tree)
  # bind the eukaryote input tree to 'Eukaryota'
  prokaryote.phylo <- ape::bind.tree(prokaryote.phylo, x@tree, where=which(tree_names(prokaryote.phylo) == '2759'))
  x@tree <- prokaryote.phylo
  x
}
