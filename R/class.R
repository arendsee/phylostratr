setOldClass("phylo")
setClass(
  "Strata",
  representation(
    tree = "phylo",
    data = "list",
    focal_species = "character" 
  )
)

Strata <- function(focal_species, tree, data=list()){
  new('Strata',
    tree = tree,
    data = data,
    focal_species = focal_species
  )
}
