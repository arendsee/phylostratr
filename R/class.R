setOldClass("phylo")
setClass(
  "Strata",
  representation(
    tree = "phylo",
    data = "list",
    focal_species = "character" 
  )
)
