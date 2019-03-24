context("tree_building.R")

data(atree)
prok <- uniprot_sample_prokaryotes()

test_that("lineages_to_phylo", {
  skip_on_travis()
  lin <- taxizedb::classification(c(3702, 9606))
  expect_equal(tree_names(lineages_to_phylo(lin)), c('3702', '9606', '2759'))
})
