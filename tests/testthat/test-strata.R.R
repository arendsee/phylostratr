context("strata.R")

strata <- Strata(
    tree = ape::read.tree(system.file('extdata', 'yeast', 'tree', package='phylostratr')),
    data = list(faa=list(
             'Saccharomyces_cerevisiae'   = 'yeast/cerevisiae.faa',
             'Saccharomyces_paradoxus'    = 'yeast/paradoxus.faa',
             'Saccharomyces_mikatae'      = 'yeast/mikatae.faa',
             'Saccharomyces_kudriavzevii' = 'yeast/kudriavzevii.faa',
             'Saccharomyces_arboricola'   = 'yeast/arboricola.faa',
             'Saccharomyces_eubayanus'    = 'yeast/eubayanus.faa',
             'Saccharomyces_uvarum'       = 'yeast/uvarum.faa'
           )),
    focal_species = 'Saccharomyces_cerevisiae'
) %>% sort_strata

test_that("is_valid_strata catches problems", {
  expect_silent(is_valid_strata(strata)) 
  # missin required field
  expect_error(is_valid_strata(strata, required='ladida')) 
  # species name in data is not in tree
  expect_error({
    bad_strata <- strata
    names(bad_strata@faa)[1] <- "Bob"
    is_valid_strata(bad_strata)
  })
  # wrong focal species
  expect_error({
    bad_strata <- strata
    bad_strata@focal_species <- "Bob"
    is_valid_strata(bad_strata)
  })
})

test_that("strata_convert works", {
  skip_on_travis()
  # name-->id-->name is not guaranteed to result in an identical strata,
  # since more than one name may be associated with a given id.
  # in this case, I get the same names back but without the underscores
  expect_equal(
    strata %>%
      strata_convert(target='tip', to='id') %>%
      strata_convert(target='tip', to='name'),
    {
      strata2 <- strata
      strata2@tree$tip.label <- gsub("_", " ", strata2@tree$tip.label)
      strata2@focal_species <- sub("_", " ", strata2@focal_species)
      names(strata2@data$faa) <- gsub("_", " ", names(strata2@data$faa))
      strata2
    }
  )
  expect_equal(
      as.numeric(strata_convert(strata, target='tip', to='id')@tree$tip.label),
      c(4932, 27291, 114525, 114524, 706196, 1080349, 230603)
  )
  expect_equal(
      as.numeric(names(strata_convert(strata, target='tip', to='id')@data$faa)),
      c(4932, 27291, 114525, 114524, 706196, 1080349, 230603)
  )
})

test_that("strata_apply over identity changes nothing", {
  ident <- function(x, ...) x
  # The tree is unchanged
  expect_true(ape::all.equal.phylo(strata_apply(strata, ident)@tree, strata@tree))
  # The data are unchaged
  expect_true(setequal(strata_apply(strata, ident)@data$faa, strata@data$faa))
  # Focal species is unchanged
  expect_equal(strata_apply(strata, ident)@focal_species, strata@focal_species)
})

test_that("GitHub issue #6", {
  skip_on_travis()
  taxids=c("4362", "16681", "1194090", "1142394")
  tree <- phylostratr::lineages_to_phylo(taxizedb::classification(taxids))
  strata <- phylostratr::Strata(focal_species=16681, tree)
  expect_equal(names(phylostratr::strata_map(strata, identity)), c("131567", "4360", "16681"))
})
