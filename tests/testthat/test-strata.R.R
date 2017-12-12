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

test_that("strata_convert works", {
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
