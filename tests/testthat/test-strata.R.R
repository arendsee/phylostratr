context("strata.R")

strata <- Strata(
    tree = ape::read.tree(system.file('extdata', 'yeast', 'tree', package='phylostratr')),
    data = list(faa=list(
             'Saccharomyces cerevisiae'   = 'yeast/cerevisiae.faa',
             'Saccharomyces paradoxus'    = 'yeast/paradoxus.faa',
             'Saccharomyces mikatae'      = 'yeast/mikatae.faa',
             'Saccharomyces kudriavzevii' = 'yeast/kudriavzevii.faa',
             'Saccharomyces arboricola'   = 'yeast/arboricola.faa',
             'Saccharomyces eubayanus'    = 'yeast/eubayanus.faa',
             'Saccharomyces uvarum'       = 'yeast/uvarum.faa'
           )),
    focal_name = 'Saccharomyces cerevisiae',
    focal_id = 4932
)

test_that("strata_convert works", {
  expect_equal(
    strata %>%
      strata_convert(target='tip', to='id') %>%
      strata_convert(target='tip', to='name'),
    strata
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

