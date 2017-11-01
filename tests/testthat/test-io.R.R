context("io.R")

test_that("The sample data can be loaded", {
  expect_warning(load_hittable(system.file('extdata', 'araport11', 'araport11_subset99.tab', package='phylostratr')))
})

capture_warnings({
  d <- load_hittable(system.file('extdata', 'araport11', 'araport11_subset99.tab', package='phylostratr'))
})
