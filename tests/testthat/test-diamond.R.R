context("blast.R")


tree <- ape::rtree(5) %>% set_node_names
tree$tip.label <- letters[1:5]
strata <- Strata(
  tree=tree,
  focal_species='a',
  data=list(faa=list(
    a = file.path('sample_data', 'a.faa'),
    b = file.path('sample_data', 'b.faa'),
    c = file.path('sample_data', 'c.faa'),
    d = file.path('sample_data', 'd.faa'),
    e = file.path('sample_data', 'e.faa')
  ))
)

strata_results <- NULL
test_that("strata_diamond", {
  expect_equal(
    {
      strata_results <<- strata_diamond(
        strata,
        blast_args=list(verbose=FALSE)
      )
      strata_results@data$blast_result
    },
    list(
      a='a.tab',
      b='b.tab',
      c='c.tab',
      d='d.tab',
      e='e.tab'
    )
  )
})