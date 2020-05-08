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

test_that("strata_besthits", {
  expect_true(
    # results must be a list of lists
    strata_besthits(strata_results)@data$besthit %>% sapply(is.list) %>% all
  )
  expect_true(
    # results must be a list of lists of data.frames
    strata_besthits(strata_results)@data$besthit %>% sapply(is.data.frame) %>% all
  )
})

test_that("merge_besthits", {
  expect_true({
    besthits <- strata_besthits(strata_results) %>% merge_besthits
    # the columns are correct
    all(c('staxid', 'qseqid', 'evalue', 'score', 'mrca', 'ps') %in% names(besthits)) &&
      # all staxa are represented
      setequal(letters[1:5], unique(besthits$staxid))
  })
})

test_that("check_hit_table works", {
  expect_error(check_hit_table(5))
  expect_error({
    strata_besthits(strata_results) %>% merge_besthits %>%
      dplyr::select(-qseqid) %>% check_hit_table
  })
  expect_error({
    strata_besthits(strata_results) %>% merge_besthits %>%
      dplyr::select(-mrca) %>%
      check_hit_table(has_mrca=TRUE)
  })
  expect_error({
    strata_besthits(strata_results) %>% merge_besthits %>%
      dplyr::select(-ps) %>%
      check_hit_table(has_ps=TRUE)
  })
})

unlink('a.tab')
unlink('b.tab')
unlink('c.tab')
unlink('d.tab')
unlink('e.tab')
unlink('blastdb', recursive=TRUE)
