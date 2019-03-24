context("blast.R")

blastresult <- 'a-vs-b.tab'
unlink(blastresult)
unlink('blastdb', recursive=TRUE)

test_that("Can convert BLAST filenames to Strata objects", {
  skip_on_travis()
  # test handling of good data
  data(arabidopsis_strata)
  tree <- subtree(arabidopsis_strata, "3700")@tree
  files <- paste0(tree$tip.label, '.tab')
  names(files) <- tree$tip.label
  x <- .filenames_to_phylo(files, "3702")
  expect_true(all.equal(tree, x@tree))
  expect_equal(x@focal_species, "3702")
  expect_equal(x@data, list(blast_result = files))
  # test handling of bad data (die screaming)
  expect_error(.filenames_to_phylo(c("troll.tab")))
})

test_that("Can make and run blast", {
  expect_equal(
    {
      db <- make_blast_database(file.path('sample_data', 'b.faa')) 
      run_blastp(
        query_fastafile = file.path('sample_data', 'a.faa'),
        subject_taxid = 1234, # currently I don't test the correctness of the taxon id
        blastdb = db,
        blastresult = blastresult,
        nthreads = 1
      )
    },
    blastresult
  )
  expect_true(
    .blastdb_exists(file.path('blastdb', 'b.faa'))
  )
  expect_false(
    .blastdb_exists("foofoo")
  )
  expect_equal(
    names(readr::read_tsv(blastresult)),
    c('qseqid', 'sseqid', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'score', 'staxid')
  )
  expect_equal(
    nrow(readr::read_tsv(blastresult)),
    563
  )
})

unlink('blastdb', recursive=TRUE)
unlink('a-vs-b.tab')

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

test_that("strata_blast", {
  expect_equal(
    {
      strata_results <<- strata_blast(
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
