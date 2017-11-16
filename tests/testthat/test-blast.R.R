context("blast.R")

blastresult <- 'a-vs-b.tab'
unlink(blastresult)
unlink('blastdb', recursive=TRUE)

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
    c('qseqid', 'staxid', 'evalue', 'score')
  )
  expect_equal(
    nrow(readr::read_tsv(blastresult)),
    563
  )
})

unlink('blastdb', recursive=TRUE)
unlink('a-vs-b.tab')


strata <- list(
  ps1 = list(a = file.path('sample_data', 'a.faa')),
  ps2 = list(),
  ps3 = list(b = file.path('sample_data', 'b.faa')),
  ps4 = list(c = file.path('sample_data', 'c.faa')),
  ps5 = list(d = file.path('sample_data', 'd.faa'),
       e = file.path('sample_data', 'e.faa')),
  ps6 = list()
)

strata_results <- NULL

test_that("strata_blast", {
  expect_equal(
    {
      strata_results <<- strata_blast(
        file.path('sample_data', 'a.faa'),
        strata,
        blast_args=list(verbose=FALSE)
      )
      strata_results
    },
    list(
      ps1 = c(a='a.tab'),
      ps2 = character(0),
      ps3 = c(b='b.tab'),
      ps4 = c(c='c.tab'),
      ps5 = c(d='d.tab', e='e.tab'),
      ps6 = character(0)
    )
  )
})

test_that("strata_besthits", {
  expect_true(
    # results must be a list of lists
    strata_besthits(strata_results) %>% sapply(is.list) %>% all
  )
  expect_true(
    # results must be a list of lists of data.frames
    strata_besthits(strata_results) %>% sapply(sapply, is.data.frame) %>% unlist %>% all
  )
})

test_that("merge_besthits", {
  expect_true({
    besthits <- strata_besthits(strata_results) %>% merge_besthits
    all(c('staxid', 'qseqid', 'evalue', 'score', 'mrca', 'ps') %in% names(besthits)) &&
    nrow(besthits) == 425
  })
})

unlink('a.tab')
unlink('b.tab')
unlink('c.tab')
unlink('d.tab')
unlink('e.tab')
unlink('blastdb', recursive=TRUE)
