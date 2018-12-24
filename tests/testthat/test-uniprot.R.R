context("uniprot-access.R")

test_that("Can access uniprot downstream IDs", {
  expect_true({
    # Get Brassicaceae species
    ids <- uniprot_downstream_ids(3700)  
    # As UniProt grows, the number of species will increase, so I just test for
    # the species present as of 2017.
    is.integer(ids) &&
    all( c(81972L, 3702L, 50452L, 3708L, 109376L, 51351L, 81985L, 72664L) %in% ids)
  })
})

# FIXME: this used to be handled gracefully with a warning, but now it raises
# an extemely cryptic error about being unable to read a temporary file (the
# file actually exists and has read permissions but is empty). This used to
# work, I'm not sure what is wrong. 
test_that("uniprot_map2pfam handles missing stuff gracefully", {
  expect_error(uniprot_map2pfam("asdf"))
})

virus <- NULL
test_that("Can retrieve UniProt proteomes", {
  expect_message({
    unlink('virus', recursive=TRUE)
    virus <<- uniprot_retrieve_proteome(128952, dir='virus', dryrun=TRUE)
  })
  expect_false({
    # the dryrun above should not have downloaded anything
    file.exists(virus)
  })
  expect_true({
    # download the ebola virus to your computer
    ebola_prot <- uniprot_retrieve_proteome(128952, dir='virus')
    file.exists(ebola_prot)
  })
  expect_true({
    # repeating the download should be super fast (just returns the existing file)
    system.time(uniprot_retrieve_proteome(128952, dir='virus'))['elapsed'] < 0.1
  })
})

unlink('virus', recursive=TRUE)
