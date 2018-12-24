context("blast.R::read_blast")

test_that("Raise error when missing required columns", {
  # NOTE: warning is expected
  expect_error(read_blast('sample_data/blast-header-missing_col.tab'))
})
test_that("Correctly handle files with or without headers", {
  good <- data.frame(
    qseqid='a', sseqid='b',
    qstart=42L, qend=82L,
    sstart=99L, send=109L,
    evalue=0.001, score=30,
    staxid="9606",
    stringsAsFactors=FALSE
  )
  good_no_taxid <- good[, -9]
  expect_equal(
    read_blast('sample_data/blast-header-no_taxid.tab', with_taxid=FALSE),
    good_no_taxid
  )
  # NOTE: warning is expected
  expect_error(
    read_blast('sample_data/blast-header-no_taxid.tab', with_taxid=TRUE)
  )
  expect_equal(
    read_blast('sample_data/blast-header-taxid.tab'),
    good
  )
  expect_equal(
    read_blast('sample_data/blast-no_header-taxid.tab', col_names=FALSE),
    good
  )
  expect_equal(
    read_blast('sample_data/blast-no_header-no_taxid.tab', with_taxid=FALSE, col_names=FALSE),
    good_no_taxid
  )
  expect_equal(
    read_blast('sample_data/blast-no_header-no_taxid.tab', col_names=FALSE, taxid=9606L),
    good
  )
})
