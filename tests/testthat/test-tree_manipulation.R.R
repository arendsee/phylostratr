context("tree_manipulation.R")

test_that("clean_phyid", {
  data(atree)
  prok <- uniprot_sample_prokaryotes()
  # edge cases are correctly handled
  expect_true(length(clean_phyid(atree, NULL)) == 0)
  expect_true(length(clean_phyid(atree, character(0))) == 0)
  # incorrect lengths are caught
  expect_error(clean_phyid(atree, character(0), len=1))
  expect_error(clean_phyid(atree, 1, len=2))
  # die on unsupported type
  expect_error(clean_phyid(atree, 1, type="yolo"))
  # automatic numeric indices work
  expect_equal(clean_phyid(atree, 1), 1)
  expect_equal(clean_phyid(atree, 1.000000001), 1)
  expect_equal(clean_phyid(atree, 19), 19)
  expect_error(clean_phyid(atree, 20))
  # automatic character indices work
  expect_equal(clean_phyid(atree, "1"), 1)
  expect_equal(clean_phyid(atree, "19"), 19)
  expect_error(clean_phyid(atree, "20"))
  # typed numeric indices work
  expect_equal(clean_phyid(atree, 1, type='index'), 1)
  expect_equal(clean_phyid(atree, 1.000000001, type='index'), 1)
  expect_equal(clean_phyid(atree, 19, type='index'), 19)
  expect_error(clean_phyid(atree, 20, type='index'))
  # automatic character indices work
  expect_equal(clean_phyid(atree, "1", type='index'), 1)
  expect_equal(clean_phyid(atree, "19", type='index'), 19)
  expect_error(clean_phyid(atree, "20", type='index'))
  # automatic character names work
  expect_equal(clean_phyid(atree, 't5'), 5)
  expect_equal(clean_phyid(atree, 'n15'), 15)
  expect_error(clean_phyid(atree, 'unicorn'))
  # typed character names work
  expect_equal(clean_phyid(atree, 't5', type='name'), 5)
  expect_equal(clean_phyid(atree, 'n15', type='name'), 15)
  expect_error(clean_phyid(atree, 'unicorn', type='name'))
  # types can distinguish between integrel names and integrel ids 
  # use names when possible ('2' is the name of a node in prok)
  expect_equal(clean_phyid(prok, 2, type='auto'), 87)
  expect_equal(clean_phyid(prok, 2, type='name'), 87)
  expect_equal(clean_phyid(prok, 2, type='index'), 2)
})
