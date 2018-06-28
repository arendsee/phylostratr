context("diversity-algo")

data(atree)

algo1 <- function(...){
  diverse_subtree(FUN=.algo1, ...)$tip.label
}

test_that("without weights", {
  expect_equal(
    lapply(1:10, function(i) algo1(atree, i)),
    list(
      c("t4"), # choose the deepest node
      c("t4", "t10"),
      c("t4", "t9", "t10"),
      c("t1", "t4", "t9", "t10"),
      c("t1", "t4", "t7", "t9", "t10"),
      c("t1", "t3", "t4", "t7", "t9", "t10"),
      c("t1", "t2", "t3", "t4", "t7", "t9", "t10"), # t2 before t6: alphabetic order 
      c("t1", "t2", "t3", "t4", "t6", "t7", "t9", "t10"),
      c("t1", "t2", "t3", "t4", "t6", "t7", "t8", "t9", "t10"),
      c("t1", "t2", "t3", "t4", "t5", "t6", "t7", "t8", "t9", "t10"))
  )
})

test_that("with weights", {
  # the first selected is the one with the highest score
  expect_equal(algo1(atree, 1, weights=c("t3"=1.001)), "t3")
  # low weight loses to high divesity 
  expect_equal(
    algo1(atree, 2, weights=c("t3"=1.001, "t6"=1.001)),
    c("t6", "t10")
  )
  # but high enough scores are chosen over diversity
  expect_equal(
    algo1(atree, 2, weights=c("t3"=10, "t6"=10)),
    c("t3", "t6")
  )
})

test_that("polytomies are well-behaved", {
  unresolved <- ape::read.tree(text="(r,(a,b,c),(e,f,g));")
  # alternate between children of polytomies
  expect_equal(
    lapply(1:7, function(i) algo1(unresolved, i)),
    list(
      c("a"),
      c("a", "e"),
      c("r", "a", "e"),
      c("r", "a", "b", "e"),
      c("r", "a", "b", "e", "f"),
      c("r", "a", "b", "c", "e", "f"),
      c("r", "a", "b", "c", "e", "f", "g"))
  )
})
