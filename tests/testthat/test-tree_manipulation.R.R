context("tree_manipulation.R")

data(atree)
prok <- uniprot_sample_prokaryotes()

test_that("tree_names", {
  expect_equal(tree_names(atree), c(paste0('t', 1:10), paste0('n', 11:19)))
})

test_that("tree_size", {
  expect_equal(tree_size(atree), 19)
})

# This is a particularly important function. It is used to validate the id
# inputs of every tree function that takes IDs. I do not test id-related corner
# cases in any of the other functions, since they all call this one.
test_that("clean_phyid", {
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

test_that("children", {
  # a leaf has no children
  expect_equal(children(atree, 't1'), integer(0))
  # can find leaf children
  expect_equal(children(atree, 'n19'), 7:8)
  # can find node children
  expect_equal(children(atree, 'n11'), c(10,12))
})

test_that("descendent_nodes", {
  expect_equal(descendent_nodes(atree, 'n11'), 1:19)
  expect_equal(descendent_nodes(atree, 'n14'), c(1,2,3,14,15))
})

test_that("descendents", {
  expect_equal(descendents(atree, 'n11'), 1:10)
  expect_equal(descendents(atree, 'n14'), c(1,2,3))
})

test_that("is_root", {
  expect_equal(which(is_root(atree)), 11)
})

test_that("leafs", {
  expect_equal(leafs(atree), 1:10)
})

test_that("lineage", {
  expect_equal(lineage(atree, 5), c(11,12,13,16,17,18,5))
  expect_equal(lineage(atree, 13), c(11,12,13))
  expect_equal(lineage(atree, 11), 11)
})

test_that("merge_phylo", {
  expect_equal({
    subtrees <- list(
      subtree(atree, 'n19'),
      subtree(atree, 'n18'),
      subtree(atree, 'n15')
    )
    tree_names(merge_phylo(atree, subtrees))
  },
    c("t1", "t2", "t4", "t5", "t7", "t8", "n13", "n15", "n16", "n18", "n19")
  )
})

test_that("nleafs", {
  expect_equal(nleafs(atree), 10)
})

test_that("nodes", {
  expect_equal(nodes(atree), 11:19)
})

test_that("parent", {
  expect_identical(parent(atree, 11), NA_integer_)
  expect_equal(parent(atree, 5), 18)
  expect_equal(parent(atree, 17), 16)
  expect_equal(parent(atree, c(5,17,11)), c(18L,16L,NA_integer_))
})

test_that("prune", {
  expect_equal(
    tree_names(prune(atree, c(14,17))),
    c("t7", "t8", "t9", "t10", "n11", "n12", "n19")
  )
  expect_equal(
    tree_names(prune(atree, 13)),
    c("t9", "t10", "n11")
  )
  expect_equal(
    tree_names(prune(atree, 12)),
    c("t10", "n11")
  )
  # pruning the root leads to an empty tree
  expect_equal(
    tree_names(prune(atree, 11)),
    character(0)
  )
})

test_that("set_node_names", {
  expect_equal(
    ape::rtree(5) %>% set_node_names %>% {.$node.label},
    c("n6", "n7", "n8", "n9")
  )
})

test_that("sisters", {
  expect_equal(sisters(atree, 17), 19)
  # root has no sisters
  expect_equal(sisters(atree, 11), integer(0))
})

test_that("sister_trees", {
  expect_equal(names(sister_trees(atree, 17)), "19")
  expect_equal(tree_names(sister_trees(atree, 17)[[1]]), c("t7", "t8", "n19"))
  expect_equal(sister_trees(atree, 11), list())
})

test_that("subtree(collapse=TRUE, descend=TRUE)", {
  # the subset of the root is the input tree
  expect_equal(subtree(atree, 11), atree)
  expect_equal(tree_names(subtree(atree, c(1,4,7))), c('t1', 't4', 't7', 'n13', 'n16'))
  # single tip
  expect_equal(tree_names(subtree(atree, 5)), c('t5', 'n18'))
  # internal nodes propagate to tips
  expect_equal(subtree(atree, 17), subtree(atree, 4:6))
  # nothing gets you nothing
  expect_equal(tree_names(subtree(atree, integer(0))), character(0))
})
test_that("subtree(collapse=FALSE, descend=TRUE)", {
  expect_equal(subtree(atree, 11, collapse=FALSE), atree)
  expect_equal(
    tree_names(subtree(atree, c(1,7), collapse=FALSE)),
    c('t1', 't7', 'n11', 'n12', 'n13', 'n14', 'n15', 'n16', 'n19'))
  expect_equal(
    tree_names(subtree(atree, 7, collapse=FALSE)),
    c('t7', 'n11', 'n12', 'n13', 'n16', 'n19'))
  expect_equal(subtree(atree, 17, collapse=FALSE), subtree(atree, 4:6, collapse=FALSE))
  expect_equal(tree_names(subtree(atree, integer(0), collapse=FALSE)), character(0))
})
test_that("subtree(collapse=TRUE, descend=FALSE)", {
  expect_equal(tree_names(subtree(atree, 11, descend=FALSE)), character(0))
  expect_equal(tree_names(subtree(atree, c(1,4,7), descend=FALSE)), c('t1', 't4', 't7', 'n13', 'n16'))
  expect_equal(tree_names(subtree(atree, 7, descend=FALSE)), c('t7', 'n19'))
  expect_identical(subtree(atree, c(7,17), descend=FALSE), subtree(atree, 7, descend=FALSE))
  expect_equal(tree_names(subtree(atree, integer(0), descend=FALSE)), character(0))
})

test_that("subtree edges are correct", {
  expect_equal(
    subtree(atree, c(1,4,7), descend=FALSE)$edge,
    matrix(c(4,4,5,5,1,5,2,3), ncol=2)
  )
  # check single tip edge
  expect_equal(subtree(atree, 't1')$edge, matrix(c(2,1), ncol=2))
  # check single tip names
  expect_equal(tree_names(subtree(atree, 't1')), c("t1", "n15"))
})

test_that("backbone trees are cool", {
  one_tree <- lineage_to_ancestor_tree(lineage(atree, 5, use_name=TRUE))
  expect_equal(tree_names(one_tree)[parent(one_tree, "t5")], "n18") 
  expect_equal(tree_names(one_tree)[parent(one_tree, "n16")], "n13") 
})
