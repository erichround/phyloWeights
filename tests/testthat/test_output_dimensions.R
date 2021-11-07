context("Output dimensions")
library(phyloWeights)

multiPhylo_ABCD <- 
  c(read.tree(text = "(((A:0.2,B:0.2,C:0.2):1.8,D:2):0.3);"),
    read.tree(text = "(((A:1,B:1,C:1):1,D:2):0.3);"), 
    read.tree(text = "(((A:1.8,B:1.8,C:1.8):0.2,D:2):0.3);"), 
    read.tree(text = "((((A:1,B:1):0.8,C:1.8):0.2,D:2):0.3);"),
    read.tree(text = "(((A:1.8,(B:1,C:1):0.8):0.2,D:2):0.3);"),
    read.tree(text = "(((A:1,B:1):1,(C:1.8,D:1.8):0.2):0.3);")
  )

data_ABCD <- data.frame(tip = c("A", "B", "C", "D"), 
                        is_SOV = c(1, 1, 1, 0), 
                        n_consonants = c(18, 20, 22, 40),
                        stringsAsFactors = FALSE)

test_that("BM has right number of values", {
  expect_equal(length(BM(multiPhylo_ABCD[[1]])), 4)
  expect_error(BM(multiPhylo_ABCD)) 
})

test_that("ACL has right number of values", {
  expect_equal(length(ACL(multiPhylo_ABCD[[1]])), 4)
  expect_error(ACL(multiPhylo_ABCD)) 
})

test_that("phylo_average has right elements", {
  a <- phylo_average(multiPhylo_ABCD[[1]], data_ABCD)
  b <- phylo_average(multiPhylo_ABCD, data_ABCD)
  expect_equal(names(a),
               c("phy", "data", "ACL_weights", "BM_weights", "ACL_averages", 
                "BM_averages"))
  expect_equal(class(a$phy), "phylo")
  expect_equal(class(a$data), "data.frame")
  expect_equal(class(a$ACL_weights), "data.frame")
  expect_equal(class(a$BM_weights), "data.frame")
  expect_equal(class(a$ACL_averages), "data.frame")
  expect_equal(class(a$BM_averages), "data.frame")
  expect_equal(nrow(a$BM_averages), 1)
  
  expect_equal(names(b),
               c("phy", "data", "ACL_weights", "BM_weights", "ACL_averages", 
                 "BM_averages"))
  expect_equal(class(b$phy), "multiPhylo")
  expect_equal(class(b$data), "data.frame")
  expect_equal(class(b$ACL_weights), "data.frame")
  expect_equal(class(b$BM_weights), "data.frame")
  expect_equal(class(b$ACL_averages), "data.frame")
  expect_equal(class(b$BM_averages), "data.frame")
  expect_equal(nrow(b$BM_averages), 6)
})