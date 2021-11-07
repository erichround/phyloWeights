context("Result values")
library(phyloWeights)

multiPhylo_ABCD <- 
  c(read.tree(text = "(((A:0.2,B:0.2,C:0.2):1.8,D:2):0.3);"),
    read.tree(text = "(((A:1,B:1,C:1):1,D:2):0.3);"), 
    read.tree(text = "(((A:1.8,(B:1,C:1):0.8):0.2,D:2):0.3);")
  )

data_ABCD <- data.frame(tip = c("A", "B", "C", "D"), 
                        is_SOV = c(1, 1, 1, 0), 
                        n_consonants = c(18, 20, 22, 40),
                        stringsAsFactors = FALSE)

test_that("BM has right values", {
  expect_equal(round(BM(multiPhylo_ABCD[[1]]), 5),
               c(A = 0.18562, B = 0.18562, C = 0.18562, D = 0.44314))
  expect_equal(round(BM(multiPhylo_ABCD[[2]]), 5),
               c(A = 0.23175, B = 0.23175, C = 0.23175, D = 0.30476))
  expect_equal(round(BM(multiPhylo_ABCD[[3]]), 5),
               c(A = 0.26797, B = 0.22886, C = 0.22886, D = 0.27431))
})

test_that("ACL has right values", {
  expect_equal(round(ACL(multiPhylo_ABCD[[1]]), 5),
               c(A = 0.17241, B = 0.17241, C = 0.17241, D = 0.48276))
  expect_equal(round(ACL(multiPhylo_ABCD[[2]]), 5),
               c(A = 0.2, B = 0.2, C = 0.2, D = 0.4))
  expect_equal(round(ACL(multiPhylo_ABCD[[3]]), 5),
               c(A = 0.28384, B = 0.19651, C = 0.19651, D = 0.32314))
})

test_that("phylo_average has right values", {
  a <- phylo_average(multiPhylo_ABCD, data_ABCD)
  expect_equal(round(a$ACL_averages$is_SOV, 5), c(0.51724, 0.6, 0.67686))
  expect_equal(round(a$BM_averages$is_SOV, 5), c(0.55686, 0.69524, 0.72569))
  expect_equal(round(a$ACL_averages$n_consonants, 5), c(29.65517, 28, 26.28821))
  expect_equal(round(a$BM_averages$n_consonants, 5), 
               c(28.8628, 26.09524, 25.40796))
})