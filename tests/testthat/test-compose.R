#R

context("Test compose peptides")

test_that("compose peptide check", {
  

 
  #table(aa_pool_x8)
  #length(aa_pool_x8)
  #table(aa_pool_1_2_9_10)
  #length(aa_pool_1_2_9_10)
  #table(aa_pool_3_8)
  #length(aa_pool_3_8)
  

  expect_true( length(get_pool_x8()) == 100 )
  expect_true( length(get_pool_1_2_9_10()) == 100 )
  expect_true( length(get_pool_3_8()) == 100 )
  #expect_true( sum(unique(x$Peptide.Sequence) %in% peptides$Peptide.Sequence) == 14)
  #expect_true( abs(0.95 * 5 * 2 * 14) <= nrow(x.sum.cv))
})