#R

context("Test compose peptides")

test_that("compose peptide check", {
  


  expect_true( length(get_pool_x8()) == 100 )
  expect_true( length(get_pool_1_2_9_10()) == 100 )
  expect_true( length(get_pool_3_8()) == 100 )
  
})