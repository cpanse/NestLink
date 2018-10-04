#R

context("Test compose peptides")

test_that("compose peptide check", {
  


  expect_true( length(get_pool_x8()) == 100 )
  expect_true( length(get_pool_1_2_9_10()) == 100 )
  expect_true( length(get_pool_3_8()) == 100 )
  
})


test_that("check ssrc function",{
  irtPeptide.ssrc <- c(7.33685 , 11.93745 , 29.22845 , 16.65045 , 20.58245 ,
                       27.03345 , 19.42845 , 22.32845 , 22.71845 , 36.67245 ,
                       25.31445 , 34.21845 , 38.108815 , 39.461215 , 40.516115 ,
                       40.488815 , 46.237915 , 41.876215 , 41.476515 , 46.025115 ,
                       7.33685 , 11.59745 , 16.65045 , 20.58245 , 19.42845 ,
                       22.32845 , 22.79845 , 25.57445 , 37.89545 , 40.488815 ,
                       46.111915)
  
  expect_equal(as.numeric(ssrc(as.character(iRTpeptides$peptide))),
               irtPeptide.ssrc, tolerance = 0.01)
})