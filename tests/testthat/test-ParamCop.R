test_that("Parametric copula fit", {

  #Check summary table
  #expect_message(fitPar)

  # check the scale and dependency parameter are positive

  expect_true(all(fitPar[,6]>=0))
  expect_equal(rownames(fitPar)[2], "frank")
})
