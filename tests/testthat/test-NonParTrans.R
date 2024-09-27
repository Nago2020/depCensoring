test_that("Nonparametric transformation model check", {

  ## Check the fit of dependent censoring model

  #Check summary table
  #expect_message(output$parameterEstimates)

  expect_true(-1<rhohat|rhohat<1)
  expect_equal(l, 5)
  expect_true(inherits(output, "fitNonPar"))
})
