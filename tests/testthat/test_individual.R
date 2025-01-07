test_that("create_individual works", {
  individual <- gpredomicsR:::create_individual()
  expect_s3_class(individual, "Individual")
  expect_true(!is.null(individual$ptr))
})