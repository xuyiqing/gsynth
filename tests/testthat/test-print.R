# Tests for print.gsynth()

data(simdata, package = "gsynth")

skip_on_cran()

# Shared fixtures
out_se <- gsynth(Y ~ D + X1 + X2, data = simdata,
                 index = c("id", "time"), se = TRUE, nboots = 50,
                 r = 2, force = "two-way", parallel = FALSE, seed = 1234)

out_nose <- gsynth(Y ~ D + X1 + X2, data = simdata,
                   index = c("id", "time"), se = FALSE, r = 2,
                   force = "two-way", parallel = FALSE)

test_that("print.gsynth() without SE shows expected output", {
  captured <- capture.output(print(out_nose))
  text <- paste(captured, collapse = "\n")

  expect_true(grepl("Call:", text))
  expect_true(grepl("Average Treatment Effect on the Treated:", text))
  expect_true(grepl("by Period", text))
  expect_true(grepl("Uncertainty estimates not available", text))
})

test_that("print.gsynth() with SE shows expected output", {
  captured <- capture.output(print(out_se))
  text <- paste(captured, collapse = "\n")

  expect_true(grepl("Call:", text))
  expect_true(grepl("Average Treatment Effect on the Treated:", text))
  expect_true(grepl("by Period", text))
  expect_false(grepl("Uncertainty estimates not available", text))
})

test_that("print.gsynth() with covariates shows coefficients", {
  # simdata has X1 and X2 covariates
  captured <- capture.output(print(out_se))
  text <- paste(captured, collapse = "\n")

  expect_true(grepl("Coefficients for the Covariates:", text))
})

test_that("print.gsynth() does not error", {
  expect_no_error(print(out_se))
  expect_no_error(print(out_nose))
})
