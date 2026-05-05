# Tests for the v1.5.0 hard-gate alignment with fect v2.2.0+ and
# the public commitment in Xu (2026), "A Response to Green and Aronow,"
# Appendix A.5: inference = "parametric" combined with IFE-EM or MC
# now errors instead of silently coercing to bootstrap.

data(simdata, package = "gsynth")

test_that("inference = 'parametric' + estimator = 'ife' errors", {
  skip_on_cran()
  expect_error(
    suppressPackageStartupMessages(
      gsynth(Y ~ D + X1 + X2, data = simdata,
             index = c("id", "time"), estimator = "ife",
             inference = "parametric", se = TRUE, nboots = 50,
             r = 2, parallel = FALSE)
    )
  )
})

test_that("inference = 'parametric' + estimator = 'mc' errors", {
  skip_on_cran()
  expect_error(
    suppressPackageStartupMessages(
      gsynth(Y ~ D + X1 + X2, data = simdata,
             index = c("id", "time"), estimator = "mc",
             inference = "parametric", se = TRUE, nboots = 50,
             parallel = FALSE)
    )
  )
})

test_that("inference = 'parametric' + EM = TRUE errors", {
  skip_on_cran()
  expect_error(
    suppressPackageStartupMessages(
      gsynth(Y ~ D + X1 + X2, data = simdata,
             index = c("id", "time"), EM = TRUE,
             inference = "parametric", se = TRUE, nboots = 50,
             r = 2, parallel = FALSE)
    )
  )
})

test_that("inference = 'bootstrap' with IFE-EM works (no hard gate)", {
  skip_on_cran()
  out <- suppressPackageStartupMessages(
    gsynth(Y ~ D + X1 + X2, data = simdata,
           index = c("id", "time"), estimator = "ife",
           inference = "bootstrap", se = TRUE, nboots = 50,
           r = 2, parallel = FALSE)
  )
  expect_s3_class(out, "gsynth")
  expect_equal(out$vartype, "bootstrap")
})

test_that("inference = 'jackknife' with IFE-EM works (no hard gate)", {
  skip_on_cran()
  out <- suppressPackageStartupMessages(
    gsynth(Y ~ D + X1 + X2, data = simdata,
           index = c("id", "time"), estimator = "ife",
           inference = "jackknife", se = TRUE,
           r = 2, parallel = FALSE)
  )
  expect_s3_class(out, "gsynth")
  expect_equal(out$vartype, "jackknife")
})

test_that("inference = 'parametric' with GSC (default) works", {
  skip_on_cran()
  out <- gsynth(Y ~ D + X1 + X2, data = simdata,
                index = c("id", "time"), estimator = "gsynth",
                inference = "parametric", se = TRUE, nboots = 50,
                r = 2, parallel = FALSE)
  expect_s3_class(out, "gsynth")
  expect_equal(out$vartype, "parametric")
})
