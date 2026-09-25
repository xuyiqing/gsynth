# Tests for the v1.5.0 hard-gate alignment with fect v2.2.0+ and
# the public commitment in Xu (2026), "A Response to Green and Aronow,"
# Appendix A.5: inference = "parametric" combined with IFE-EM or MC
# now errors instead of silently coercing to bootstrap.

data("gsynth", package = "gsynth")  # loads simdata and turnout

## G5 (#47): gsynth stops before fect, with a message in gsynth's
## argument names (fect's message names vartype, method and
## time.component.from, which gsynth users never set).
.expect_gate_message <- function(err) {
  msg <- conditionMessage(err)
  expect_match(msg, "inference = \"parametric\"", fixed = TRUE)
  expect_match(msg, "nonparametric", fixed = TRUE)
  expect_match(msg, "jackknife", fixed = TRUE)
  expect_false(grepl("vartype", msg, fixed = TRUE))
  expect_false(grepl("time.component.from", msg, fixed = TRUE))
}

test_that("inference = 'parametric' + estimator = 'ife' errors", {
  skip_on_cran()
  err <- expect_error(
    suppressPackageStartupMessages(
      gsynth(Y ~ D + X1 + X2, data = simdata,
             index = c("id", "time"), estimator = "ife",
             inference = "parametric", se = TRUE, nboots = 50,
             CV = FALSE, r = 2, parallel = FALSE)
    )
  )
  .expect_gate_message(err)
})

test_that("inference = 'parametric' + estimator = 'mc' errors", {
  skip_on_cran()
  err <- expect_error(
    suppressPackageStartupMessages(
      gsynth(Y ~ D + X1 + X2, data = simdata,
             index = c("id", "time"), estimator = "mc",
             inference = "parametric", se = TRUE, nboots = 50,
             parallel = FALSE)
    )
  )
  .expect_gate_message(err)
})

test_that("inference = 'parametric' + EM = TRUE errors", {
  skip_on_cran()
  err <- expect_error(
    suppressPackageStartupMessages(
      gsynth(Y ~ D + X1 + X2, data = simdata,
             index = c("id", "time"), EM = TRUE,
             inference = "parametric", se = TRUE, nboots = 50,
             CV = FALSE, r = 2, parallel = FALSE)
    )
  )
  .expect_gate_message(err)
})

test_that("inference = 'bootstrap' with IFE-EM works (no hard gate)", {
  skip_on_cran()
  out <- suppressPackageStartupMessages(
    gsynth(Y ~ D + X1 + X2, data = simdata,
           index = c("id", "time"), estimator = "ife",
           inference = "bootstrap", se = TRUE, nboots = 50,
           CV = FALSE, r = 2, parallel = FALSE)
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
           CV = FALSE, r = 2, parallel = FALSE)
  )
  expect_s3_class(out, "gsynth")
  expect_equal(out$vartype, "jackknife")
})

test_that("inference = 'parametric' with GSC (default) works", {
  skip_on_cran()
  out <- gsynth(Y ~ D + X1 + X2, data = simdata,
                index = c("id", "time"), estimator = "gsynth",
                inference = "parametric", se = TRUE, nboots = 50,
                CV = FALSE, r = 2, parallel = FALSE)
  expect_s3_class(out, "gsynth")
  expect_equal(out$vartype, "parametric")
})

## G5: calls that fect allows stay allowed (the gate is not widened).
test_that("parametric + estimator = 'ife' + nevertreated runs with se = TRUE", {
  skip_on_cran()
  out <- suppressPackageStartupMessages(
    gsynth(Y ~ D + X1 + X2, data = simdata,
           index = c("id", "time"), estimator = "ife",
           time.component.from = "nevertreated",
           inference = "parametric", se = TRUE, nboots = 20,
           CV = FALSE, r = 2, parallel = FALSE, seed = 1)
  )
  expect_s3_class(out, "gsynth")
  expect_equal(out$vartype, "parametric")
  expect_false(is.null(out$est.avg))
})

test_that("parametric + estimator = 'gsynth' runs even with notyettreated", {
  skip_on_cran()
  ## fect always fits estimator = "gsynth" on never-treated units
  out <- gsynth(Y ~ D + X1 + X2, data = simdata,
                index = c("id", "time"), estimator = "gsynth",
                time.component.from = "notyettreated",
                inference = "parametric", se = TRUE, nboots = 20,
                CV = FALSE, r = 2, parallel = FALSE, seed = 1)
  expect_s3_class(out, "gsynth")
  expect_equal(out$vartype, "parametric")
})

test_that("parametric + estimator = 'ife'/'mc' runs with se = FALSE", {
  skip_on_cran()
  out_ife <- suppressPackageStartupMessages(
    gsynth(Y ~ D + X1 + X2, data = simdata,
           index = c("id", "time"), estimator = "ife",
           inference = "parametric", se = FALSE,
           CV = FALSE, r = 2, parallel = FALSE)
  )
  expect_s3_class(out_ife, "gsynth")
  out_mc <- suppressPackageStartupMessages(
    gsynth(Y ~ D + X1 + X2, data = simdata,
           index = c("id", "time"), estimator = "mc",
           inference = "parametric", se = FALSE,
           parallel = FALSE, seed = 1)
  )
  expect_s3_class(out_mc, "gsynth")
})
