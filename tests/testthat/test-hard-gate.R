# Tests for the v1.5.0 hard-gate alignment with fect v2.2.0+ and
# the public commitment in Xu (2026), "A Response to Green and Aronow,"
# Appendix A.5: inference = "parametric" combined with IFE-EM or MC
# now errors instead of silently coercing to bootstrap. gsynth stops
# before fect does, with a message in gsynth's arguments (#47).

data("gsynth", package = "gsynth")  # loads simdata and turnout

## Evaluates `expr`, expects an error, and checks that the message names
## gsynth's arguments and remedies, not fect's (vartype, method, ...).
.expect_gate_error <- function(expr) {
  err <- tryCatch(suppressPackageStartupMessages(expr),
                  error = function(e) e)
  expect_s3_class(err, "error")
  msg <- if (inherits(err, "error")) conditionMessage(err) else ""
  expect_match(msg, "inference = \"parametric\"", fixed = TRUE)
  expect_match(msg, "nonparametric", fixed = TRUE)
  expect_match(msg, "jackknife", fixed = TRUE)
  expect_false(grepl("vartype", msg, fixed = TRUE))
  expect_false(grepl("time.component.from", msg, fixed = TRUE))
  expect_false(grepl("method =", msg, fixed = TRUE))
}

test_that("inference = 'parametric' + estimator = 'ife' errors", {
  skip_on_cran()
  .expect_gate_error(
    gsynth(Y ~ D + X1 + X2, data = simdata,
           index = c("id", "time"), estimator = "ife",
           inference = "parametric", se = TRUE, nboots = 50,
           r = 2, CV = FALSE, parallel = FALSE)
  )
})

test_that("inference = 'parametric' + estimator = 'mc' errors", {
  skip_on_cran()
  .expect_gate_error(
    gsynth(Y ~ D + X1 + X2, data = simdata,
           index = c("id", "time"), estimator = "mc",
           inference = "parametric", se = TRUE, nboots = 50,
           parallel = FALSE)
  )
})

test_that("inference = 'parametric' + EM = TRUE errors", {
  skip_on_cran()
  .expect_gate_error(
    gsynth(Y ~ D + X1 + X2, data = simdata,
           index = c("id", "time"), EM = TRUE,
           inference = "parametric", se = TRUE, nboots = 50,
           r = 2, CV = FALSE, parallel = FALSE)
  )
})

test_that("inference = 'bootstrap' with IFE-EM works (no hard gate)", {
  skip_on_cran()
  out <- suppressPackageStartupMessages(
    gsynth(Y ~ D + X1 + X2, data = simdata,
           index = c("id", "time"), estimator = "ife",
           inference = "bootstrap", se = TRUE, nboots = 50,
           r = 2, CV = FALSE, parallel = FALSE)
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
           r = 2, CV = FALSE, parallel = FALSE)
  )
  expect_s3_class(out, "gsynth")
  expect_equal(out$vartype, "jackknife")
})

test_that("inference = 'parametric' with GSC (default) works", {
  skip_on_cran()
  out <- gsynth(Y ~ D + X1 + X2, data = simdata,
                index = c("id", "time"), estimator = "gsynth",
                inference = "parametric", se = TRUE, nboots = 50,
                r = 2, CV = FALSE, parallel = FALSE)
  expect_s3_class(out, "gsynth")
  expect_equal(out$vartype, "parametric")
})

## Calls that fect allows stay allowed (the gate mirrors fect's condition)

test_that("parametric + estimator = 'ife' with never-treated controls works", {
  skip_on_cran()
  out <- suppressPackageStartupMessages(
    gsynth(Y ~ D + X1 + X2, data = simdata,
           index = c("id", "time"), estimator = "ife",
           time.component.from = "nevertreated",
           inference = "parametric", se = TRUE, nboots = 20,
           r = 2, CV = FALSE, parallel = FALSE, seed = 1)
  )
  expect_s3_class(out, "gsynth")
  expect_equal(out$vartype, "parametric")
})

test_that("parametric + estimator = 'ife' or 'mc' with se = FALSE works", {
  skip_on_cran()
  out_ife <- suppressPackageStartupMessages(
    gsynth(Y ~ D + X1 + X2, data = simdata,
           index = c("id", "time"), estimator = "ife",
           inference = "parametric", se = FALSE,
           r = 2, CV = FALSE, parallel = FALSE)
  )
  out_mc <- suppressPackageStartupMessages(
    gsynth(Y ~ D + X1 + X2, data = simdata,
           index = c("id", "time"), estimator = "mc",
           inference = "parametric", se = FALSE,
           parallel = FALSE, seed = 1)
  )
  expect_s3_class(out_ife, "gsynth")
  expect_s3_class(out_mc, "gsynth")
  expect_null(out_ife$vartype)
  expect_null(out_mc$vartype)
})

## ci.method = "basic" with the jackknife. fect refuses this pairing
## whatever `se` is; gsynth stops first, in its own arguments, and says
## that "normal" is the jackknife's only interval.

## Evaluates `expr`, expects an error, and checks the message wording.
.expect_jackknife_basic_error <- function(expr) {
  err <- tryCatch(suppressPackageStartupMessages(expr),
                  error = function(e) e)
  expect_s3_class(err, "error")
  msg <- if (inherits(err, "error")) conditionMessage(err) else ""
  expect_match(msg, "ci.method = \"basic\"", fixed = TRUE)
  expect_match(msg, "inference = \"jackknife\"", fixed = TRUE)
  expect_match(msg, "ci.method = \"normal\" is its only interval", fixed = TRUE)
  expect_false(grepl("vartype", msg, fixed = TRUE))
}

test_that("inference = 'jackknife' + ci.method = 'basic' errors", {
  skip_on_cran()
  .expect_jackknife_basic_error(
    gsynth(Y ~ D + X1 + X2, data = simdata,
           index = c("id", "time"), inference = "jackknife",
           ci.method = "basic", se = TRUE,
           r = 2, CV = FALSE, parallel = FALSE)
  )
  ## also with se = FALSE (fect refuses the pairing either way)
  .expect_jackknife_basic_error(
    gsynth(Y ~ D + X1 + X2, data = simdata,
           index = c("id", "time"), inference = "jackknife",
           ci.method = "basic", se = FALSE,
           r = 2, CV = FALSE, parallel = FALSE)
  )
  ## and for estimator = "ife"
  .expect_jackknife_basic_error(
    gsynth(Y ~ D + X1 + X2, data = simdata,
           index = c("id", "time"), estimator = "ife",
           inference = "jackknife", ci.method = "basic", se = TRUE,
           r = 2, CV = FALSE, parallel = FALSE)
  )
})

test_that("the jackknife with ci.method = 'normal' works", {
  skip_on_cran()
  out <- gsynth(Y ~ D + X1 + X2, data = simdata,
                index = c("id", "time"), inference = "jackknife",
                ci.method = "normal", se = TRUE,
                r = 2, CV = FALSE, parallel = FALSE)
  expect_s3_class(out, "gsynth")
  expect_equal(out$vartype, "jackknife")
})

test_that("ci.method = 'basic' with parametric or nonparametric inference works", {
  skip_on_cran()
  out_par <- gsynth(Y ~ D + X1 + X2, data = simdata,
                    index = c("id", "time"), inference = "parametric",
                    ci.method = "basic", se = TRUE, nboots = 20,
                    r = 2, CV = FALSE, parallel = FALSE, seed = 1)
  out_np <- gsynth(Y ~ D + X1 + X2, data = simdata,
                   index = c("id", "time"), inference = "nonparametric",
                   ci.method = "basic", se = TRUE, nboots = 20,
                   r = 2, CV = FALSE, parallel = FALSE, seed = 1)
  expect_s3_class(out_par, "gsynth")
  expect_s3_class(out_np, "gsynth")
  expect_equal(out_par$vartype, "parametric")
  expect_equal(out_np$vartype, "bootstrap")
})

test_that("calls that stop elsewhere keep their messages", {
  skip_on_cran()
  d <- simdata
  d$w <- 1 + d$id %% 3
  d$w2 <- 1 + d$id %% 4
  ## gsynth's weight check comes first
  expect_error(
    gsynth(Y ~ D + X1 + X2, data = d, index = c("id", "time"),
           inference = "jackknife", ci.method = "basic", se = TRUE,
           weight = "w", W.agg = "w2", r = 2, CV = FALSE, parallel = FALSE),
    "`weight` and `W.agg` name different columns", fixed = TRUE
  )
  ## ci.method values other than "normal" and "basic" are fect's to reject
  expect_error(
    gsynth(Y ~ D + X1 + X2, data = simdata, index = c("id", "time"),
           inference = "jackknife", ci.method = "percentile", se = TRUE,
           r = 2, CV = FALSE, parallel = FALSE),
    "ci.method = \"percentile\" is not supported in fect()", fixed = TRUE
  )
})
