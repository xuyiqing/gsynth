# Tests for the v1.5.0 placeboTest argument on gsynth()
# (forwards fect's placebo machinery). gsynth does not expose
# carryoverTest because gsynth does not allow treatment reversals,
# which fect's carryover test requires.

data(simdata, package = "gsynth")

test_that("placeboTest = TRUE produces placebo slots on the fit", {
  skip_on_cran()
  out <- gsynth(Y ~ D + X1 + X2, data = simdata,
                index = c("id", "time"), force = "two-way",
                CV = FALSE, r = 2, se = TRUE,
                inference = "bootstrap", nboots = 50,
                placeboTest = TRUE, placebo.period = c(-2, 0),
                parallel = FALSE)
  expect_true(isTRUE(out$placeboTest))
  expect_false(is.null(out$est.placebo))
  expect_false(is.null(out$att.placebo))
})

test_that("placeboTest = FALSE leaves placebo slots empty", {
  skip_on_cran()
  out <- gsynth(Y ~ D + X1 + X2, data = simdata,
                index = c("id", "time"), force = "two-way",
                CV = FALSE, r = 2, se = TRUE,
                inference = "bootstrap", nboots = 50,
                parallel = FALSE)
  expect_false(isTRUE(out$placeboTest))
})

test_that("placebo.period defaults to NULL when placeboTest = FALSE", {
  skip_on_cran()
  ## Just verifying no error fires when only one of the two args is set.
  out <- gsynth(Y ~ D + X1 + X2, data = simdata,
                index = c("id", "time"), force = "two-way",
                CV = FALSE, r = 2, se = FALSE,
                placebo.period = c(-2, 0),
                parallel = FALSE)
  expect_s3_class(out, "gsynth")
})

test_that("placebo.period accepts a length-1 integer", {
  skip_on_cran()
  out <- gsynth(Y ~ D + X1 + X2, data = simdata,
                index = c("id", "time"), force = "two-way",
                CV = FALSE, r = 2, se = TRUE,
                inference = "bootstrap", nboots = 50,
                placeboTest = TRUE, placebo.period = -2,
                parallel = FALSE)
  expect_true(isTRUE(out$placeboTest))
})

test_that("carryoverTest is not exposed on gsynth()", {
  skip_on_cran()
  ## gsynth() does not accept carryoverTest in its formals (gsynth requires
  ## staggered no-reversal panels; carryover test requires reversals).
  fmls <- names(formals(gsynth))
  expect_false("carryoverTest" %in% fmls)
  expect_false("carryover.period" %in% fmls)
})

test_that("placeboTest is forwarded only when explicitly TRUE", {
  skip_on_cran()
  ## Default fit (placeboTest = FALSE) and explicit FALSE produce
  ## numerically identical fits.
  set.seed(42)
  out_default <- gsynth(Y ~ D + X1 + X2, data = simdata,
                        index = c("id", "time"), force = "two-way",
                        CV = FALSE, r = 2, se = FALSE,
                        parallel = FALSE)
  set.seed(42)
  out_explicit <- gsynth(Y ~ D + X1 + X2, data = simdata,
                         index = c("id", "time"), force = "two-way",
                         CV = FALSE, r = 2, se = FALSE,
                         placeboTest = FALSE,
                         parallel = FALSE)
  expect_equal(out_default$att, out_explicit$att)
})

test_that("placebo fit populates est.placebo with point/SE/CI columns", {
  skip_on_cran()
  out <- gsynth(Y ~ D + X1 + X2, data = simdata,
                index = c("id", "time"), force = "two-way",
                CV = FALSE, r = 2, se = TRUE,
                inference = "bootstrap", nboots = 50,
                placeboTest = TRUE, placebo.period = c(-2, 0),
                parallel = FALSE)
  expect_false(is.null(out$est.placebo))
  expect_true(any(grepl("ATT|p.value|S.E.|CI",
                        colnames(as.matrix(out$est.placebo)))))
})
