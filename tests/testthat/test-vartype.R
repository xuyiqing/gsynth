# G2 (#37, #60, #75): the fit records the inference it used in `vartype`
# (one string), and the stored call is the call as typed. Before, gsynth
# copied `inference` into the call as `vartype`: NULL when inference was not
# typed (so fect's effect() treated parametric draws as bootstrap draws), a
# symbol when passed through a variable, and an argument that gsynth() does
# not have (so update() failed and print() showed it).

data("gsynth", package = "gsynth")  # loads simdata and turnout

.fit_se <- function(...) {
  suppressMessages(
    gsynth(Y ~ D + X1 + X2, data = simdata, index = c("id", "time"),
           force = "two-way", CV = FALSE, r = 2, parallel = FALSE,
           seed = 1, ...)
  )
}

test_that("vartype holds the inference used, as one string", {
  skip_on_cran()
  expect_identical(.fit_se(se = TRUE, nboots = 20)$vartype, "parametric")
  expect_identical(
    .fit_se(se = TRUE, nboots = 20, inference = "nonparametric")$vartype,
    "bootstrap"
  )
  expect_identical(.fit_se(se = TRUE, inference = "jackknife")$vartype,
                   "jackknife")
  expect_null(.fit_se(se = FALSE)$vartype)
})

test_that("the stored call is the call as typed (no vartype)", {
  skip_on_cran()
  typed <- .fit_se(se = TRUE, nboots = 20, inference = "nonparametric")
  untyped <- .fit_se(se = TRUE, nboots = 20)
  for (f in list(typed, untyped)) {
    expect_false("vartype" %in% names(f$call))
  }
  expect_identical(typed$call$inference, "nonparametric")

  ## print() shows the call as typed
  printed <- capture.output(print(typed))
  expect_false(any(grepl("vartype", printed, fixed = TRUE)))
  expect_true(any(grepl("inference = \"nonparametric\"", printed,
                        fixed = TRUE)))
})

test_that("update() works on se = TRUE fits, with or without inference", {
  skip_on_cran()
  typed <- .fit_se(se = TRUE, nboots = 20, inference = "nonparametric")
  untyped <- .fit_se(se = TRUE, nboots = 20)

  up_typed <- suppressMessages(update(typed, nboots = 10))
  expect_s3_class(up_typed, "gsynth")
  expect_identical(up_typed$vartype, "bootstrap")
  expect_equal(up_typed$call$nboots, 10)

  up_untyped <- suppressMessages(update(untyped, nboots = 10))
  expect_identical(up_untyped$vartype, "parametric")
})

test_that("an inference passed through a variable is resolved", {
  skip_on_cran()
  ## Called directly (through a wrapper's `...`, match.call() records a
  ## symbol as `..N`); update() looks up `inf` where it is called.
  inf <- "parametric"
  by_var <- suppressMessages(
    gsynth(Y ~ D + X1 + X2, data = simdata, index = c("id", "time"),
           force = "two-way", CV = FALSE, r = 2, parallel = FALSE,
           seed = 1, se = TRUE, nboots = 20, inference = inf)
  )
  expect_identical(by_var$vartype, "parametric")
  expect_false("vartype" %in% names(by_var$call))
  up_var <- suppressMessages(update(by_var, nboots = 10))
  expect_identical(up_var$vartype, "parametric")
})

test_that("print() returns the fit invisibly", {
  skip_on_cran()
  for (f in list(.fit_se(se = TRUE, nboots = 20), .fit_se(se = FALSE))) {
    out <- NULL
    capture.output(out <- withVisible(print(f)))
    expect_false(out$visible)
    expect_identical(out$value, f)
  }
})

## Needs fect >= 2.4.6 with its fix A2 (effect() reads the inference from
## the fit's vartype). Fails with fect dev 412d7ae, which reads it from the
## call (now always without vartype): parametric intervals are then built
## as if from bootstrap draws and are centered near zero.
test_that("effect() intervals contain the estimates (fect A2)", {
  skip_on_cran()
  for (f in list(.fit_se(se = TRUE, nboots = 50),
                 .fit_se(se = TRUE, nboots = 50, inference = "parametric"))) {
    est <- suppressMessages(effect(f, cumu = TRUE))$effect.est.att
    expect_true(all(est[, "CI.lower"] <= est[, "ATT"] &
                      est[, "ATT"] <= est[, "CI.upper"]))
  }
})
