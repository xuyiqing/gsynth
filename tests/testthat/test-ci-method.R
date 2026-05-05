# Tests for the v1.5.0 ci.method argument on gsynth()
# (forwards fect's ci.method, added in fect v2.4.2).

data(simdata, package = "gsynth")

test_that("ci.method = 'normal' produces Wald CIs (point +/- z * SE)", {
  skip_on_cran()
  set.seed(11)
  out <- gsynth(Y ~ D + X1 + X2, data = simdata,
                index = c("id", "time"), force = "two-way",
                CV = FALSE, r = 2, se = TRUE,
                inference = "bootstrap", nboots = 100,
                ci.method = "normal", parallel = FALSE)
  ## CI bounds should be symmetric around the point estimate up to
  ## numerical noise: lo + hi - 2 * point ~ 0.
  est <- out$est.att
  midpoint <- (est[, "CI.lower"] + est[, "CI.upper"]) / 2
  expect_true(all(abs(midpoint - est[, "ATT"]) < 1e-8))
})

test_that("ci.method = 'basic' produces non-Wald CIs that differ from normal", {
  skip_on_cran()
  set.seed(11)
  out_normal <- gsynth(Y ~ D + X1 + X2, data = simdata,
                       index = c("id", "time"), force = "two-way",
                       CV = FALSE, r = 2, se = TRUE,
                       inference = "bootstrap", nboots = 100,
                       ci.method = "normal", parallel = FALSE)
  set.seed(11)
  out_basic  <- suppressWarnings(
    gsynth(Y ~ D + X1 + X2, data = simdata,
           index = c("id", "time"), force = "two-way",
           CV = FALSE, r = 2, se = TRUE,
           inference = "bootstrap", nboots = 100,
           ci.method = "basic", parallel = FALSE)
  )
  ## Same point estimates and SEs (those don't depend on ci.method).
  expect_equal(out_normal$est.att[, "ATT"],  out_basic$est.att[, "ATT"])
  expect_equal(out_normal$est.att[, "S.E."], out_basic$est.att[, "S.E."])
  ## CIs should differ on at least one event time.
  expect_false(isTRUE(all.equal(out_normal$est.att[, "CI.lower"],
                                out_basic$est.att[, "CI.lower"])))
})

test_that("ci.method default is 'normal' (matches fect default)", {
  skip_on_cran()
  set.seed(11)
  out_default <- gsynth(Y ~ D + X1 + X2, data = simdata,
                        index = c("id", "time"), force = "two-way",
                        CV = FALSE, r = 2, se = TRUE,
                        inference = "bootstrap", nboots = 100,
                        parallel = FALSE)
  set.seed(11)
  out_normal  <- gsynth(Y ~ D + X1 + X2, data = simdata,
                        index = c("id", "time"), force = "two-way",
                        CV = FALSE, r = 2, se = TRUE,
                        inference = "bootstrap", nboots = 100,
                        ci.method = "normal", parallel = FALSE)
  expect_equal(out_default$est.att, out_normal$est.att)
})

test_that("invalid ci.method is caught by fect", {
  skip_on_cran()
  expect_error(
    gsynth(Y ~ D + X1 + X2, data = simdata,
           index = c("id", "time"), force = "two-way",
           CV = FALSE, r = 2, se = TRUE,
           inference = "bootstrap", nboots = 50,
           ci.method = "not-a-method", parallel = FALSE)
  )
})

test_that("ci.method does not change SE / point estimate", {
  skip_on_cran()
  set.seed(11)
  out_normal <- gsynth(Y ~ D + X1 + X2, data = simdata,
                       index = c("id", "time"), force = "two-way",
                       CV = FALSE, r = 2, se = TRUE,
                       inference = "bootstrap", nboots = 100,
                       ci.method = "normal", parallel = FALSE)
  set.seed(11)
  out_basic <- suppressWarnings(
    gsynth(Y ~ D + X1 + X2, data = simdata,
           index = c("id", "time"), force = "two-way",
           CV = FALSE, r = 2, se = TRUE,
           inference = "bootstrap", nboots = 100,
           ci.method = "basic", parallel = FALSE)
  )
  expect_equal(out_normal$att, out_basic$att)
})
