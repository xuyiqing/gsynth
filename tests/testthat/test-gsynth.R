# Tests for the main gsynth() function

data(simdata, package = "gsynth")

# --- Basic functionality ---

test_that("gsynth() default method returns valid object without SE", {
  skip_on_cran()
  out <- gsynth(Y ~ D + X1 + X2, data = simdata,
                index = c("id", "time"), se = FALSE, r = 2,
                force = "two-way", parallel = FALSE)

  expect_s3_class(out, "gsynth")
  expect_equal(out$method, "gsynth")
  expect_true(is.numeric(out$att.avg))
  expect_length(out$att.avg, 1)
  expect_true(is.numeric(out$att))
  expect_true(length(out$att) > 0)
  expect_s3_class(out$data, "data.frame")
})

test_that("gsynth() EM method returns valid object", {
  skip_on_cran()
  out <- gsynth(Y ~ D + X1 + X2, data = simdata,
                index = c("id", "time"), EM = TRUE, se = FALSE,
                r = 2, parallel = FALSE)

  expect_s3_class(out, "gsynth")
  expect_equal(out$method, "ife")
})

test_that("gsynth() MC method returns valid object", {
  skip_on_cran()
  out <- gsynth(Y ~ D + X1 + X2, data = simdata,
                index = c("id", "time"), estimator = "mc",
                se = FALSE, parallel = FALSE)

  expect_s3_class(out, "gsynth")
  expect_equal(out$method, "mc")
})

test_that("gsynth() with parametric SE returns uncertainty estimates", {
  skip_on_cran()
  out <- gsynth(Y ~ D + X1 + X2, data = simdata,
                index = c("id", "time"), se = TRUE, nboots = 50,
                r = 2, parallel = FALSE)

  expect_false(is.null(out$est.avg))
  expect_true(is.matrix(out$est.att))
  expect_true(all(c("ATT", "S.E.", "CI.lower", "CI.upper") %in%
                    colnames(out$est.att)))
  expect_equal(out$vartype, "parametric")
})

test_that("gsynth() with EM + SE uses bootstrap by default", {
  skip_on_cran()
  out <- expect_no_warning(
    gsynth(Y ~ D + X1 + X2, data = simdata,
           index = c("id", "time"), EM = TRUE, se = TRUE,
           nboots = 50, r = 2, parallel = FALSE)
  )

  expect_equal(out$vartype, "bootstrap")
})

test_that("gsynth() EM + parametric inference produces override warning", {
  skip_on_cran()
  expect_warning(
    out <- gsynth(Y ~ D + X1 + X2, data = simdata,
                  index = c("id", "time"), EM = TRUE, se = TRUE,
                  inference = "parametric", nboots = 50, r = 2,
                  parallel = FALSE),
    "nonparametric bootstrap"
  )
  expect_equal(out$vartype, "bootstrap")
})

test_that("gsynth() nonparametric inference is remapped", {
  skip_on_cran()
  out <- gsynth(Y ~ D + X1 + X2, data = simdata,
                index = c("id", "time"), inference = "nonparametric",
                se = TRUE, nboots = 50, r = 2, parallel = FALSE)

  expect_true(out$vartype %in% c("bootstrap", "jackknife"))
})

test_that("gsynth() cross-validation selects factors", {
  skip_on_cran()
  out <- gsynth(Y ~ D + X1 + X2, data = simdata,
                index = c("id", "time"), CV = TRUE, r = c(0, 3),
                se = FALSE, parallel = FALSE)

  expect_s3_class(out, "gsynth")
  expect_false(is.null(out$r.cv))
})

test_that("gsynth() formula interface works with named argument", {
  skip_on_cran()
  out <- gsynth(formula = Y ~ D + X1 + X2, data = simdata,
                index = c("id", "time"), se = FALSE, r = 2,
                force = "two-way", parallel = FALSE)

  expect_s3_class(out, "gsynth")
  expect_true(is.numeric(out$att.avg))
})

test_that("gsynth() Y/D interface works without formula", {
  skip_on_cran()
  out <- gsynth(data = simdata, Y = "Y", D = "D",
                X = c("X1", "X2"), index = c("id", "time"),
                se = FALSE, r = 2, force = "two-way",
                parallel = FALSE)

  expect_s3_class(out, "gsynth")
  expect_true(is.numeric(out$att.avg))
})

# --- Call storage ---

test_that("gsynth() stores call object correctly", {
  skip_on_cran()
  out <- gsynth(Y ~ D + X1 + X2, data = simdata,
                index = c("id", "time"), se = FALSE, r = 2,
                parallel = FALSE)

  expect_true(is.call(out$call))
})

# --- Edge cases ---

test_that("gsynth() force = 'none' works", {
  skip_on_cran()
  out <- gsynth(Y ~ D + X1 + X2, data = simdata,
                index = c("id", "time"), se = FALSE, r = 2,
                force = "none", parallel = FALSE)

  expect_s3_class(out, "gsynth")
})

test_that("gsynth() force = 'time' works", {
  skip_on_cran()
  out <- gsynth(Y ~ D + X1 + X2, data = simdata,
                index = c("id", "time"), se = FALSE, r = 2,
                force = "time", parallel = FALSE)

  expect_s3_class(out, "gsynth")
})

test_that("gsynth() r = 0 (no factors) works", {
  skip_on_cran()
  out <- gsynth(Y ~ D + X1 + X2, data = simdata,
                index = c("id", "time"), se = FALSE, r = 0,
                parallel = FALSE)

  expect_s3_class(out, "gsynth")
})

test_that("gsynth() seed produces reproducible results", {
  skip_on_cran()
  out1 <- gsynth(Y ~ D + X1 + X2, data = simdata,
                 index = c("id", "time"), se = FALSE, r = 2,
                 parallel = FALSE, seed = 42)
  out2 <- gsynth(Y ~ D + X1 + X2, data = simdata,
                 index = c("id", "time"), se = FALSE, r = 2,
                 parallel = FALSE, seed = 42)

  expect_identical(out1$att.avg, out2$att.avg)
})

test_that("gsynth() normalize = TRUE works", {
  skip_on_cran()
  out <- gsynth(Y ~ D + X1 + X2, data = simdata,
                index = c("id", "time"), se = FALSE, r = 2,
                normalize = TRUE, parallel = FALSE)

  expect_s3_class(out, "gsynth")
})
