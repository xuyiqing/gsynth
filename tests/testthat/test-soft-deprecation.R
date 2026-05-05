# Tests for v1.5.0 soft-deprecation messages on IFE-EM, MC, EM = TRUE,
# and effect(). Each message fires once per session via the
# .gsynth_deprecation_seen closure in R/default.R.

data(simdata, package = "gsynth")

# Helper: peel into the package internals to reset the deprecation cache,
# so each test starts from a clean slate. Tests run in the package's own
# environment, so we can access the closure directly.
.reset_deprecation_seen <- function() {
  env <- environment(gsynth:::.gsynth_deprecation_seen)
  env$seen <- character(0)
  invisible(NULL)
}

test_that("estimator = 'ife' emits a deprecation message", {
  skip_on_cran()
  .reset_deprecation_seen()
  msg <- capture.output(
    out <- gsynth(Y ~ D + X1 + X2, data = simdata,
                  index = c("id", "time"), estimator = "ife",
                  se = FALSE, r = 2, parallel = FALSE),
    type = "message"
  )
  expect_match(paste(msg, collapse = "\n"),
               "soft-deprecated", fixed = FALSE)
  expect_match(paste(msg, collapse = "\n"),
               "fect::fect", fixed = TRUE)
  expect_s3_class(out, "gsynth")
})

test_that("estimator = 'mc' emits a deprecation message", {
  skip_on_cran()
  .reset_deprecation_seen()
  msg <- capture.output(
    out <- gsynth(Y ~ D + X1 + X2, data = simdata,
                  index = c("id", "time"), estimator = "mc",
                  se = FALSE, parallel = FALSE),
    type = "message"
  )
  expect_match(paste(msg, collapse = "\n"),
               "soft-deprecated", fixed = FALSE)
  expect_match(paste(msg, collapse = "\n"),
               "fect::fect", fixed = TRUE)
  expect_s3_class(out, "gsynth")
})

test_that("EM = TRUE emits a deprecation message", {
  skip_on_cran()
  .reset_deprecation_seen()
  msg <- capture.output(
    out <- gsynth(Y ~ D + X1 + X2, data = simdata,
                  index = c("id", "time"), EM = TRUE,
                  se = FALSE, r = 2, parallel = FALSE),
    type = "message"
  )
  combined <- paste(msg, collapse = "\n")
  expect_match(combined, "soft-deprecated", fixed = FALSE)
  expect_match(combined, "fect::fect", fixed = TRUE)
  expect_s3_class(out, "gsynth")
})

test_that("deprecation message fires only once per session per key", {
  skip_on_cran()
  .reset_deprecation_seen()

  ## First call: message expected
  msg1 <- capture.output(
    gsynth(Y ~ D + X1 + X2, data = simdata,
           index = c("id", "time"), estimator = "ife",
           se = FALSE, r = 2, parallel = FALSE),
    type = "message"
  )
  ## Second call: message should NOT fire again
  msg2 <- capture.output(
    gsynth(Y ~ D + X1 + X2, data = simdata,
           index = c("id", "time"), estimator = "ife",
           se = FALSE, r = 2, parallel = FALSE),
    type = "message"
  )

  expect_true(any(grepl("soft-deprecated", msg1)))
  expect_false(any(grepl("soft-deprecated", msg2)))
})

test_that("effect() emits a deprecation message pointing to estimand()", {
  skip_on_cran()
  .reset_deprecation_seen()

  ## fect::effect() requires bootstrap/jackknife results; fit with se = TRUE.
  out <- gsynth(Y ~ D + X1 + X2, data = simdata,
                index = c("id", "time"), se = TRUE, nboots = 50, r = 2,
                parallel = FALSE)

  msg <- capture.output(
    cumu <- effect(out, cumu = TRUE, plot = FALSE),
    type = "message"
  )
  combined <- paste(msg, collapse = "\n")
  expect_match(combined, "soft-deprecated", fixed = FALSE)
  expect_match(combined, "estimand", fixed = TRUE)
})
