# Smoke tests for v1.5.0 modern-theme pass-through args on plot.gsynth():
# legacy.style, highlight, highlight.fill. Also verifies that
# type = "loading.overlap" routes through to fect cleanly.

data(simdata, package = "gsynth")

test_that("plot(out, type = 'gap', legacy.style = TRUE) returns a plot", {
  skip_on_cran()
  out <- gsynth(Y ~ D + X1 + X2, data = simdata,
                index = c("id", "time"), se = FALSE, r = 2,
                parallel = FALSE)
  p <- plot(out, type = "gap", legacy.style = TRUE)
  expect_true(inherits(p, "ggplot") || inherits(p, "gtable") ||
              is.list(p))
})

test_that("plot(out, type = 'gap', legacy.style = FALSE) returns a plot", {
  skip_on_cran()
  out <- gsynth(Y ~ D + X1 + X2, data = simdata,
                index = c("id", "time"), se = FALSE, r = 2,
                parallel = FALSE)
  p <- plot(out, type = "gap", legacy.style = FALSE)
  expect_true(inherits(p, "ggplot") || inherits(p, "gtable") ||
              is.list(p))
})

test_that("plot(out, type = 'loading.overlap') runs without error", {
  skip_on_cran()
  out <- gsynth(Y ~ D + X1 + X2, data = simdata,
                index = c("id", "time"), se = FALSE, r = 2,
                parallel = FALSE)
  expect_error(
    plot(out, type = "loading.overlap"),
    NA  # expect no error
  )
})

test_that("plot(out, highlight = NULL, highlight.fill = FALSE) is the default path", {
  skip_on_cran()
  out <- gsynth(Y ~ D + X1 + X2, data = simdata,
                index = c("id", "time"), se = FALSE, r = 2,
                parallel = FALSE)
  p <- plot(out, type = "gap",
            highlight = NULL, highlight.fill = FALSE)
  expect_true(inherits(p, "ggplot") || inherits(p, "gtable") ||
              is.list(p))
})

test_that("plot(out, highlight.fill = TRUE) runs without error", {
  skip_on_cran()
  out <- gsynth(Y ~ D + X1 + X2, data = simdata,
                index = c("id", "time"), se = FALSE, r = 2,
                parallel = FALSE)
  expect_error(
    plot(out, type = "gap", highlight.fill = TRUE),
    NA
  )
})
