# Tests for the effect() function

data(simdata, package = "gsynth")

# Shared fixture: gsynth object with SE (required for effect())
# Created once at file scope to avoid repeated slow computation
skip_on_cran()

out <- gsynth(Y ~ D + X1 + X2, data = simdata,
              index = c("id", "time"), se = TRUE, nboots = 50,
              r = 2, force = "two-way", parallel = FALSE, seed = 1234)

test_that("effect() cumulative effect returns valid output", {
  eff <- effect(out, cumu = TRUE)
  expect_false(is.null(eff))
})

test_that("effect() non-cumulative returns valid output", {
  eff <- effect(out, cumu = FALSE)
  expect_false(is.null(eff))
})

test_that("effect() with period range returns valid output", {
  eff <- effect(out, period = c(1, 5))
  expect_false(is.null(eff))
})

test_that("effect() with plot = TRUE does not error", {
  expect_no_error(effect(out, plot = TRUE))
})
