# Tests for plot.gsynth()

data(simdata, package = "gsynth")

skip_on_cran()

# Shared fixtures
out_se <- gsynth(Y ~ D + X1 + X2, data = simdata,
                 index = c("id", "time"), se = TRUE, nboots = 50,
                 r = 2, force = "two-way", parallel = FALSE, seed = 1234)

out_nose <- gsynth(Y ~ D + X1 + X2, data = simdata,
                   index = c("id", "time"), se = FALSE, r = 2,
                   force = "two-way", parallel = FALSE)

# --- fect-delegated plot types ---

test_that("plot.gsynth() type = 'gap' with SE works", {
  expect_no_error(plot(out_se, type = "gap"))
})

test_that("plot.gsynth() type = 'gap' without SE works", {
  expect_no_error(plot(out_nose, type = "gap"))
})

test_that("plot.gsynth() type = 'ct' with SE works", {
  expect_no_error(plot(out_se, type = "ct"))
})

test_that("plot.gsynth() type = 'ct' without SE triggers known fect bug", {
  # Known upstream bug: fect's plot.R checks x$vartype which is NULL when se=FALSE
  # This will be fixed in a future fect release
  expect_error(plot(out_nose, type = "ct"))
})

test_that("plot.gsynth() type = 'factors' works", {
  expect_no_error(plot(out_se, type = "factors", nfactors = 1))
})

test_that("plot.gsynth() type = 'loadings' works", {
  expect_no_error(plot(out_se, type = "loadings", nfactors = 1))
})

# --- panelView plot types ---

test_that("plot.gsynth() type = 'raw' works", {
  expect_no_error(plot(out_se, type = "raw"))
})

test_that("plot.gsynth() type = 'missing' works", {
  expect_no_error(plot(out_se, type = "missing"))
})

# --- Plot parameters ---

test_that("plot.gsynth() custom labels work", {
  expect_no_error(plot(out_se, type = "gap", xlab = "Time", ylab = "Effect"))
})

test_that("plot.gsynth() legendOff works", {
  expect_no_error(plot(out_se, type = "gap", legendOff = TRUE))
})

test_that("plot.gsynth() theme.bw = FALSE works", {
  expect_no_error(plot(out_se, type = "gap", theme.bw = FALSE))
})
