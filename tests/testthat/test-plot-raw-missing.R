# plot(type = "raw") and plot(type = "missing") are drawn by panelView from
# the data and variable names stored in the fit. Before, they rebuilt the
# model from x$call$formula (failing after a Y =/D = call, with a formula
# in a variable, or with as.formula(...)), ignored `id`, and drew the plot
# inside the call, so `plot(fit, type = "raw") + theme_dark()` drew twice
# (#58, #26).

data("gsynth", package = "gsynth")  # loads simdata and turnout

skip_on_cran()

# Shared fixtures: the same model given four ways
fml <- Y ~ D + X1 + X2
f_frm <- suppressMessages(
  gsynth(Y ~ D + X1 + X2, data = simdata, index = c("id", "time"),
         r = 2, CV = FALSE, parallel = FALSE)
)
f_str <- suppressMessages(
  gsynth(data = simdata, Y = "Y", D = "D", X = c("X1", "X2"),
         index = c("id", "time"), r = 2, CV = FALSE, parallel = FALSE)
)
f_var <- suppressMessages(
  gsynth(fml, data = simdata, index = c("id", "time"),
         r = 2, CV = FALSE, parallel = FALSE)
)
f_asf <- suppressMessages(
  gsynth(stats::as.formula("Y ~ D + X1 + X2"), data = simdata,
         index = c("id", "time"), r = 2, CV = FALSE, parallel = FALSE)
)

## Runs `fun()` with a new PDF file as the current device. Returns the
## value of `fun()` (or its error), the number of pages drawn in the file,
## and whether the device list and the current device were the same after
## `fun()` as before. The PDF device is closed and the file deleted even
## on failure.
.draw_test <- function(fun) {
  tf <- tempfile(fileext = ".pdf")
  old <- grDevices::dev.cur()
  grDevices::pdf(tf)
  dev <- grDevices::dev.cur()
  on.exit({
    if (dev %in% grDevices::dev.list()) grDevices::dev.off(dev)
    if (old > 1L && old %in% grDevices::dev.list()) grDevices::dev.set(old)
    unlink(tf)
  }, add = TRUE)
  before <- list(grDevices::dev.list(), grDevices::dev.cur())
  value <- tryCatch(fun(), error = function(e) e)
  after <- list(grDevices::dev.list(), grDevices::dev.cur())
  grDevices::dev.off(dev)
  lines <- readLines(tf, warn = FALSE)
  pages <- sum(grepl("/Type /Page ", lines, fixed = TRUE, useBytes = TRUE))
  list(value = value, pages = pages, same_devices = identical(before, after))
}

.types <- c("raw", "missing")

test_that("raw and missing plots work however the model was given", {
  fits <- list(formula = f_frm, strings = f_str, variable = f_var,
               as.formula = f_asf)
  for (nm in names(fits)) {
    for (type in .types) {
      res <- .draw_test(function() plot(fits[[nm]], type = type))
      expect_s3_class(res$value, "ggplot")
    }
  }
})

test_that("the string, variable and as.formula fits give the formula fit's plot", {
  for (type in .types) {
    p_frm <- .draw_test(function() plot(f_frm, type = type))$value
    for (f in list(f_str, f_var, f_asf)) {
      p <- .draw_test(function() plot(f, type = type))$value
      expect_equal(p$data, p_frm$data)
    }
  }
})

test_that("id restricts the units in raw and missing plots", {
  ids <- c(101, 102, 110)
  for (type in .types) {
    p_all <- .draw_test(function() plot(f_frm, type = type))$value
    p_id <- .draw_test(function() plot(f_frm, type = type, id = ids))$value
    expect_lt(nrow(p_id$data), nrow(p_all$data))
  }
  ## the same plot data as panelview() itself gives for these units
  pv_mis <- .draw_test(function() {
    panelView::panelview(data = simdata, Y = "Y", D = "D", X = c("X1", "X2"),
                         index = c("id", "time"), pre.post = TRUE, id = ids,
                         main = "Treatment Status and Missing Data")
  })$value
  p_mis <- .draw_test(function() plot(f_frm, type = "missing", id = ids))$value
  expect_equal(p_mis$data, pv_mis$data)
  pv_raw <- .draw_test(function() {
    panelView::panelview(data = simdata, Y = "Y", D = "D", X = c("X1", "X2"),
                         index = c("id", "time"), type = "outcome", id = ids,
                         main = "Raw Data")
  })$value
  p_raw <- .draw_test(function() plot(f_frm, type = "raw", id = ids))$value
  expect_equal(p_raw$data, pv_raw$data)
})

test_that("raw and missing plots are drawn only when printed", {
  for (type in .types) {
    ## nothing is drawn while the result is assigned
    res <- .draw_test(function() {
      p <- plot(f_frm, type = type)
      invisible(p)
    })
    expect_equal(res$pages, 0)
    expect_true(res$same_devices)
    ## adding a theme and printing draws exactly one page
    res <- .draw_test(function() {
      print(plot(f_frm, type = type) + ggplot2::theme_dark())
    })
    expect_equal(res$pages, 1)
    expect_true(res$same_devices)
  }
})

test_that("raw and missing plots return the ggplot visibly", {
  for (type in .types) {
    v <- .draw_test(function() withVisible(plot(f_frm, type = type)))$value
    expect_true(v$visible)
    expect_s3_class(v$value, "ggplot")
  }
})

test_that("the device list and current device are left as they were", {
  ## with the devices open at this point (none in a fresh session)
  before <- list(grDevices::dev.list(), grDevices::dev.cur())
  p <- plot(f_frm, type = "missing")
  expect_identical(list(grDevices::dev.list(), grDevices::dev.cur()), before)

  ## also when panelview() stops with an error
  bad <- f_frm
  bad$index <- c("id", "no_such_column")
  res <- .draw_test(function() plot(bad, type = "raw"))
  expect_s3_class(res$value, "error")
  expect_true(res$same_devices)
})
