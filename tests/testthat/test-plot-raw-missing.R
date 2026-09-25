# G3 (#58, #26): plot(type = "raw") and plot(type = "missing") are built
# from the data and variable names the fit stores, use `id`, and return the
# plot without drawing it. Before, they rebuilt the model from
# x$call$formula (failing after a Y = / D = call, with a formula stored in a
# variable, or with as.formula()), ignored `id`, and drew the plot during
# the call (#58: plot(fit, type = "raw") + theme_dark() drew twice).

data("gsynth", package = "gsynth")  # loads simdata and turnout

## Number of pages drawn on a PDF file device while `expr` is evaluated.
.pages_drawn <- function(expr) {
  tf <- tempfile(fileext = ".pdf")
  grDevices::pdf(tf, compress = FALSE)
  dv <- grDevices::dev.cur()
  on.exit({
    if (dv %in% grDevices::dev.list()) grDevices::dev.off(dv)
    unlink(tf)
  })
  force(expr)
  grDevices::dev.off(dv)
  txt <- readLines(tf, warn = FALSE)
  ## page objects are "/Type /Page"; the page tree is "/Type /Pages"
  sum(grepl("/Type /Page([^s]|$)", txt, useBytes = TRUE))
}

.dev_state <- function() list(grDevices::dev.list(), grDevices::dev.cur())

test_that("raw and missing plots work for every way of naming the variables", {
  skip_on_cran()
  fit_f <- gsynth(Y ~ D + X1 + X2, data = simdata, index = c("id", "time"),
                  force = "two-way", CV = FALSE, r = 2, parallel = FALSE)
  fit_s <- gsynth(Y = "Y", D = "D", X = c("X1", "X2"), data = simdata,
                  index = c("id", "time"),
                  force = "two-way", CV = FALSE, r = 2, parallel = FALSE)
  f <- Y ~ D + X1 + X2
  fit_v <- gsynth(f, data = simdata, index = c("id", "time"),
                  force = "two-way", CV = FALSE, r = 2, parallel = FALSE)
  fit_a <- gsynth(stats::as.formula("Y ~ D + X1 + X2"), data = simdata,
                  index = c("id", "time"),
                  force = "two-way", CV = FALSE, r = 2, parallel = FALSE)

  for (ty in c("raw", "missing")) {
    ref <- plot(fit_f, type = ty)
    expect_true(inherits(ref, "ggplot"))
    for (x in list(fit_s, fit_v, fit_a)) {
      p <- plot(x, type = ty)
      expect_true(inherits(p, "ggplot"))
      ## the same plot as for the formula fit
      expect_identical(p$data, ref$data)
    }
  }
})

test_that("id restricts the units shown in raw and missing plots", {
  skip_on_cran()
  fit <- gsynth(Y ~ D + X1 + X2, data = simdata, index = c("id", "time"),
                force = "two-way", CV = FALSE, r = 2, parallel = FALSE)
  ids <- c(101, 102, 110)

  ## raw: only the outcomes of the chosen units are drawn
  p_all <- plot(fit, type = "raw")
  p_id <- plot(fit, type = "raw", id = ids)
  expect_lt(nrow(p_id$data), nrow(p_all$data))
  expect_true(setequal(p_id$data$outcome, simdata$Y[simdata$id %in% ids]))

  ## missing: one row of cells per chosen unit
  m_all <- plot(fit, type = "missing")
  m_id <- plot(fit, type = "missing", id = ids)
  expect_equal(length(unique(m_all$data$units)), length(unique(simdata$id)))
  expect_equal(length(unique(m_id$data$units)), length(ids))
  expect_equal(nrow(m_id$data), sum(simdata$id %in% ids))
})

test_that("raw and missing plots are returned, not drawn, and draw once", {
  skip_on_cran()
  fit <- gsynth(Y = "Y", D = "D", X = c("X1", "X2"), data = simdata,
                index = c("id", "time"),
                force = "two-way", CV = FALSE, r = 2, parallel = FALSE)
  for (ty in c("raw", "missing")) {
    v <- NULL
    ## nothing is drawn while the plot is built and assigned
    expect_equal(.pages_drawn(v <- withVisible(plot(fit, type = ty))), 0)
    ## returned visibly, like the other plot types
    expect_true(v$visible)
    expect_true(inherits(v$value, "ggplot"))
    ## adding a theme and printing draws one page (#58)
    expect_equal(
      .pages_drawn(print(plot(fit, type = ty) + ggplot2::theme_dark())), 1
    )
  }
})

test_that("raw and missing plots leave the graphics devices as they were", {
  skip_on_cran()
  fit <- gsynth(Y ~ D + X1 + X2, data = simdata, index = c("id", "time"),
                force = "two-way", CV = FALSE, r = 2, parallel = FALSE)

  ## no device open (when none is open now): none is left open
  if (is.null(grDevices::dev.list())) {
    had_rplots <- file.exists("Rplots.pdf")
    p <- plot(fit, type = "raw")
    p <- plot(fit, type = "missing")
    opened <- grDevices::dev.list()
    ## (an old build draws on a new default device: close it)
    for (d in opened) grDevices::dev.off(d)
    if (!had_rplots) unlink("Rplots.pdf")
    expect_null(opened)
  }

  ## a device open: it stays open and current
  tf <- tempfile(fileext = ".pdf")
  grDevices::pdf(tf)
  dv <- grDevices::dev.cur()
  on.exit({
    if (dv %in% grDevices::dev.list()) grDevices::dev.off(dv)
    unlink(tf)
  }, add = TRUE)
  before <- .dev_state()
  p <- plot(fit, type = "raw")
  expect_identical(.dev_state(), before)
  p <- plot(fit, type = "missing")
  expect_identical(.dev_state(), before)

  ## ... also when panelView stops with an error
  bad <- fit
  bad[["index"]] <- c("no_such_column", "time")
  expect_error(plot(bad, type = "raw"))
  expect_identical(.dev_state(), before)
})

test_that("a fit without stored data gets a clear error", {
  skip_on_cran()
  fit <- gsynth(Y ~ D + X1 + X2, data = simdata, index = c("id", "time"),
                force = "two-way", CV = FALSE, r = 2, parallel = FALSE)
  fit[["data"]] <- NULL
  expect_error(plot(fit, type = "raw"), "does not store the data",
               fixed = TRUE)
})
