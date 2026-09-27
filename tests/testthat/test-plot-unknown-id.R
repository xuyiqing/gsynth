# plot(type = "raw") and plot(type = "missing") with an `id` that is not a
# unit of the data stop with a message naming it (gsynth #107). Before, the
# id was passed to panelView, which printed "List of units removed from
# dataset" and then stopped with "dim(X) must have a positive length", or
# drew the plot without the unknown unit.

test_that("G13: raw and missing plots stop on an unknown id and name it", {
  skip_on_cran()
  e <- new.env()
  utils::data("gsynth", package = "gsynth", envir = e)
  fit <- suppressMessages(gsynth(Y ~ D + X1 + X2, data = e$simdata,
                                 index = c("id", "time"), force = "two-way",
                                 r = 2, CV = FALSE, se = FALSE, seed = 1,
                                 parallel = FALSE))
  grDevices::pdf(file = NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  for (type in c("raw", "missing")) {
    expect_error(plot(fit, type = type, id = 999),
                 "Unit(s) not in the data: 999.", fixed = TRUE)
    ## a known and an unknown id: stop too, naming only the unknown one
    expect_error(plot(fit, type = type, id = c(101, 999)),
                 "Unit(s) not in the data: 999.", fixed = TRUE)
  }
  expect_error(plot(fit, type = "missing", id = "Wisconsin"),
               "Unit(s) not in the data: Wisconsin.", fixed = TRUE)
  ## known ids still draw
  expect_s3_class(plot(fit, type = "raw", id = c(101, 102)), "ggplot")
  expect_s3_class(plot(fit, type = "missing", id = 101), "ggplot")
})
