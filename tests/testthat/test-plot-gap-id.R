# plot(type = "gap", id = ...) draws the chosen treated units' own effects,
# as gsynth 1.2.x did for one unit (gsynth #106). plot.gsynth() passes `id`
# to fect's plot.fect(), which draws this plot from fect 2.4.7 on (fect
# #162). Before, `id` was ignored and the gap plot showed the average effect
# over all treated units.

## The estimates the plot draws: x, y, ymin and ymax of its point-range
## layers, sorted by x.
.gg_points <- function(p) {
  geoms <- vapply(p$layers, function(l) class(l$geom)[1], "")
  d <- do.call(rbind, lapply(which(geoms == "GeomPointrange"), function(i)
    ggplot2::layer_data(p, i)[, c("x", "y", "ymin", "ymax")]))
  d <- d[order(d$x), ]
  rownames(d) <- NULL
  d
}

## By hand from the fit: at each time relative to the treatment onset
## (fit$T.on), the mean of fit$eff over the chosen units' cells.
.gg_by_hand <- function(fit, ids) {
  j <- match(as.character(ids), as.character(fit$id))
  rel <- fit$T.on[, j, drop = FALSE]
  eff <- fit$eff[, j, drop = FALSE]
  ok <- !is.na(rel) & !is.na(eff)
  data.frame(x = sort(unique(rel[ok])),
             y = as.numeric(tapply(eff[ok], rel[ok], mean)))
}

test_that("G14: the gap plot with a treated id draws that unit's effects, as in gsynth 1.2.x", {
  skip_on_cran()
  e <- new.env()
  utils::data("gsynth", package = "gsynth", envir = e)
  fit <- suppressMessages(gsynth(Y ~ D + X1 + X2, data = e$simdata,
                                 index = c("id", "time"), force = "two-way",
                                 r = 2, CV = FALSE, se = FALSE, seed = 1,
                                 parallel = FALSE))
  msgs <- character(0)
  p <- withCallingHandlers(plot(fit, type = "gap", id = 101),
                           message = function(m) {
                             msgs <<- c(msgs, conditionMessage(m))
                             invokeRestart("muffleMessage")
                           })
  pts <- .gg_points(p)
  j <- which(fit$id == 101)
  expect_equal(pts$x, as.numeric(fit$T.on[, j]))
  expect_equal(pts$y, as.numeric(fit$eff[, j]), tolerance = 1e-12)
  ## unit 101's gap in its first three treated periods (the example of #106)
  expect_equal(round(pts$y[pts$x %in% 1:3], 3), c(0.342, 4.332, 4.308))
  expect_true(all(is.na(pts$ymin)))
  expect_identical(p$labels$title, "id = 101")
  expect_length(msgs, 1)
  expect_match(msgs, "Unit-level uncertainty is not shown", fixed = TRUE)
  ## the default type of plot.gsynth() is "gap"
  expect_equal(.gg_points(suppressMessages(plot(fit, id = 101))), pts)
  ## several ids: their average at each relative time
  p3 <- suppressMessages(plot(fit, type = "gap", id = c(101, 102, 103)))
  ref <- .gg_by_hand(fit, c(101, 102, 103))
  expect_equal(.gg_points(p3)$x, ref$x)
  expect_equal(.gg_points(p3)$y, ref$y, tolerance = 1e-12)
  expect_identical(p3$labels$title, "Average over 3 units")
  ## id = NULL: the average effect over all treated units, as before
  p0 <- suppressMessages(plot(fit, type = "gap"))
  expect_equal(.gg_points(p0)$y, as.numeric(fit$att), tolerance = 1e-12)
  expect_identical(p0$labels$title, "Estimated Dynamic Treatment Effects")
})

test_that("G15: with SEs the gap plot for an id draws no interval; control and unknown ids stop", {
  skip_on_cran()
  e <- new.env()
  utils::data("gsynth", package = "gsynth", envir = e)
  fit <- suppressMessages(gsynth(Y ~ D + X1 + X2, data = e$simdata,
                                 index = c("id", "time"), force = "two-way",
                                 r = 2, CV = FALSE, se = TRUE, nboots = 20,
                                 seed = 1, parallel = FALSE))
  ## the average gap plot has intervals, the plot for one unit none
  expect_true(all(is.finite(.gg_points(suppressMessages(plot(fit, type = "gap")))$ymin)))
  pts <- .gg_points(suppressMessages(plot(fit, type = "gap", id = 102)))
  expect_equal(pts$y, as.numeric(fit$eff[, fit$id == 102]), tolerance = 1e-12)
  expect_true(all(is.na(pts$ymin)) && all(is.na(pts$ymax)))
  expect_error(plot(fit, type = "gap", id = 106),
               "Unit(s) in \"id\" never treated (control units): 106.", fixed = TRUE)
  expect_error(plot(fit, type = "gap", id = 999),
               "Unit(s) in \"id\" not in the data: 999.", fixed = TRUE)
  expect_error(plot(fit, type = "gap", id = c(101, 999)),
               "Unit(s) in \"id\" not in the data: 999.", fixed = TRUE)
})
