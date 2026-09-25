# G4 (#101): `weight` weights only the averaged treatment effects (fect's
# W.agg), as ?gsynth documents and as gsynth <= 1.2.1 did. Before, gsynth
# passed it as fect's W, which also weights the model fit.

data("gsynth", package = "gsynth")  # loads simdata and turnout

## simdata with unit-level weights. With `unbalanced = TRUE`, control units
## 110-115 miss periods 1-4: then a weighted model fit differs from the
## unweighted one (on a complete control panel the fit ignores weights).
.wdata <- function(unbalanced = TRUE) {
  d <- simdata
  d$w <- 1 + (d$id %% 3)
  d$w2 <- 2 * d$w
  if (unbalanced) d <- d[!(d$id %in% 110:115 & d$time %in% 1:4), ]
  d
}

.fit_w <- function(dat, ...) {
  suppressMessages(
    gsynth(Y ~ D + X1 + X2, data = dat, index = c("id", "time"),
           force = "two-way", parallel = FALSE, ...)
  )
}

test_that("weight weights the averaged effects, not the model fit", {
  skip_on_cran()
  ub <- .wdata()
  unw <- .fit_w(ub, r = 2, CV = FALSE)
  wt <- .fit_w(ub, r = 2, CV = FALSE, weight = "w")
  wa <- .fit_w(ub, r = 2, CV = FALSE, W.agg = "w")

  ## the model is the unweighted one
  expect_equal(wt$beta, unw$beta)
  expect_equal(wt$eff, unw$eff)
  ## the averages are the W.agg ones ...
  expect_equal(wt$att.avg, wa$att.avg)
  expect_equal(wt$att, wa$att)
  ## ... that is, the weighted mean of the unweighted effects over the
  ## treated cells
  ids <- as.numeric(as.character(unw$id))
  w <- matrix(1 + (ids %% 3), nrow(unw$eff), ncol(unw$eff), byrow = TRUE)
  treated <- unw$D.dat == 1 & !is.na(unw$eff)
  expect_equal(wt$att.avg,
               sum(w[treated] * unw$eff[treated]) / sum(w[treated]))
  expect_false(isTRUE(all.equal(wt$att.avg, unw$att.avg)))
  ## fect records the roles
  expect_false(wt$W.in.fit)
  expect_true(wt$W.in.agg)
})

test_that("weight and a different W.agg stop; the same column runs", {
  skip_on_cran()
  ub <- .wdata()
  expect_error(
    .fit_w(ub, r = 2, CV = FALSE, weight = "w", W.agg = "w2"),
    "`weight` and `W.agg` name different columns", fixed = TRUE
  )
  same <- .fit_w(ub, r = 2, CV = FALSE, weight = "w", W.agg = "w")
  wa <- .fit_w(ub, r = 2, CV = FALSE, W.agg = "w")
  expect_equal(same$beta, wa$beta)
  expect_equal(same$att.avg, wa$att.avg)
})

## Needs fect >= 2.4.6 with its fix B4 (with CV on, the aggregation weight
## leaked into the model fit). Fails with fect dev 412d7ae.
test_that("weight with CV on leaves the model unweighted (fect B4)", {
  skip_on_cran()
  ub <- .wdata()
  unw <- .fit_w(ub, r = c(0, 5), seed = 1)
  wt <- .fit_w(ub, r = c(0, 5), weight = "w", seed = 1)
  expect_equal(wt$r.cv, unw$r.cv)
  expect_equal(wt$beta, unw$beta)
  expect_equal(wt$eff, unw$eff)
})

## Needs fect >= 2.4.6 with its fix A4 (weighted fits with parametric
## standard errors stopped). Fails with fect dev 412d7ae.
test_that("weight with se = TRUE and parametric inference runs (fect A4)", {
  skip_on_cran()
  d <- .wdata(unbalanced = FALSE)
  wt <- .fit_w(d, r = 2, CV = FALSE, weight = "w",
               se = TRUE, nboots = 20, seed = 1)
  expect_equal(wt$vartype, "parametric")
  expect_equal(unname(wt$est.avg[1, 1]), wt$att.avg)
})
