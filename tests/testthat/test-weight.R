# `weight` weights only the averaged treatment effects (the same role as
# `W.agg`), as documented and as in gsynth 1.2.1. It does not enter the
# model fit (`W.est` does). Before, gsynth passed it as fect's `W`, which
# fills both roles (#101).

data("gsynth", package = "gsynth")  # loads simdata and turnout

skip_on_cran()

## simdata with weights that vary by unit and some missing control cells:
## for the gsynth estimator, fit weights change the estimates only when
## control cells are missing
.wdata <- local({
  d <- simdata
  d$w <- 1 + d$id %% 3     # weights 1, 2, 3 by unit
  d$w2 <- 1 + d$id %% 4    # a different weight column
  d[!(d$id %in% 110:115 & d$time %in% 1:4), ]   # 24 control cells missing
})

.fit_w <- function(...) {
  suppressMessages(
    gsynth(Y ~ D + X1 + X2, data = .wdata, index = c("id", "time"),
           force = "two-way", parallel = FALSE, seed = 1, ...)
  )
}

## weighted mean of a fit's unit-period effects over the treated cells
.weighted_att <- function(fit, d, wcol) {
  W <- matrix(NA_real_, nrow(fit$Y.dat), ncol(fit$Y.dat))
  W[cbind(match(d$time, fit$rawtime), match(d$id, fit$id))] <- d[[wcol]]
  cells <- which(fit$D.dat == 1 & fit$I.dat == 1)
  sum(W[cells] * fit$eff[cells]) / sum(W[cells])
}

test_that("`weight` averages the effects of the unweighted model fit", {
  f_u <- .fit_w(r = 2, CV = FALSE, se = FALSE)
  f_weight <- .fit_w(r = 2, CV = FALSE, se = FALSE, weight = "w")
  f_wagg <- .fit_w(r = 2, CV = FALSE, se = FALSE, W.agg = "w")
  treated <- which(f_u$D.dat == 1)

  ## the model fit is the unweighted one
  expect_equal(f_weight$beta, f_u$beta)
  expect_equal(f_weight$eff[treated], f_u$eff[treated])
  ## the averages are weighted, exactly as with W.agg
  expect_equal(f_weight$att.avg, f_wagg$att.avg)
  expect_equal(f_weight$att, f_wagg$att)
  expect_equal(f_weight$att.avg, .weighted_att(f_u, .wdata, "w"))
  ## and the weights matter here
  expect_false(isTRUE(all.equal(f_weight$att.avg, f_u$att.avg)))
})

test_that("`weight` is recorded as an averaging weight, not a fit weight", {
  f_weight <- .fit_w(r = 2, CV = FALSE, se = FALSE, weight = "w")
  expect_false(f_weight$W.in.fit)
  expect_true(f_weight$W.in.agg)
})

## the fitted object, or the error condition if the fit stops
.try_fit <- function(...) tryCatch(.fit_w(...), error = function(e) e)
.error_message <- function(x) if (inherits(x, "error")) conditionMessage(x) else ""

test_that("`weight` and a different `W.agg` stop; the same column runs", {
  f_conflict <- .try_fit(r = 2, CV = FALSE, se = FALSE, weight = "w",
                         W.agg = "w2")
  expect_s3_class(f_conflict, "error")
  expect_match(.error_message(f_conflict),
               "`weight` and `W.agg` name different columns", fixed = TRUE)
  expect_match(.error_message(f_conflict), "W.est", fixed = TRUE)
  f_same <- .fit_w(r = 2, CV = FALSE, se = FALSE, weight = "w", W.agg = "w")
  f_weight <- .fit_w(r = 2, CV = FALSE, se = FALSE, weight = "w")
  expect_equal(f_same$att.avg, f_weight$att.avg)
  expect_equal(f_same$beta, f_weight$beta)
})

## Needs fect >= 2.4.6 with the cross-validation weight fix (fect run B4):
## with fect dev 412d7ae, CV fits the model with the averaging weight.
test_that("`weight` stays out of the model fit under cross-validation", {
  f_u <- .fit_w(r = c(2, 2), CV = TRUE, se = FALSE)
  f_weight <- .fit_w(r = c(2, 2), CV = TRUE, se = FALSE, weight = "w")
  expect_equal(f_weight$beta, f_u$beta)
  expect_equal(f_weight$att.avg, .weighted_att(f_u, .wdata, "w"))
})

## Needs fect >= 2.4.6 with the parametric-bootstrap weight fix (fect run
## A4): with fect dev 412d7ae, weights + parametric inference stop with
## "number of items to replace is not a multiple of replacement length".
test_that("`weight` works with se = TRUE and parametric inference", {
  f_weight <- .try_fit(r = 2, CV = FALSE, weight = "w", se = TRUE,
                       inference = "parametric", nboots = 20)
  expect_identical(.error_message(f_weight), "")
  expect_s3_class(f_weight, "gsynth")
  if (inherits(f_weight, "gsynth")) {
    expect_identical(f_weight$vartype, "parametric")
    expect_equal(unname(f_weight$est.avg[1, 1]), f_weight$att.avg)
  }
})
