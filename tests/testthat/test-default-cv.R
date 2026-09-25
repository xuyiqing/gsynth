# G1 (#26, #34, #102; #17/#82): with CV = TRUE (the default) and one number
# for r, gsynth() cross-validates the number of factors from r to 5, as
# ?gsynth documents and as gsynth <= 1.2.1 did. From 1.3.1 on, one r reached
# fect as a one-candidate search, so the default call fitted r = 0 (the
# two-way fixed-effects model) without any search.

data("gsynth", package = "gsynth")  # loads simdata and turnout

## Fit simdata with two-way fixed effects; return the fit and the messages
## printed while fitting.
.fit_cv <- function(...) {
  rec <- new.env()
  rec$msgs <- character(0)
  fit <- withCallingHandlers(
    gsynth(Y ~ D + X1 + X2, data = simdata, index = c("id", "time"),
           force = "two-way", se = FALSE, parallel = FALSE, seed = 1, ...),
    message = function(m) {
      rec$msgs <- c(rec$msgs, conditionMessage(m))
      invokeRestart("muffleMessage")
    }
  )
  list(fit = fit, msgs = rec$msgs)
}

test_that("the default call cross-validates r from 0 to 5, like r = c(0, 5)", {
  skip_on_cran()
  dflt <- .fit_cv()
  rng <- .fit_cv(r = c(0, 5))

  expect_equal(unname(dflt$fit[["CV.out"]][, "r"]), 0:5)
  expect_identical(dflt$fit$r.cv, rng$fit$r.cv)
  expect_identical(dflt$fit$att.avg, rng$fit$att.avg)
  expect_identical(dflt$fit[["CV.out"]], rng$fit[["CV.out"]])
  expect_identical(dflt$fit$eff, rng$fit$eff)
  ## implied weights exist whenever factors are chosen (#17/#82)
  if (dflt$fit$r.cv > 0) {
    expect_false(is.null(dflt$fit$wgt.implied))
  }
  ## no "too few pre-treatment records" or one-candidate message
  expect_false(any(grepl("cannot be performed", dflt$msgs)))
  expect_false(any(grepl("Only one candidate", dflt$msgs)))
})

test_that("one r below 5 with CV searches r to 5, like c(r, 5)", {
  skip_on_cran()
  one <- .fit_cv(r = 1)$fit
  rng <- .fit_cv(r = c(1, 5))$fit

  expect_equal(unname(one[["CV.out"]][, "r"]), 1:5)
  expect_identical(one$r.cv, rng$r.cv)
  expect_identical(one$att.avg, rng$att.avg)
  expect_identical(one[["CV.out"]], rng[["CV.out"]])
})

test_that("r >= 5 with CV, and CV = FALSE, still fit a single model", {
  skip_on_cran()
  ## r = 5 with CV: one candidate, as before
  five <- .fit_cv(r = 5)$fit
  expect_equal(nrow(five[["CV.out"]]), 1L)
  expect_equal(unname(five[["CV.out"]][, "r"]), 5)
  expect_equal(five$r.cv, 5)

  ## CV = FALSE: r is used as given, no search
  fixed <- .fit_cv(r = 2, CV = FALSE)$fit
  expect_null(fixed[["CV.out"]])
  expect_equal(fixed$r.cv, 2)

  none <- .fit_cv(CV = FALSE)$fit
  expect_null(none[["CV.out"]])
  expect_equal(none$r.cv, 0)
})

test_that("estimator = 'ife' without r also searches 0 to 5", {
  skip_on_cran()
  dflt <- .fit_cv(estimator = "ife")$fit
  rng <- .fit_cv(estimator = "ife", r = c(0, 5))$fit

  ## ife fits keep their table in CV.out.ife ($CV.out would partially match)
  expect_equal(unname(dflt[["CV.out.ife"]][, "r"]), 0:5)
  expect_identical(dflt$r.cv, rng$r.cv)
  expect_identical(dflt$att.avg, rng$att.avg)
  expect_identical(dflt[["CV.out.ife"]], rng[["CV.out.ife"]])
})
