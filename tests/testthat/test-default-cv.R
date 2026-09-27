# gsynth() with CV = TRUE (the default) and one number for r chooses the
# number of factors by cross-validation from r to 5, as ?gsynth documents
# and as gsynth 1.2.1 did. From 1.3.1 to 1.4.0 the default call fitted
# r = 0 without a search (#26, #34, #102).

data("gsynth", package = "gsynth")  # loads simdata and turnout

skip_on_cran()

.fit_cv <- function(...) {
  suppressMessages(
    gsynth(Y ~ D + X1 + X2, data = simdata, index = c("id", "time"),
           force = "two-way", se = FALSE, parallel = FALSE, seed = 1, ...)
  )
}

## r values in a cross-validation table (NULL when there is no table)
.cv_r <- function(cv.out) {
  if (is.null(cv.out)) return(NULL)
  unname(cv.out[, "r"])
}

test_that("the default call cross-validates r from 0 to 5", {
  msgs <- capture_messages(
    f_def <- gsynth(Y ~ D + X1 + X2, data = simdata, index = c("id", "time"),
                    force = "two-way", se = FALSE, parallel = FALSE, seed = 1)
  )
  f_rng <- .fit_cv(r = c(0, 5))

  expect_equal(.cv_r(f_def[["CV.out"]]), c(0, 1, 2, 3, 4, 5))
  expect_identical(f_def[["CV.out"]], f_rng[["CV.out"]])
  expect_identical(f_def$r.cv, f_rng$r.cv)
  expect_identical(f_def$att.avg, f_rng$att.avg)
  ## simdata has two factors, so the search picks r > 0 and the implied
  ## weights exist (they are NULL at r = 0: #17, #82)
  expect_gt(f_def$r.cv, 0)
  expect_false(is.null(f_def$wgt.implied))
  ## no "CV cannot be performed" / single-candidate notice
  expect_false(any(grepl("cannot be performed", msgs)))
  expect_false(any(grepl("Only one candidate", msgs)))
})

test_that("one r below 5 with CV = TRUE searches r to 5", {
  f_one <- .fit_cv(r = 1, k = 5)
  f_rng <- .fit_cv(r = c(1, 5), k = 5)

  expect_equal(.cv_r(f_one[["CV.out"]]), c(1, 2, 3, 4, 5))
  expect_identical(f_one[["CV.out"]], f_rng[["CV.out"]])
  expect_identical(f_one$r.cv, f_rng$r.cv)
  expect_identical(f_one$att.avg, f_rng$att.avg)
})

test_that("r >= 5 with CV and a fixed r with CV = FALSE are unchanged", {
  ## r = 5 with CV: one candidate, the same fit as r = 5 without CV
  f_5 <- .fit_cv(r = 5)
  f_5_fixed <- .fit_cv(r = 5, CV = FALSE)
  expect_equal(.cv_r(f_5[["CV.out"]]), 5)
  expect_equal(f_5$r.cv, 5)
  expect_equal(f_5$att.avg, f_5_fixed$att.avg)

  ## CV = FALSE: r is used as given, no cross-validation table
  f_2 <- .fit_cv(r = 2, CV = FALSE)
  expect_equal(f_2$r.cv, 2)
  expect_null(f_2[["CV.out"]])

  ## a range with CV = FALSE fits the smaller value
  f_13 <- .fit_cv(r = c(1, 3), CV = FALSE)
  expect_equal(f_13$r.cv, 1)
  expect_null(f_13[["CV.out"]])
})

test_that("estimator = 'ife' without r also searches r from 0 to 5", {
  ## fit[["CV.out.ife"]]: fit$CV.out would partially match it
  f_def <- .fit_cv(estimator = "ife", k = 5)
  f_rng <- .fit_cv(estimator = "ife", r = c(0, 5), k = 5)

  expect_equal(.cv_r(f_def[["CV.out.ife"]]), c(0, 1, 2, 3, 4, 5))
  expect_identical(f_def[["CV.out.ife"]], f_rng[["CV.out.ife"]])
  expect_identical(f_def$r.cv, f_rng$r.cv)
  expect_identical(f_def$att.avg, f_rng$att.avg)
})
