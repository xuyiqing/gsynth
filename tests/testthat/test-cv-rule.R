## gsynth() passes cv.rule, cv.nobs and cv.donut to fect (added in 1.5.0).
## Before, gsynth(..., cv.rule = "min") stopped with "unused argument".

skip_on_cran()

.rule_panel <- function(panel) {
  set.seed(panel)
  N <- 30; TT <- 16; Ntr <- 10
  F <- matrix(rnorm(TT * 2), TT, 2)
  L <- matrix(rnorm(N * 2), N, 2); L[, 2] <- 0.3 * L[, 2]
  a <- rnorm(N); xi <- rnorm(TT); X <- matrix(rnorm(N * TT), TT, N)
  T0 <- rep(Inf, N); T0[seq_len(Ntr)] <- sample(11:13, Ntr, replace = TRUE)
  D <- sapply(seq_len(N), function(i) as.integer(seq_len(TT) >= T0[i]))
  Y <- outer(xi, rep(1, N)) + outer(rep(1, TT), a) + F %*% t(L) + 0.5 * X +
    2 * D + matrix(rnorm(N * TT), TT, N)
  data.frame(id = rep(seq_len(N), each = TT), time = rep(seq_len(TT), N),
             Y = c(Y), D = c(D), X = c(X))
}

.fit_g <- function(d, ...) {
  set.seed(1)
  suppressMessages(suppressWarnings(
    gsynth(Y ~ D + X, data = d, index = c("id", "time"), force = "two-way",
           CV = TRUE, r = c(0, 4), se = FALSE, parallel = FALSE, ...)))
}

test_that("gsynth() passes cv.rule to fect", {
  skip_if(utils::packageVersion("fect") < "2.4.7",
          "fect < 2.4.7 ignores cv.rule for never-treated CV")
  d <- .rule_panel(2)
  g_min <- .fit_g(d, cv.rule = "min")
  g_1se <- .fit_g(d, cv.rule = "1se")
  r_argmin <- unname(g_min$CV.out[which.min(g_min$CV.out[, "MSPE"]), "r"])
  skip_if(r_argmin == g_1se$r.cv, "this panel no longer separates the rules")
  expect_equal(g_min$r.cv, r_argmin)
  expect_equal(.fit_g(d)$r.cv, g_1se$r.cv)  # "1se" is the default
})

test_that("gsynth() passes cv.nobs and cv.donut to block CV", {
  d <- .rule_panel(2)
  expect_no_error(.fit_g(d, cv.method = "block", cv.nobs = 5, cv.donut = 2))
  ## the settings reach fect's block folds: a different table than the defaults
  a <- .fit_g(d, cv.method = "block")
  b <- .fit_g(d, cv.method = "block", cv.nobs = 5, cv.donut = 2)
  expect_false(isTRUE(all.equal(a$CV.out, b$CV.out)))
})
