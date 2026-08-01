test_that("weight is passed through to the estimator", {
  skip_on_cran()

  # Treatment effects differ across treated units; the weights favour the
  # large-effect half heavily. A weighted ATT must therefore differ from the
  # unweighted one. Before `weight` was forwarded these were bit-identical.
  set.seed(42)
  N <- 30; T <- 16; T0 <- 12; ntr <- 10

  unit <- rep(seq_len(N), each = T)
  time <- rep(seq_len(T), times = N)
  f1 <- rnorm(T); l1 <- rnorm(N)

  treated <- seq_len(ntr)
  D <- as.integer(unit %in% treated & time > T0)

  eff <- rep(0, N)
  eff[treated[1:(ntr / 2)]] <- 1
  eff[treated[(ntr / 2 + 1):ntr]] <- 9

  Y <- l1[unit] * f1[time] + eff[unit] * D + rnorm(N * T, sd = 0.2)
  wt <- ifelse(unit %in% treated[(ntr / 2 + 1):ntr], 50, 1)
  dat <- data.frame(unit, time, Y, D, wt)

  args <- list(formula = Y ~ D, data = dat, index = c("unit", "time"),
               force = "two-way", r = 1, CV = FALSE, se = FALSE,
               parallel = FALSE)

  plain <- suppressWarnings(suppressMessages(do.call(gsynth, args)))
  wtd <- suppressWarnings(suppressMessages(
    do.call(gsynth, c(args, list(weight = "wt")))))

  expect_false(isTRUE(all.equal(plain$att.avg, wtd$att.avg)))

  # and it should move toward the weighted truth, not merely differ
  truth_unwt <- mean(eff[treated])
  truth_wtd <- weighted.mean(eff[treated], ifelse(seq_len(ntr) > ntr / 2, 50, 1))
  expect_lt(abs(wtd$att.avg - truth_wtd), abs(wtd$att.avg - truth_unwt))
})
