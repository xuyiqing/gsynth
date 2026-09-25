# gsynth() stores the call as typed and records the inference it used in
# fit$vartype. Before, it copied a typed `inference` into the stored call
# as `vartype`: print() showed an argument the user never typed, update()
# failed with "unused argument (vartype = ...)", and fect's effect(), which
# read call$vartype, got NULL after a default fit and built its intervals
# from placebo draws centered at zero (#37, #60, #75).

data("gsynth", package = "gsynth")  # loads simdata and turnout

skip_on_cran()

# Shared fixtures (called directly, not through a helper, so that the
# stored calls can be re-evaluated by update())
f_def <- suppressMessages(
  gsynth(Y ~ D + X1 + X2, data = simdata, index = c("id", "time"),
         force = "two-way", r = 2, CV = FALSE, se = TRUE, nboots = 50,
         parallel = FALSE, seed = 3)
)
f_par <- suppressMessages(
  gsynth(Y ~ D + X1 + X2, data = simdata, index = c("id", "time"),
         force = "two-way", r = 2, CV = FALSE, se = TRUE, nboots = 50,
         inference = "parametric", parallel = FALSE, seed = 3)
)
f_np <- suppressMessages(
  gsynth(Y ~ D + X1 + X2, data = simdata, index = c("id", "time"),
         force = "two-way", r = 2, CV = FALSE, se = TRUE, nboots = 20,
         inference = "nonparametric", parallel = FALSE, seed = 3)
)
f_jk <- suppressMessages(
  gsynth(Y ~ D + X1 + X2, data = simdata, index = c("id", "time"),
         force = "two-way", r = 2, CV = FALSE, se = TRUE,
         inference = "jackknife", parallel = FALSE, seed = 3)
)
inf_var <- "jackknife"
f_var <- suppressMessages(
  gsynth(Y ~ D + X1 + X2, data = simdata, index = c("id", "time"),
         force = "two-way", r = 2, CV = FALSE, se = TRUE,
         inference = inf_var, parallel = FALSE, seed = 3)
)
f_nose <- suppressMessages(
  gsynth(Y ~ D + X1 + X2, data = simdata, index = c("id", "time"),
         force = "two-way", r = 2, CV = FALSE, se = FALSE, parallel = FALSE)
)
f_nose_typed <- suppressMessages(
  gsynth(Y ~ D + X1 + X2, data = simdata, index = c("id", "time"),
         force = "two-way", r = 2, CV = FALSE, se = FALSE,
         inference = "parametric", parallel = FALSE)
)

test_that("se = TRUE fits record the inference used as one string in vartype", {
  expect_identical(f_def$vartype, "parametric")   # default for the gsynth estimator
  expect_identical(f_par$vartype, "parametric")
  expect_identical(f_np$vartype, "bootstrap")     # "nonparametric" is the bootstrap
  expect_identical(f_jk$vartype, "jackknife")
  expect_identical(f_var$vartype, "jackknife")    # value of the variable, not its name
})

test_that("se = FALSE fits leave vartype NULL", {
  expect_null(f_nose$vartype)
  expect_null(f_nose_typed$vartype)
})

test_that("the stored call is the call as typed, with no vartype element", {
  for (f in list(f_def, f_par, f_np, f_jk, f_var, f_nose, f_nose_typed)) {
    expect_true(is.call(f$call))
    expect_false("vartype" %in% names(f$call))
  }
  expect_identical(f_par$call$inference, "parametric")
  expect_identical(f_np$call$inference, "nonparametric")
  expect_identical(f_var$call$inference, as.name("inf_var"))
  expect_false("inference" %in% names(f_def$call))
})

test_that("update() refits with typed, untyped and variable-passed inference", {
  refit <- function(f, ...) {
    tryCatch(suppressMessages(update(f, ...)), error = function(e) e)
  }
  u_def <- refit(f_def, nboots = 10)
  u_par <- refit(f_par, nboots = 10)
  u_var <- refit(f_var, se = FALSE)
  u_nose <- refit(f_nose_typed, r = 1)

  expect_s3_class(u_def, "gsynth")
  expect_s3_class(u_par, "gsynth")
  expect_s3_class(u_var, "gsynth")
  expect_s3_class(u_nose, "gsynth")
  expect_identical(u_par$vartype, "parametric")
  expect_identical(u_def$vartype, "parametric")
  expect_equal(u_par$att.avg, f_par$att.avg)   # same point estimate
  expect_equal(u_var$att.avg, f_var$att.avg)
  expect_equal(u_nose$r.cv, 1)
})

test_that("print() shows the call as typed and returns the fit invisibly", {
  out_par <- capture.output(res_par <- withVisible(print(f_par)))
  expect_false(any(grepl("vartype", out_par)))
  expect_true(any(grepl("inference = \"parametric\"", out_par, fixed = TRUE)))
  expect_false(res_par$visible)
  expect_identical(res_par$value, f_par)

  out_nose <- capture.output(res_nose <- withVisible(print(f_nose_typed)))
  expect_false(any(grepl("vartype", out_nose)))
  expect_false(res_nose$visible)
  expect_identical(res_nose$value, f_nose_typed)
})

## Needs fect >= 2.4.6 with the effect() fix (it reads the inference from
## fit$vartype). With fect dev 412d7ae, effect() still reads call$vartype
## and this block fails, before and after the gsynth fix.
test_that("effect() intervals cover the ATT after a default-inference fit", {
  e_def <- suppressMessages(effect(f_def, cumu = TRUE))
  m <- e_def$effect.est.att
  expect_true(is.matrix(m))
  inside <- m[, "ATT"] >= m[, "CI.lower"] & m[, "ATT"] <= m[, "CI.upper"]
  expect_true(all(inside))
  ## the same intervals as for the fit with inference typed
  e_par <- suppressMessages(effect(f_par, cumu = TRUE))
  expect_equal(e_def$effect.est.att, e_par$effect.est.att)
})
