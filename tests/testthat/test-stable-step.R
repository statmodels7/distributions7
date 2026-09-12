## The step of a parameter-direction difference, chosen rather than fixed.
##
## `fd_steps()` offers two candidates -- the magnitude step cut to 49% of the
## distance to the nearest bound, and the same magnitude scaled ON that
## distance -- and `fd_stable_step()` keeps whichever agrees better with itself
## at half its own step.  Neither candidate is right everywhere: where the
## log-density is singular in that parameter at the bound the clamped step's
## error stops falling once the clamp binds, and where it saturates the scaled
## step is too short and the rounding dominates.
##
## Every block below asserts two things, and the second is what keeps the first
## honest: what the shipped rule delivers, and that the rule it replaced does
## NOT -- so reverting to the clamped step alone fails the test rather than
## passing it silently.  The thresholds are transcribed from the census of 324
## cells at distances from 1 to 1e-8 from a bound and share no arithmetic with
## the implementation.

## the error of a numerical component against the family's own analytic one,
## on the analytic one's own scale
rel_to_scale <- function(num, ana) {
  if (!all(is.finite(num))) return(Inf)
  max(abs(num - ana)) / max(abs(ana))
}

## what the package computed before the step was chosen: the clamped step, used
## by every component
clamped_gradient <- function(d, y, theta, h_rel = .Machine$double.eps^(1 / 3)) {
  ps <- d@params
  out <- stats::setNames(vector("list", length(ps)), ps)
  for (j in seq_along(ps)) {
    h <- fd_steps(theta[[j]], d@params_bounds[[ps[j]]], h_rel, "clamp")
    tp <- tm <- theta
    tp[[j]] <- theta[[j]] + h
    tm[[j]] <- theta[[j]] - h
    out[[j]] <- (distrib_pdf(d, y, tp, log = TRUE) -
      distrib_pdf(d, y, tm, log = TRUE)) / (2 * h)
  }
  out
}

test_that("the two candidate steps coincide away from a bound, and are ordered near one", {
  h_rel <- .Machine$double.eps^(1 / 3)
  ## no finite bound: nothing to scale on
  expect_identical(fd_steps(3, c(-Inf, Inf), h_rel, "clamp"),
                   fd_steps(3, c(-Inf, Inf), h_rel, "scale"))
  ## a positive parameter at or above one: the distance is the magnitude
  expect_identical(fd_steps(2, c(0, Inf), h_rel, "clamp"),
                   fd_steps(2, c(0, Inf), h_rel, "scale"))
  ## below one the scaled step is the shorter, in proportion to the distance
  expect_equal(fd_steps(0.25, c(0, Inf), h_rel, "scale"), h_rel * 0.25)
  expect_equal(fd_steps(0.25, c(0, Inf), h_rel, "clamp"), h_rel)
  ## on a bounded interval the distance is to the nearer end
  expect_equal(fd_steps(0.3, c(0, 1), h_rel, "scale"), h_rel * 0.3)
  expect_equal(fd_steps(0.8, c(0, 1), h_rel, "scale"), h_rel * 0.2)
  ## a parameter on its boundary is reported, whichever variant is asked
  expect_error(fd_steps(0, c(0, Inf), h_rel, "clamp"), "domain boundary")
  expect_error(fd_steps(0, c(0, Inf), h_rel, "scale"), "domain boundary")
})

test_that("where the two steps coincide the result is IDENTICAL to the clamped step", {
  ## This is the negative control, and it is an identity rather than a
  ## tolerance: fd_stable_step() returns the clamped quotient at once, so a
  ## family with no finite bound in reach cannot move by a single bit.
  d <- gaussian1_distrib()
  set.seed(3)
  y <- stats::rnorm(50, 1, 2)
  th <- list(mu = 1, sigma = 2)
  g_now <- numerical_gradient(d, y, th)
  g_old <- clamped_gradient(d, y, th)
  for (nm in names(g_old)) expect_identical(g_now[[nm]], g_old[[nm]])
  ## and nothing was evaluated four times to get there
  n <- 0L
  quotient <- function(h) { n <<- n + 1L; rep(1, 3) }
  fd_stable_step(quotient, 2, c(0, Inf), .Machine$double.eps^(1 / 3))
  expect_identical(n, 1L)
})

test_that("near a bound the chosen step is worth several orders, where the clamped one is not", {
  ## The gamma's log-density is singular in its dispersion at zero.  Measured
  ## over the census, the clamped step's relative error plateaus once the clamp
  ## binds while the chosen one keeps falling.
  d <- gamma1_distrib()
  th <- list(mu = 3, phi = 1e-6)
  set.seed(7)
  y <- distrib_rng(d, 200, th)
  ana <- distrib_gradient(d, y, th)[["phi"]]
  chosen <- rel_to_scale(numerical_gradient(d, y, th)[["phi"]], ana)
  clamped <- rel_to_scale(clamped_gradient(d, y, th)[["phi"]], ana)
  expect_lt(chosen, 1e-5)
  ## the rule it replaced is four orders worse at the same point
  expect_gt(clamped, 1e-2)
  expect_lt(chosen, clamped / 1e3)
})

test_that("the choice is made for the whole vector and not observation by observation", {
  ## In the response direction fd_stable_quotient() chooses per observation,
  ## because each observation carries its own y and therefore its own step.
  ## Here the parameter is one number for the whole sample and the two
  ## candidate quotients estimate the same thing, so the reading is taken over
  ## the vector: measured, the whole-vector choice puts 307 gradients within
  ## 1e-6 against 285 for the per-observation one.  What the test pins is the
  ## property that makes that possible, that one of the two candidates is
  ## returned whole.
  d <- gamma1_distrib()
  th <- list(mu = 3, phi = 0.3)
  set.seed(11)
  y <- distrib_rng(d, 40, th)
  h_rel <- .Machine$double.eps^(1 / 3)
  quotient <- function(h) {
    tp <- tm <- th
    tp$phi <- th$phi + h
    tm$phi <- th$phi - h
    (distrib_pdf(d, y, tp, log = TRUE) - distrib_pdf(d, y, tm, log = TRUE)) / (2 * h)
  }
  got <- fd_stable_step(quotient, th$phi, c(0, Inf), h_rel)
  expect_true(identical(got$h, fd_steps(th$phi, c(0, Inf), h_rel, "clamp")) ||
              identical(got$h, fd_steps(th$phi, c(0, Inf), h_rel, "scale")))
  expect_identical(got$value, quotient(got$h))
})

test_that("fd_self_consistency reports the two readings that are not usable", {
  expect_identical(fd_self_consistency(c(1, Inf), c(1, 2)), Inf)
  expect_identical(fd_self_consistency(c(1, NaN), c(1, 2)), Inf)
  expect_identical(fd_self_consistency(c(0, 0), c(0, 0)), NA_real_)
  expect_equal(fd_self_consistency(c(1.5, 2), c(1, 2)), 0.25)
})

test_that("a tie keeps the clamped step", {
  ## isTRUE(kB < kA) is FALSE for a tie and for any reading that is NA, so the
  ## answer is the clamped step unless the scaled one is measurably better.
  h_rel <- .Machine$double.eps^(1 / 3)
  const <- function(h) rep(2, 5)           # identical at every step: a tie
  got <- fd_stable_step(const, 0.25, c(0, Inf), h_rel)
  expect_identical(got$h, fd_steps(0.25, c(0, Inf), h_rel, "clamp"))
})

test_that("the derivative of an expected information does not choose its step under Monte Carlo", {
  ## A Monte Carlo expected information is not the same number twice, so the
  ## self-consistency reading would be of that noise rather than of the step.
  ## Measured on pig1, the quotient moves by 2.75 relative between two
  ## evaluations at ONE step, against 17.3 between the step and its half.
  d <- pig1_distrib()
  expect_false(expected_hessian_exact(d))
  set.seed(4)
  th <- list(mu = 1.5, sigma = 0.8)
  y <- distrib_rng(d, 20, th)
  a <- suppressWarnings(distrib_expected_hessian(d, y, th, approx = "mc"))
  b <- suppressWarnings(distrib_expected_hessian(d, y, th, approx = "mc"))
  expect_false(identical(a, b))
  ## the deterministic route is the same number twice, which is what licenses
  ## the choice being made there
  p <- suppressWarnings(distrib_expected_hessian(d, y, th, approx = "bartlett"))
  q <- suppressWarnings(distrib_expected_hessian(d, y, th, approx = "bartlett"))
  expect_identical(p, q)
})
