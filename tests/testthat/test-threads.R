# The thread policy: the parallel kernels decompose over the elements of
# their output, so each observation's derivatives are written in full by one
# thread and no reduction is ever split. The twins therefore ask for
# identical() and not a tolerance.
#
# They asked for 1e-13 until 2026-08-21, on the reading that the last-bit
# differences the Windows CI runner reported out of its own polygamma path
# were the runtime's and unbindable by package code. That was measured
# again and the reading was wrong twice over. The differences are not
# deterministic -- gamma1's deriv3 at phi = 1/19 returned six distinct
# results over six identical calls at one thread count -- and they ARE
# bindable: a worker that installs the calling thread's floating-point
# environment reproduces the sequential value exactly (d7_par.h). What the
# tolerance had been hiding is therefore a fit whose answer moved between
# two runs, which is what identical() is here to keep out.
# Tests ask for 2 threads, which is what CRAN's check machines have.

test_that("converted kernels agree at 1 and 2 threads", {
  set.seed(1)
  n <- 40000                       # above every internal threshold
  y <- rnorm(n, 1, 2)
  mu <- rep(1, n)
  sg <- rep(2, n)
  for (fn in list(gaussian_gradient_cpp, gaussian_hessian_cpp,
                  gaussian_expected_hessian_cpp, gaussian_deriv3_cpp,
                  gaussian_deriv4_cpp)) {
    expect_identical(fn(y, mu, sg, 1L), fn(y, mu, sg, 2L))
  }
  yg <- rgamma(n, shape = 2, scale = 1.5)
  ph <- rep(0.5, n)
  for (fn in list(gamma1_gradient_cpp, gamma1_hessian_cpp,
                  gamma1_expected_hessian_cpp, gamma1_deriv3_cpp,
                  gamma1_deriv3_expected_cpp, gamma1_deriv4_cpp,
                  gamma1_deriv4_expected_cpp)) {
    expect_identical(fn(yg, rep(3, n), ph, 1L), fn(yg, rep(3, n), ph, 2L))
  }
})

test_that("the transcendental families' kernels agree too", {
  set.seed(3)
  n <- 20000                       # above kMinCostly for every family here
  both <- function(fn, ...) expect_identical(fn(..., threads = 1L),
                                            fn(..., threads = 2L))
  yp <- rpois(n, 4)
  for (fn in list(poisson_gradient_cpp, poisson_hessian_cpp,
                  poisson_expected_hessian_cpp, poisson_deriv3_cpp,
                  poisson_deriv3_expected_cpp, poisson_deriv4_cpp,
                  poisson_deriv4_expected_cpp)) {
    both(fn, as.numeric(yp), rep(4, n))
  }
  yn <- rnbinom(n, mu = 5, size = 2)
  for (fn in list(negbin_gradient_cpp, negbin_hessian_cpp,
                  negbin_expected_hessian_cpp, negbin_deriv3_cpp,
                  negbin_deriv3_expected_cpp, negbin_deriv4_cpp,
                  negbin_deriv4_expected_cpp)) {
    both(fn, as.numeric(yn), rep(5, n), rep(2, n))
  }
  for (fn in list(negbin1_gradient_cpp, negbin1_hessian_cpp,
                  negbin1_expected_hessian_cpp)) {
    both(fn, as.numeric(yn), rep(5, n), rep(2, n))
  }
  yb <- rbeta(n, 2, 3)
  for (fn in list(beta_gradient_cpp, beta_hessian_cpp,
                  beta_expected_hessian_cpp, beta_deriv3_cpp,
                  beta_deriv4_cpp)) {
    both(fn, yb, rep(0.4, n), rep(5, n))
  }
  yt <- rt(n, df = 6)
  for (fn in list(student_t_gradient_cpp, student_t_hessian_cpp,
                  student_t_expected_hessian_cpp, student_t_deriv3_cpp,
                  student_t_deriv4_cpp)) {
    both(fn, yt, rep(0, n), rep(1.2, n), rep(6, n))
  }
  yi <- 1 / rgamma(n, shape = 3, rate = 2)
  for (fn in list(invgauss_gradient_cpp, invgauss_hessian_cpp,
                  invgauss_expected_hessian_cpp, invgauss_deriv3_cpp,
                  invgauss_deriv3_expected_cpp, invgauss_deriv4_cpp,
                  invgauss_deriv4_expected_cpp)) {
    both(fn, yi, rep(1.5, n), rep(0.5, n))
  }
  ybb <- rbinom(n, 10, 0.4)
  for (fn in list(betabinom_gradient_cpp, betabinom_hessian_cpp,
                  betabinom_expected_hessian_cpp)) {
    both(fn, as.numeric(ybb), rep(0.4, n), rep(0.3, n), 10)
  }
  yg2 <- rgamma(n, shape = 2, scale = 1.5)
  for (fn in list(gamma_gradient_cpp, gamma_hessian_cpp,
                  gamma_expected_hessian_cpp, gamma_deriv3_cpp,
                  gamma_deriv3_expected_cpp, gamma_deriv4_cpp,
                  gamma_deriv4_expected_cpp)) {
    both(fn, yg2, rep(3, n), rep(4.5, n))
  }
})

test_that("a fit does not depend on threads", {
  set.seed(2)
  y <- rnorm(30000, 2, 3)
  f1 <- fit_distrib(gaussian1_distrib(), y)
  f2 <- fit_distrib(gaussian1_distrib(), y,
                    threads = numericals7::n_threads(2))
  expect_equal(coef(f1), coef(f2), tolerance = 1e-10)
  expect_equal(f1@loglik, f2@loglik, tolerance = 1e-10)
  expect_equal(vcov(f1), vcov(f2), tolerance = 1e-10)

  yg <- rgamma(6000, shape = 2, scale = 1.5)
  g1 <- fit_distrib(gamma1_distrib(), yg)
  g2 <- fit_distrib(gamma1_distrib(), yg,
                    threads = numericals7::n_threads(2))
  expect_equal(coef(g1), coef(g2), tolerance = 1e-10)
  expect_equal(g1@loglik, g2@loglik, tolerance = 1e-10)
})

test_that("threads must be the n_threads() object", {
  expect_error(fit_distrib(gaussian1_distrib(), rnorm(50), threads = 2),
               "n_threads")
})

test_that("the process-level setting is restored after a fit", {
  old <- Sys.getenv("RCPP_PARALLEL_NUM_THREADS", unset = NA_character_)
  invisible(fit_distrib(gaussian1_distrib(), rnorm(200),
                        threads = numericals7::n_threads(2)))
  expect_identical(Sys.getenv("RCPP_PARALLEL_NUM_THREADS",
                              unset = NA_character_), old)
})

# The regression for the floating-point environment of d7_par.h's worker.
# These are the two arguments the sweep of 2026-08-21 found: R's psigamma
# returns one value on the calling thread and another, one ulp away, on a
# TBB worker at x = 19 and x = 40, so a kernel reading it there gave a
# different answer according to which chunk an element landed in. gamma1
# takes its shape from 1/phi and negbin2's third derivative reads
# psigamma(y + theta, 2), which sweeps the support -- between them they
# cover the parameter-scalar and the data-dependent argument.
#
# Repeatability AT ONE COUNT is the assertion that matters: the failure was
# not that two counts disagreed but that two identical calls did, so a
# 1-against-2 comparison alone could pass by luck.
test_that("a kernel reading psigamma is reproducible and count-free", {
  set.seed(5)
  n <- 20000
  yg <- rgamma(n, shape = 19, scale = 1)
  for (phi in c(1 / 19, 1 / 40)) {                  # shape 19 and 40
    args <- list(yg, rep(2, n), rep(phi, n))
    seq3 <- do.call(gamma1_deriv3_cpp, c(args, 1L))
    for (i in 1:4) {
      expect_identical(do.call(gamma1_deriv3_cpp, c(args, 2L)), seq3)
    }
  }

  yn <- rnbinom(n, size = 2, mu = 3)
  args <- list(as.numeric(yn), rep(3, n), rep(2, n))
  seq3 <- do.call(negbin_deriv3_cpp, c(args, 1L))
  for (i in 1:4) {
    expect_identical(do.call(negbin_deriv3_cpp, c(args, 2L)), seq3)
  }
})

# `threads` says HOW MANY, not merely whether: par_for() hands the count to
# parallelFor(), whose resolveValue() prefers an explicit positive value to
# RCPP_PARALLEL_NUM_THREADS. Until 2026-08-21 the count was left to that
# variable, so a caller outside a fit -- which is every caller that does not
# go through numericals7::local_threads() -- ran on every core the machine
# had whatever it asked for. The assertion is on the ANSWER rather than on a
# timing, which is what a test can hold: what a wrong count would break is
# not the value but the promise, so the value is checked across four counts
# and the process-level setting is asserted untouched.
test_that("a kernel honours its own count and leaves the process alone", {
  set.seed(6)
  n <- 20000
  before <- Sys.getenv("RCPP_PARALLEL_NUM_THREADS", unset = NA_character_)
  yg <- rgamma(n, shape = 3, scale = 1.5)
  ref <- gamma1_gradient_cpp(yg, rep(3, n), rep(0.5, n), 1L)
  for (k in c(2L, 3L, 4L)) {
    expect_identical(gamma1_gradient_cpp(yg, rep(3, n), rep(0.5, n), k), ref)
  }
  after <- Sys.getenv("RCPP_PARALLEL_NUM_THREADS", unset = NA_character_)
  expect_identical(after, before)
})

# The families converted in 0.30.0. The arguments VARY BY OBSERVATION, which
# is the only setting where the trap this conversion had to avoid shows: a
# scalar hoisted out of the loop and written inside it is shared state once
# the iterations are split, and with a scalar parameter every thread would
# write the same value and the answer would come out right by accident. Three
# families had it (chisq, exponential, geometric) and the identity below is
# what caught them.
test_that("the newly converted kernels agree at 1 and 2 threads", {
  set.seed(7)
  n <- 40000                      # above every internal threshold
  vp <- runif(n, 0.5, 3)          # a parameter per observation
  vq <- runif(n, 0.5, 3)
  both <- function(fn, ...) expect_identical(fn(..., threads = 1L),
                                             fn(..., threads = 2L))
  yr <- rnorm(n, 1, 2); yp <- rgamma(n, 3, 1); yc <- rpois(n, 3)

  for (fn in list(chisq_gradient_cpp, chisq_hessian_cpp,
                  chisq_deriv3_cpp, chisq_deriv4_cpp)) {
    both(fn, yp, vp + 1)
  }
  for (fn in list(exponential_gradient_cpp, exponential_hessian_cpp,
                  exponential_expected_hessian_cpp, exponential_deriv3_cpp,
                  exponential_deriv4_cpp)) {
    both(fn, yp, vp)
  }
  for (fn in list(geometric_gradient_cpp, geometric_hessian_cpp,
                  geometric_expected_hessian_cpp, geometric_deriv3_cpp,
                  geometric_deriv4_cpp)) {
    both(fn, as.numeric(rgeom(n, 0.3)), vp)
  }
  for (fn in list(cauchy_gradient_cpp, cauchy_hessian_cpp,
                  cauchy_expected_hessian_cpp, cauchy_deriv3_cpp,
                  cauchy_deriv4_cpp)) {
    both(fn, yr, vp, vq)
  }
  for (fn in list(gaussian2_gradient_cpp, gaussian2_hessian_cpp,
                  gaussian2_deriv3_cpp, gaussian2_deriv4_cpp)) {
    both(fn, yr, vp, vq)
  }
  for (fn in list(gaussian3_gradient_cpp, gaussian3_hessian_cpp,
                  gaussian3_deriv3_cpp, gaussian3_deriv4_cpp)) {
    both(fn, yr, vp, vq)
  }
  for (fn in list(logistic_gradient_cpp, logistic_hessian_cpp,
                  logistic_deriv3_cpp, logistic_deriv4_cpp)) {
    both(fn, yr, vp, vq)
  }
  for (fn in list(lognormal_gradient_cpp, lognormal_hessian_cpp,
                  lognormal_deriv3_cpp, lognormal_deriv4_cpp)) {
    both(fn, yp, vp, vq)
  }
  for (fn in list(weibull_deriv3_cpp, weibull_deriv4_cpp)) both(fn, yp, vp, vq)
  for (fn in list(gumbel_deriv3_cpp, gumbel_deriv4_cpp)) both(fn, yr, vp, vq)
  for (fn in list(laplace_deriv3_cpp, laplace_deriv4_cpp)) both(fn, yr, vp, vq)
  for (fn in list(invgauss2_gradient_cpp, invgauss2_hessian_cpp)) {
    both(fn, yp, vp, vq)
  }
  # these two carry the family's own `expected` flag before the count
  for (fn in list(invgauss2_deriv3_cpp, invgauss2_deriv4_cpp)) {
    both(fn, yp, vp, vq, FALSE)
    both(fn, yp, vp, vq, TRUE)
  }
  both(bernoulli_gradient_cpp, as.numeric(rbinom(n, 1, .4)), runif(n, .2, .8))
  both(binomial_gradient_cpp, as.numeric(rbinom(n, 10, .4)), runif(n, .2, .8), 10)
  both(gengamma_gradient_cpp, yp, vp, vq, vp + 0.5)
  both(gpd_gradient_cpp, yp, vp, runif(n, 0.05, 0.4))
  both(skewnormal_deriv3_cpp, yr, vp, vq, runif(n, -1, 1))
  both(skewnormal_deriv4_cpp, yr, vp, vq, runif(n, -1, 1))
})

# The generalized Pareto's near-zero branch is a polynomial in u = xi z with
# scalar coefficients, evaluated by a kernel: the two elementwise powers the
# R loop raised per iteration are algebraically one, and what is left is a
# scalar recursion. Measured, the derivative went from 660 ms at n = 20000 to
# 17, so the parallel region is a small part of what that closed -- the test
# is that the two thread counts agree and the answer is the one the Leibniz
# branch gives where both are accurate.
test_that("the gpd series kernel is count-free and matches the other branch", {
  set.seed(8)
  n <- 20000
  y <- rexp(n)
  th <- list(sigma = 1.5, xi = 0.2)
  for (ord in 3:4) {
    a <- gpd_components(y, th, ord, threads = 1L)
    expect_identical(gpd_components(y, th, ord, threads = 2L), a)
  }
  # the two branches are compared against each other in test-gpd.R, at
  # arguments chosen so that both are accurate: the Leibniz form's accuracy
  # is eps (xi z)^-b, so an overlap reaching down to a small xi z would be
  # testing that form's cancellation and not this kernel
})
