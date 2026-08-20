# The thread policy: the parallel kernels decompose over the elements of
# their output, so each observation's derivatives are written in full by one
# thread and no reduction is ever split. The twins compare 1 against 2
# threads at a tolerance of 1e-13 rather than identical(): the Windows CI
# runner's R runtime returns per-thread LAST BITS from its own polygamma
# path -- one ulp, deterministic, the same value at the same index across
# three independent binaries (with the sequential branch routed through the
# worker's own function, and again with that function noinline), which is
# the opposite of a race signature and not something package code can bind.
# The tolerance still fails a split reduction or a data race by ten orders.
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
    expect_equal(fn(y, mu, sg, 1L), fn(y, mu, sg, 2L), tolerance = 1e-13)
  }
  yg <- rgamma(n, shape = 2, scale = 1.5)
  ph <- rep(0.5, n)
  for (fn in list(gamma1_gradient_cpp, gamma1_hessian_cpp,
                  gamma1_expected_hessian_cpp, gamma1_deriv3_cpp,
                  gamma1_deriv3_expected_cpp, gamma1_deriv4_cpp,
                  gamma1_deriv4_expected_cpp)) {
    expect_equal(fn(yg, rep(3, n), ph, 1L), fn(yg, rep(3, n), ph, 2L),
                 tolerance = 1e-13)
  }
})

test_that("the transcendental families' kernels agree too", {
  set.seed(3)
  n <- 20000                       # above kMinCostly for every family here
  both <- function(fn, ...) expect_equal(fn(..., threads = 1L),
                                         fn(..., threads = 2L),
                                         tolerance = 1e-13)
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
