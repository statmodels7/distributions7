# Zero-inflated and zero-adjusted wrappers: probability functions, moments,
# and derivatives (finite differences + expectation-based expected Hessian).

zi_za_cases <- function() {
  list(
    zip = list(d = zero_inflated(poisson_distrib()), theta = list(mu = 3, zi = 0.25)),
    zinb = list(d = zero_inflated(negbin2_distrib()), theta = list(mu = 4, theta = 1.5, zi = 0.3)),
    zap = list(d = zero_adjusted(poisson_distrib()), theta = list(mu = 3, za = 0.4)),
    zagamma = list(d = zero_adjusted(gamma2_distrib()), theta = list(mu = 3, sigma2 = 2, za = 0.3))
  )
}

test_that("wrapped pmfs/pdfs are normalized and place the right mass at zero", {
  # ZIP: P(0) = zi + (1-zi) dpois(0, mu)
  zip <- zero_inflated(poisson_distrib())
  th <- list(mu = 3, zi = 0.25)
  expect_equal(distrib_pdf(zip, 0, th), 0.25 + 0.75 * dpois(0, 3))
  expect_equal(distrib_pdf(zip, 2, th), 0.75 * dpois(2, 3))
  # direct finite sums are the independent normalization reference: at mu = 3
  # the mass beyond k = 200 is far below the tolerance
  expect_equal(sum(distrib_pdf(zip, 0:200, th)), 1, tolerance = 1e-9)

  # ZAP (hurdle): P(0) = za exactly, positives renormalized
  zap <- zero_adjusted(poisson_distrib())
  th <- list(mu = 3, za = 0.4)
  expect_equal(distrib_pdf(zap, 0, th), 0.4)
  expect_equal(distrib_pdf(zap, 2, th), 0.6 * dpois(2, 3) / (1 - dpois(0, 3)))
  expect_equal(sum(distrib_pdf(zap, 0:200, th)), 1, tolerance = 1e-9)

  # ZA gamma: P(0) = za, continuous part scaled by (1-za)
  zag <- zero_adjusted(gamma2_distrib())
  th <- list(mu = 3, sigma2 = 2, za = 0.3)
  expect_equal(distrib_pdf(zag, 0, th), 0.3)
  expect_equal(
    stats::integrate(function(t) distrib_pdf(zag, t, th), 1e-12, Inf)$value,
    0.7,
    tolerance = 1e-6
  )
})

test_that("wrapped cdf/quantile are mutually consistent", {
  for (nm in names(zi_za_cases())) {
    case <- zi_za_cases()[[nm]]
    d <- case$d
    th <- case$theta
    p <- c(0.05, 0.3, 0.7, 0.95)
    q <- distrib_quantile(d, p, th)

    if (S7::S7_inherits(d, discrete_distrib)) {
      expect_true(all(distrib_cdf(d, q, th) >= p - 1e-12), label = paste(nm, "F(q) >= p"))
      expect_true(all(distrib_cdf(d, q - 1, th) < p), label = paste(nm, "F(q-1) < p"))
    } else {
      pos <- q > 0
      expect_equal(distrib_cdf(d, q[pos], th), p[pos], tolerance = 1e-7, label = nm)
    }
  }
})

test_that("wrapped moments match closed forms", {
  # ZIP
  zip <- zero_inflated(poisson_distrib())
  th <- list(mu = 3, zi = 0.25)
  expect_equal(mean(zip, th), 0.75 * 3, tolerance = 1e-7)
  expect_equal(variance(zip, th), 0.75 * 3 + 0.25 * 0.75 * 9, tolerance = 1e-6)

  # ZAP (hurdle): E[Y] = (1-za) mu / (1 - exp(-mu))
  zap <- zero_adjusted(poisson_distrib())
  th <- list(mu = 3, za = 0.4)
  expect_equal(mean(zap, th), 0.6 * 3 / (1 - dpois(0, 3)), tolerance = 1e-7)

  # ZA gamma (GAMLSS 9.2): E[Y] = (1-p) mu, V = (1-p) s2 + p(1-p) mu^2
  zag <- zero_adjusted(gamma2_distrib())
  th <- list(mu = 3, sigma2 = 2, za = 0.3)
  expect_equal(mean(zag, th), 0.7 * 3, tolerance = 1e-6)
  expect_equal(variance(zag, th), 0.7 * 2 + 0.3 * 0.7 * 9, tolerance = 1e-5)
})

test_that("wrapped gradients and hessians match finite differences", {
  set.seed(41)
  for (nm in names(zi_za_cases())) {
    case <- zi_za_cases()[[nm]]
    d <- case$d
    th <- case$theta
    y <- distrib_rng(d, 25, th)
    stopifnot(any(y == 0)) # make sure both branches are exercised

    a_grad <- distrib_gradient(d, y, th)
    n_grad <- fd_gradient_ref(d, y, th)
    for (p in names(th)) {
      expect_equal(a_grad[[p]], n_grad[[p]], tolerance = 1e-5,
        label = sprintf("%s grad %s", nm, p))
    }

    a_hess <- distrib_hessian(d, y, th)
    n_hess <- fd_hessian_ref(d, y, th)
    expect_setequal(names(a_hess), names(n_hess))
    for (p in names(n_hess)) {
      expect_equal(a_hess[[p]], n_hess[[p]], tolerance = 1e-4,
        label = sprintf("%s hess %s", nm, p))
    }
  }
})

test_that("wrapped expected hessians equal the expectation of the observed hessians", {
  for (nm in names(zi_za_cases())) {
    case <- zi_za_cases()[[nm]]
    d <- case$d
    th <- case$theta

    a_exp <- distrib_expected_hessian(d, 0, th)
    for (p in names(a_exp)) {
      e_num <- expectation(
        d,
        function(y, theta) distrib_hessian(d, y, theta)[[p]],
        th
      )
      expect_equal(a_exp[[p]][1], e_num, tolerance = 1e-5,
        label = sprintf("%s expected hess %s", nm, p),
        expected.label = "numerical expectation")
    }
  }
})

test_that("zero_inflated() requires a discrete parent with mass at zero", {
  # A continuous distribution has P(Y = 0) = 0: there is nothing to inflate, and
  # the model the user is after is the zero-adjusted one.
  expect_error(zero_inflated(gaussian1_distrib()), "zero_adjusted")
  expect_error(zero_inflated(gamma2_distrib()), "zero_adjusted")

  # A discrete distribution whose support starts above zero has no mass to add to.
  d <- poisson_distrib()
  d@bounds <- c(1, Inf)
  expect_error(zero_inflated(d), "requires 0 in the support")
  expect_error(zero_adjusted(d), "requires 0 in the support")
})

test_that("the wrappers refuse a support too small to identify the extra parameter", {
  expect_error(zero_inflated(bernoulli_distrib()), "not\\s+identified")
  expect_error(zero_adjusted(bernoulli_distrib()), "not\\s+identified")
  expect_error(zero_inflated(binomial_distrib(size = 1)), "not\\s+identified")
  expect_error(zero_adjusted(binomial_distrib(size = 1)), "not\\s+identified")

  # Three support points are enough for one parent parameter plus the mixture one
  expect_true(S7::S7_inherits(zero_inflated(binomial_distrib(size = 2)), ZeroInflatedDistrib))
  expect_true(S7::S7_inherits(zero_adjusted(binomial_distrib(size = 2)), ZeroAdjustedDiscreteDistrib))
})

test_that("the refusal above is not pedantry: the parameters really are lost", {
  # Built by hand, bypassing the constructor, to show what it protects against.
  zab <- ZeroAdjustedDiscreteDistrib(
    parent_distrib = bernoulli_distrib(),
    distrib_name = "hand-built zero-adjusted bernoulli",
    dimension = "univariate", bounds = c(0, 1),
    params = c("mu", "za"),
    params_interpretation = c(mu = "prob. of success", za = "prob. of zero"),
    n_params = 2,
    params_bounds = list(mu = c(0, 1), za = c(0, 1)),
    link_params = list(mu = linkfunctions7::logit_link(), za = linkfunctions7::logit_link())
  )
  # Truncating {0, 1} away from 0 leaves everything on {1}: mu cannot be seen.
  expect_equal(distrib_pdf(zab, 0:1, list(mu = 0.2, za = 0.3)),
               distrib_pdf(zab, 0:1, list(mu = 0.9, za = 0.3)))
  expect_equal(distrib_gradient(zab, c(0, 1), list(mu = 0.2, za = 0.3))$mu, c(0, 0))
})

test_that("the wrappers cannot be stacked on one another", {
  zip <- zero_inflated(poisson_distrib())
  zap <- zero_adjusted(poisson_distrib())
  zag <- zero_adjusted(gamma2_distrib())

  expect_error(zero_adjusted(zip), "already models the probability of a zero")
  expect_error(zero_inflated(zap), "already models the probability of a zero")
  expect_error(zero_inflated(zip), "already models the probability of a zero")
  expect_error(zero_adjusted(zap), "already models the probability of a zero")
  expect_error(zero_adjusted(zag), "already models the probability of a zero")

  # And, again, what that refusal is protecting: zero-truncating a zero-inflated
  # distribution cancels zi out of the likelihood, leaving a zero score.
  stacked <- ZeroAdjustedDiscreteDistrib(
    parent_distrib = zip,
    distrib_name = "hand-built hurdle over a zero-inflated poisson",
    dimension = "univariate", bounds = c(0, Inf),
    params = c("mu", "zi", "za"),
    params_interpretation = c(mu = "mean", zi = "structural zero", za = "prob. of zero"),
    n_params = 3,
    params_bounds = list(mu = c(0, Inf), zi = c(0, 1), za = c(0, 1)),
    link_params = list(mu = linkfunctions7::log_link(), zi = linkfunctions7::logit_link(),
                       za = linkfunctions7::logit_link())
  )
  expect_equal(distrib_pdf(stacked, 0:6, list(mu = 3, zi = 0.2, za = 0.4)),
               distrib_pdf(stacked, 0:6, list(mu = 3, zi = 0.7, za = 0.4)))
  expect_equal(distrib_gradient(stacked, 0:4, list(mu = 3, zi = 0.2, za = 0.4))$zi,
               rep(0, 5))
})

test_that("a continuous parent whose support misses zero is flagged", {
  shifted <- transformation(gamma2_distrib(), affine_transform(1, 2))
  expect_warning(zero_adjusted(shifted), "disconnected")
  # ... but is still built, and is a perfectly valid distribution
  d <- suppressWarnings(zero_adjusted(shifted))
  expect_equal(distrib_pdf(d, 0, list(mu = 3, sigma2 = 2, za = 0.3)), 0.3)

  # A parent that reaches zero, boundary or interior, passes without comment
  expect_silent(zero_adjusted(gamma2_distrib()))
  expect_silent(zero_adjusted(gaussian1_distrib()))
})

test_that("the wrappers carry the parent's smoothness flags", {
  # The Laplace has a kink in mu; a wrapper that forgot to say so would make
  # check_distrib() and fit_distrib() treat the observed Hessian as usable.
  expect_equal(param_smoothness(zero_adjusted(laplace_distrib())),
               c(mu = FALSE, b = TRUE, za = TRUE))
  expect_equal(param_smoothness(transformation(laplace_distrib(), affine_transform(2, 1))),
               c(mu = FALSE, b = TRUE))
  expect_equal(param_smoothness(zero_inflated(poisson_distrib())),
               c(mu = TRUE, zi = TRUE))
})

test_that("only the zero-adjusted continuous wrapper declares an atom", {
  expect_equal(distrib_atoms(gamma2_distrib(), list(mu = 2, sigma2 = 1)),
               list(y = numeric(0), p = numeric(0)))
  expect_equal(distrib_atoms(zero_adjusted(poisson_distrib()), list(mu = 2, za = 0.3)),
               list(y = numeric(0), p = numeric(0)))
  expect_equal(distrib_atoms(zero_adjusted(gamma2_distrib()), list(mu = 2, sigma2 = 1, za = 0.3)),
               list(y = 0, p = 0.3))
})

test_that("response derivatives of a mixed distribution stop at the atom", {
  zag <- zero_adjusted(gamma2_distrib())
  th <- list(mu = 3, sigma2 = 2, za = 0.3)
  y <- c(0, 0.5, 1, 4)
  parent_g <- distrib_grad_y(gamma2_distrib(), y[-1], list(mu = 3, sigma2 = 2))
  parent_h <- distrib_hess_y(gamma2_distrib(), y[-1], list(mu = 3, sigma2 = 2))

  # The (1 - za) factor does not depend on y, so away from zero nothing changes
  expect_equal(distrib_grad_y(zag, y, th), c(NaN, parent_g))
  expect_equal(distrib_hess_y(zag, y, th), c(NaN, parent_h))
})

test_that("check_distrib validates every wrapper, atoms included", {
  skip_on_cran()
  set.seed(2024)
  cases <- list(
    list(d = zero_inflated(poisson_distrib()), th = list(mu = 3, zi = 0.25)),
    list(d = zero_adjusted(poisson_distrib()), th = list(mu = 3, za = 0.4)),
    list(d = zero_adjusted(gamma2_distrib()), th = list(mu = 3, sigma2 = 2, za = 0.3)),
    list(d = zero_adjusted(gaussian1_distrib()), th = list(mu = 1, sigma = 2, za = 0.3))
  )
  for (cs in cases) {
    out <- check_distrib(cs$d, theta = cs$th, n = 40, nsim = 5e4,
                         orders = 1:2, verbose = FALSE)
    expect_equal(out$status, rep("OK", nrow(out)),
                 label = paste(cs$d@distrib_name, ":", paste(out$check[out$status == "FAIL"],
                                                             collapse = ", ")))
  }
})

test_that("atom-awareness has not blunted check_distrib", {
  skip_on_cran()
  set.seed(99)
  # A deliberately wrong gradient on the mixed distribution must still be caught.
  broken <- zero_adjusted(gamma2_distrib())
  # S7::method(gen, cls) <- fn mutates the generic in place, so the original has
  # to be put back or the rest of the suite runs against the broken one.
  old <- S7::method(distrib_gradient, ZeroAdjustedContinuousDistrib)
  on.exit(suppressMessages(
    S7::method(distrib_gradient, ZeroAdjustedContinuousDistrib) <- old
  ), add = TRUE)

  suppressMessages(
    S7::method(distrib_gradient, ZeroAdjustedContinuousDistrib) <-
      function(distrib, y, theta, scale = c("parameter", "link"), ...) {
        g <- old(distrib, y, theta)
        g[[1]] <- g[[1]] * 1.05
        g
      }
  )

  out <- check_distrib(broken, theta = list(mu = 3, sigma2 = 2, za = 0.3),
                       n = 40, nsim = 5e4, orders = 1:2, verbose = FALSE)
  expect_true("gradient vs finite differences" %in% out$check[out$status == "FAIL"])
})
