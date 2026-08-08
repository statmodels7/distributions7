# Truncation: probability functions against closed forms, derivatives against
# finite differences, and the constructor's refusals.

trunc_cases <- function() {
  list(
    ztp   = list(d = truncated(poisson_distrib(), lower = 1), theta = list(mu = 2.5)),
    tpois = list(d = truncated(poisson_distrib(), upper = 8), theta = list(mu = 3)),
    tnb   = list(d = truncated(negbin2_distrib(), 1, 15), theta = list(mu = 4, theta = 1.5)),
    tnorm = list(d = truncated(gaussian1_distrib(), -1, 2), theta = list(mu = 0.5, sigma = 1.5)),
    tgam  = list(d = truncated(gamma2_distrib(), 0.5, 8), theta = list(mu = 3, sigma2 = 2))
  )
}

test_that("the zero-truncated Poisson matches its closed form", {
  ztp <- truncated(poisson_distrib(), lower = 1)
  th <- list(mu = 2)

  expect_equal(distrib_pdf(ztp, 0, th), 0)
  expect_equal(distrib_pdf(ztp, 1:5, th), dpois(1:5, 2) / (1 - dpois(0, 2)))
  # a direct finite sum is the independent normalization reference: at mu = 2
  # the mass beyond k = 200 is far below the tolerance
  expect_equal(sum(distrib_pdf(ztp, 1:200, th)), 1, tolerance = 1e-10)
  expect_equal(mean(ztp, th), 2 / (1 - exp(-2)), tolerance = 1e-8)

  # The support really starts at 1
  expect_equal(ztp@bounds, c(1, Inf))
  expect_true(all(distrib_rng(ztp, 500, th) >= 1))
})

test_that("the truncated Gaussian matches its closed form", {
  lo <- -1; up <- 2; mu <- 0.5; sg <- 1.5
  tn <- truncated(gaussian1_distrib(), lower = lo, upper = up)
  th <- list(mu = mu, sigma = sg)

  a <- (lo - mu) / sg
  b <- (up - mu) / sg
  Z <- pnorm(b) - pnorm(a)

  expect_equal(distrib_pdf(tn, 0.3, th), dnorm(0.3, mu, sg) / Z)
  expect_equal(distrib_pdf(tn, c(-2, 3), th), c(0, 0))
  expect_equal(stats::integrate(function(t) distrib_pdf(tn, t, th), lo, up)$value, 1,
               tolerance = 1e-9)

  # E[Y] = mu + sigma (phi(a) - phi(b)) / Z
  expect_equal(mean(tn, th), mu + sg * (dnorm(a) - dnorm(b)) / Z, tolerance = 1e-7)
  # Var(Y) = sigma^2 [1 + (a phi(a) - b phi(b))/Z - ((phi(a)-phi(b))/Z)^2]
  expect_equal(variance(tn, th),
               sg^2 * (1 + (a * dnorm(a) - b * dnorm(b)) / Z - ((dnorm(a) - dnorm(b)) / Z)^2),
               tolerance = 1e-6)
})

test_that("truncated cdf and quantile invert each other", {
  for (nm in names(trunc_cases())) {
    case <- trunc_cases()[[nm]]
    d <- case$d
    th <- case$theta
    p <- c(0.05, 0.3, 0.5, 0.7, 0.95)
    q <- distrib_quantile(d, p, th)

    expect_true(all(q >= d@bounds[1] & q <= d@bounds[2]), label = paste(nm, "in support"))
    if (S7::S7_inherits(d, discrete_distrib)) {
      expect_true(all(distrib_cdf(d, q, th) >= p - 1e-12), label = paste(nm, "F(q) >= p"))
      expect_true(all(distrib_cdf(d, q - 1, th) < p), label = paste(nm, "F(q-1) < p"))
    } else {
      expect_equal(distrib_cdf(d, q, th), p, tolerance = 1e-7, label = nm)
    }
    # Outside the truncation interval the cdf saturates
    expect_equal(distrib_cdf(d, d@bounds[1] - 1, th), 0)
  }
})

test_that("truncated gradients and hessians match finite differences", {
  set.seed(19)
  for (nm in names(trunc_cases())) {
    case <- trunc_cases()[[nm]]
    d <- case$d
    th <- case$theta
    y <- distrib_rng(d, 20, th)

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

test_that("the truncated expected hessian is the expectation of the observed one", {
  # E[H_T] = -Cov_T(score) is derived independently of the observed Hessian, so
  # agreement between the two tests both.
  for (nm in names(trunc_cases())) {
    case <- trunc_cases()[[nm]]
    d <- case$d
    th <- case$theta
    a_exp <- distrib_expected_hessian(d, 0, th)
    for (p in names(a_exp)) {
      e_num <- expectation(d, function(y, theta) distrib_hessian(d, y, theta)[[p]], th)
      expect_equal(a_exp[[p]][1], e_num, tolerance = 1e-6,
                   label = sprintf("%s expected hess %s", nm, p),
                   expected.label = "numerical expectation")
    }
  }
})

test_that("the truncated score has zero mean, as the first Bartlett identity requires", {
  for (nm in names(trunc_cases())) {
    case <- trunc_cases()[[nm]]
    d <- case$d
    th <- case$theta
    for (p in names(th)) {
      m <- expectation(d, function(y, theta) distrib_gradient(d, y, theta)[[p]], th)
      expect_equal(m, 0, tolerance = 1e-7, label = sprintf("%s E[score %s]", nm, p))
    }
  }
})

test_that("truncation refuses points that remove no mass", {
  # The motivating case: the Gamma lives on (0, Inf), so truncating at -2 does
  # nothing at all and the user meant something else.
  expect_error(truncated(gamma2_distrib(), lower = -2), "removes no probability mass")
  expect_error(truncated(gamma2_distrib(), lower = 0), "removes no probability mass")
  expect_error(truncated(poisson_distrib(), lower = 0), "removes no probability mass")
  expect_error(truncated(binomial_distrib(size = 10), upper = 10), "removes no probability mass")
  expect_error(truncated(beta1_distrib(), lower = 0.2, upper = 1.5), "removes no probability mass")

  # ... and points that would leave nothing
  expect_error(truncated(beta1_distrib(), lower = 0.2, upper = 0.1), "strictly less")
  expect_error(truncated(poisson_distrib(), 3, 3), "strictly less")
})

test_that("truncation requires at least one endpoint, and support points when discrete", {
  expect_error(truncated(poisson_distrib()), "at least one")
  expect_error(truncated(gaussian1_distrib()), "at least one")

  expect_error(truncated(poisson_distrib(), lower = 1.5), "not a point of the support")
  expect_error(truncated(negbin2_distrib(), upper = 7.2), "not a point of the support")
  # A continuous parent has no such restriction
  expect_no_error(truncated(gaussian1_distrib(), lower = 1.5))
})

test_that("a discrete truncation must leave enough support to identify the parameters", {
  # Two points carry one free probability; the negative binomial has two parameters.
  expect_error(truncated(negbin2_distrib(), 3, 4), "not\\s+be\\s+identified")
  expect_no_error(truncated(negbin2_distrib(), 3, 5))
  # One parameter needs only two points
  expect_no_error(truncated(poisson_distrib(), 3, 4))
})

test_that("truncation cannot remove zero from a distribution that models P(Y = 0)", {
  # Zero-truncating a zero-inflated law cancels zi out of the likelihood: the same
  # defect as stacking the two zero wrappers.
  expect_error(truncated(zero_inflated(poisson_distrib()), lower = 1),
               "models the\\s+probability of a zero")
  expect_error(truncated(zero_adjusted(poisson_distrib()), lower = 1),
               "models the\\s+probability of a zero")
  expect_error(truncated(zero_adjusted(gamma2_distrib()), lower = 0.5),
               "models the\\s+probability of a zero")

  # Truncating elsewhere keeps zero, and is fine
  expect_no_error(truncated(zero_inflated(poisson_distrib()), upper = 6))
  expect_no_error(truncated(zero_adjusted(gamma2_distrib()), upper = 5))
})

test_that("truncating a mixed distribution keeps its atom, rescaled", {
  zag <- zero_adjusted(gamma2_distrib())
  th <- list(mu = 3, sigma2 = 2, za = 0.3)
  tz <- truncated(zag, upper = 5)

  # The retained mass is P(Y <= 5) under the parent
  Z <- distrib_cdf(zag, 5, th)
  at <- distrib_atoms(tz, th)
  expect_equal(at$y, 0)
  expect_equal(at$p, 0.3 / Z)
  expect_equal(distrib_pdf(tz, 0, th), 0.3 / Z)

  # ... and the density plus the atom still make one
  total <- at$p + stats::integrate(function(t) distrib_pdf(tz, t, th), 1e-12, 5)$value
  expect_equal(total, 1, tolerance = 1e-7)
})

test_that("nested truncation collapses to the intersection", {
  d <- truncated(truncated(gaussian1_distrib(), -1, 3), upper = 2)
  expect_equal(c(d@lower, d@upper), c(-1, 2))
  # one level of wrapping, not two
  expect_true(S7::S7_inherits(d@parent_distrib, Gaussian1Distrib))

  d2 <- truncated(truncated(poisson_distrib(), lower = 2), upper = 9)
  expect_equal(c(d2@lower, d2@upper), c(2, 9))
  expect_true(S7::S7_inherits(d2@parent_distrib, PoissonDistrib))
})

test_that("truncation carries the parent's parameters and smoothness unchanged", {
  d <- truncated(laplace_distrib(), lower = -2, upper = 4)
  expect_equal(d@params, c("mu", "sigma"))
  expect_equal(d@n_params, 2)
  expect_equal(param_smoothness(d), c(mu = FALSE, sigma = TRUE))
  # truncation adds no parameter
  expect_equal(d@params, laplace_distrib()@params)
})

test_that("response derivatives stop at the truncation points", {
  tg <- truncated(gamma2_distrib(), 0.5, 8)
  th <- list(mu = 3, sigma2 = 2)
  y <- c(0.2, 1, 4, 9)   # 0.2 and 9 both fall outside [0.5, 8]
  parent <- distrib_grad_y(gamma2_distrib(), c(1, 4), th)

  expect_equal(distrib_grad_y(tg, y, th), c(NaN, parent, NaN))
  expect_equal(distrib_hess_y(tg, y, th),
               c(NaN, distrib_hess_y(gamma2_distrib(), c(1, 4), th), NaN))
  expect_true(is.nan(distrib_hess_y(tg, 0.2, th)))
})

test_that("truncated derivatives are correct with vectorized parameters", {
  d <- truncated(gaussian1_distrib(), -1, 2)
  mu <- c(0, 0.5, 1)
  th <- list(mu = mu, sigma = 1.5)
  y <- c(0.1, 0.4, 1.1)

  a <- distrib_gradient(d, y, th)
  # against the same computation one observation at a time
  for (i in 1:3) {
    one <- distrib_gradient(d, y[i], list(mu = mu[i], sigma = 1.5))
    expect_equal(a$mu[i], one$mu, tolerance = 1e-10)
    expect_equal(a$sigma[i], one$sigma, tolerance = 1e-10)
  }
})

test_that("check_distrib validates the truncated wrappers", {
  skip_on_cran()
  set.seed(404)
  cases <- list(
    list(d = truncated(poisson_distrib(), lower = 1), th = list(mu = 2.5)),
    list(d = truncated(gaussian1_distrib(), -1, 2), th = list(mu = 0.5, sigma = 1.5)),
    list(d = truncated(gamma2_distrib(), 0.5, 8), th = list(mu = 3, sigma2 = 2)),
    list(d = truncated(zero_adjusted(gamma2_distrib()), upper = 5),
         th = list(mu = 3, sigma2 = 2, za = 0.3))
  )
  for (cs in cases) {
    out <- check_distrib(cs$d, theta = cs$th, n = 40, nsim = 3e4,
                         orders = 1:2, verbose = FALSE)
    expect_equal(out$status, rep("OK", nrow(out)),
                 label = paste(cs$d@distrib_name, ":",
                               paste(out$check[out$status == "FAIL"], collapse = ", ")))
  }
})
