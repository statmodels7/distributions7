# The closed-form cdf derivatives added beyond the location-scale families,
# and the chain that carries a parent's onto a reparametrization.
#
# Every check runs against the partial-expectation integral
#   d^I F(q) = int_{lower}^{q} d^I f,
# assembled from the score and the density. That route shares no code with any
# of the implementations under test: they differentiate the distribution
# function, it integrates the derivative of the density. Checking a closed form
# against finite differences of the same cdf would compare an implementation
# with a cruder version of itself.

# the integrand is set to zero where the score overflows against a density that
# has already underflowed -- the Gumbel's left tail meets Inf * 0 there
cdf_ref1 <- function(d, q, th, p) {
  lo <- d@bounds[1]
  vapply(q, function(qq) {
    stats::integrate(function(y) {
      v <- distrib_gradient(d, y, th)[[p]] * distrib_pdf(d, y, th)
      v[!is.finite(v)] <- 0
      v
    }, lo, qq, rel.tol = 1e-10, subdivisions = 2000L)$value
  }, numeric(1))
}

cdf_ref2 <- function(d, q, th, nm) {
  lo <- d@bounds[1]
  pr <- distributions7:::hess_pairs(d@params)[[nm]]
  vapply(q, function(qq) {
    stats::integrate(function(y) {
      h <- distrib_hessian(d, y, th)[[nm]]
      g1 <- distrib_gradient(d, y, th)[[pr[1]]]
      g2 <- distrib_gradient(d, y, th)[[pr[2]]]
      v <- (h + g1 * g2) * distrib_pdf(d, y, th)
      v[!is.finite(v)] <- 0
      v
    }, lo, qq, rel.tol = 1e-10, subdivisions = 2000L)$value
  }, numeric(1))
}

expect_cdf_deriv <- function(d, th, q, order, tol, which = NULL) {
  got <- if (order == 1L) {
    distrib_grad_cdf(d, q, th, lower.tail = TRUE, log = FALSE)
  } else {
    distrib_hess_cdf(d, q, th, lower.tail = TRUE, log = FALSE)
  }
  nms <- if (is.null(which)) names(got) else which
  for (nm in nms) {
    r <- if (order == 1L) cdf_ref1(d, q, th, nm) else cdf_ref2(d, q, th, nm)
    expect_equal(got[[nm]], r, tolerance = tol,
                 label = sprintf("%s component %s at order %d",
                                 d@distrib_name, nm, order))
  }
}


test_that("the location-scale closed forms extend to the gumbel", {
  d <- gumbel_distrib()
  th <- list(mu = 0.7, sigma = 1.3)
  q <- c(-0.6, 0.4, 1.1, 2.5)
  expect_cdf_deriv(d, th, q, 1L, 1e-12)
  expect_cdf_deriv(d, th, q, 2L, 1e-11)
})


test_that("the positive families with an elementary cdf are closed at both orders", {
  d <- exponential_distrib()
  th <- list(mu = 1.6)
  q <- c(0.2, 1.1, 4)
  expect_cdf_deriv(d, th, q, 1L, 1e-12)
  expect_cdf_deriv(d, th, q, 2L, 1e-12)

  d <- weibull1_distrib()
  for (th in list(list(mu = 1.4, sigma = 1.9), list(mu = 0.6, sigma = 0.8))) {
    expect_cdf_deriv(d, th, q, 1L, 1e-11)
    expect_cdf_deriv(d, th, q, 2L, 1e-10)
  }
})


test_that("the generalized Pareto is closed across the shape, the small-shape series included", {
  d <- gpd_distrib()
  q <- c(0.2, 1.1, 2.5)
  # the ordinary case, the negative shape with its finite endpoint, and the
  # exponential limit where log(t)/xi^2 and z/(xi t) cancel to leading order
  expect_cdf_deriv(d, list(sigma = 1.2, xi = 0.35), q, 1L, 1e-11)
  expect_cdf_deriv(d, list(sigma = 2.0, xi = -0.25), q, 1L, 1e-11)
  expect_cdf_deriv(d, list(sigma = 1.2, xi = 1e-6), q, 1L, 1e-8)

  # the series and the general expression agree where both are computable,
  # which is what pins the switch: at 1e-3 the cancellation has not yet bitten
  th <- list(sigma = 1.2, xi = 1e-3)
  got <- distrib_grad_cdf(d, q, th, lower.tail = TRUE, log = FALSE)
  z <- q / 1.2
  t <- 1 + 1e-3 * z
  S <- 1 - distrib_cdf(d, q, th)
  direct <- -S * (log(t) / 1e-6 - z / (1e-3 * t))
  expect_equal(got$xi, direct, tolerance = 1e-6)
})


test_that("the lognormal is closed at second order", {
  d <- lognormal1_distrib()
  th <- list(mu = 0.2, sigma2 = 0.6)
  q <- c(0.3, 1.1, 3.2)
  expect_cdf_deriv(d, th, q, 2L, 1e-11)
})


test_that("a family location-scale in only two parameters closes that block", {
  q <- c(-0.4, 0.5, 1.8)
  d <- skewnormal1_distrib()
  th <- list(mu = 0.3, sigma = 1.2, alpha = 1.5)
  # the closed components; alpha is differenced and is not asserted here
  expect_cdf_deriv(d, th, q, 1L, 1e-10, which = c("mu", "sigma"))
  expect_cdf_deriv(d, th, q, 2L, 1e-9,
                   which = c("mu_mu", "sigma_sigma", "mu_sigma"))

  d <- student_t1_distrib()
  th <- list(mu = 0.1, sigma = 1.1, nu = 6)
  expect_cdf_deriv(d, th, q, 2L, 1e-10,
                   which = c("mu_mu", "sigma_sigma", "mu_sigma"))
})


test_that("the hand-written second parametrizations chain the parent's closed forms", {
  q <- c(-0.3, 0.6, 2.1)
  expect_cdf_deriv(gaussian2_distrib(), list(mu = 0.5, sigma2 = 1.7), q, 1L, 1e-12)
  expect_cdf_deriv(gaussian2_distrib(), list(mu = 0.5, sigma2 = 1.7), q, 2L, 1e-11)
  expect_cdf_deriv(gaussian3_distrib(), list(mu = 0.5, tau = 0.8), q, 1L, 1e-12)
  expect_cdf_deriv(gaussian3_distrib(), list(mu = 0.5, tau = 0.8), q, 2L, 1e-11)
  expect_cdf_deriv(invgauss2_distrib(), list(mu = 1.3, lambda = 2.2),
                   c(0.4, 1.2, 3), 1L, 1e-10)
})


test_that("reparametrize() carries the parent's cdf derivatives by the chain rule", {
  q <- c(0.4, 1.2, 3)
  expect_cdf_deriv(lognormal2_distrib(), list(mean = 1.4, var = 0.7), q, 1L, 1e-11)
  expect_cdf_deriv(weibull3_distrib(), list(mean = 1.5, sigma = 1.8), q, 1L, 1e-11)
  # the t's location and scale close; nu is differenced in the parent already
  expect_cdf_deriv(student_t2_distrib(), list(mu = 0.2, sigma = 1.2, nu = 7),
                   c(-0.4, 0.5, 1.8), 1L, 1e-10, which = c("mu", "sigma"))
})


test_that("the chain refuses to pretend when the parent has no closed form", {
  # gengamma1 differences its own cdf, so gengamma2's chain must fall back
  # rather than carry an approximation and call it exact
  expect_false(distributions7:::has_exact_cdf_deriv(gengamma1_distrib(), 1L))
  d <- gengamma2_distrib()
  th <- list(mean = 1.2, d = 2, p = 1.5)
  q <- c(0.5, 1.4)
  got <- distrib_grad_cdf(d, q, th, lower.tail = TRUE, log = FALSE)
  expect_equal(got$mean, cdf_ref1(d, q, th, "mean"), tolerance = 1e-6)
})


test_that("a wrong map derivative is caught by the chain", {
  # the guard on the guard: with the sqrt map's derivative 5% wrong the
  # gaussian2 gradient must disagree with the integral it is checked against
  q <- c(-0.3, 0.6, 2.1)
  th <- list(mu = 0.5, sigma2 = 1.7)
  parent <- gaussian1_distrib()
  bad <- distributions7:::md_gaussian2(th)
  bad[[2]][["2"]] <- bad[[2]][["2"]] * 1.05
  d1 <- distributions7:::chain_cdf_deriv(
    parent, q, list(mu = th[[1]], sigma = sqrt(th[[2]])), bad,
    c("mu", "sigma2"), 1L
  )
  r <- cdf_ref1(gaussian2_distrib(), q, th, "sigma2")
  expect_gt(max(abs(d1$sigma2 - r) / pmax(1e-8, abs(r))), 1e-3)
})
