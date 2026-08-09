# The elastic-net density: the product of a Laplace and a Gaussian at the
# same location, normalized. It exists so that the elastic-net penalty is
# the same construction as ridge and lasso rather than a branch of its own.

test_that("check_distrib passes on the elastic net", {
  res <- check_distrib(enet_distrib(),
                       theta = list(mu = 0, lambda = 2, alpha = 0.6),
                       verbose = FALSE)
  expect_true(all(res$status == "OK"),
              info = paste(res$check[res$status != "OK"], collapse = ", "))
})

test_that("the density integrates to one over a grid of both rates", {
  d <- enet_distrib()
  for (lam in c(0.3, 2, 20)) for (al in c(0.05, 0.4, 0.9, 0.995)) {
    th <- list(mu = 0.3, lambda = lam, alpha = al)
    f <- function(u) distrib_pdf(d, u, th)
    # split at the location: the density has a kink there and one call
    # across it is the quadrature's problem, not the density's
    v <- stats::integrate(f, -Inf, 0.3, rel.tol = 1e-12)$value +
      stats::integrate(f, 0.3, Inf, rel.tol = 1e-12)$value
    expect_equal(v, 1, tolerance = 1e-9,
                 info = sprintf("lambda %g alpha %g", lam, al))
  }
})

test_that("the two ends are the Laplace and the Gaussian", {
  d <- enet_distrib()
  yv <- c(-1, 0.2, 2)
  th <- list(mu = 0.2, lambda = 3)
  # alpha -> 1: the normalizing constant tends to 2/a, and x = a/sqrt(c)
  # reaches 1e6, where adding x^2/2 to a log-probability of the same size
  # would lose every digit
  expect_equal(distrib_pdf(d, yv, c(th, list(alpha = 1 - 1e-12))),
               distrib_pdf(laplace2_distrib(), yv, th), tolerance = 1e-10)
  expect_equal(distrib_pdf(d, yv, c(th, list(alpha = 1e-12))),
               stats::dnorm(yv, 0.2, 1 / sqrt(3)), tolerance = 1e-10)
})

test_that("cdf and quantile are inverse, and hold where the tail underflows", {
  d <- enet_distrib()
  pp <- c(1e-8, 0.01, 0.25, 0.5, 0.75, 0.99, 1 - 1e-8)
  for (lam in c(0.3, 20)) for (al in c(0.05, 0.9, 0.995)) {
    th <- list(mu = 0.3, lambda = lam, alpha = al)
    qq <- distrib_quantile(d, pp, th)
    expect_equal(distrib_cdf(d, qq, th), pp, tolerance = 1e-9,
                 info = sprintf("lambda %g alpha %g", lam, al))
    expect_true(all(is.finite(qq)))
  }
  # at lambda 20 and alpha 0.995 the argument of the normal tail is 63,
  # where pnorm(-x) is exactly zero and a ratio taken on the natural
  # scale would be NaN
  th <- list(mu = 0, lambda = 20, alpha = 0.995)
  expect_true(all(is.finite(distrib_cdf(d, c(-3, 0, 3), th))))
})

test_that("the closed variance is the second moment", {
  d <- enet_distrib()
  for (lam in c(0.3, 2)) for (al in c(0.1, 0.6)) {
    th <- list(mu = 0.3, lambda = lam, alpha = al)
    f <- function(u) (u - 0.3)^2 * distrib_pdf(d, u, th)
    num <- stats::integrate(f, -Inf, 0.3, rel.tol = 1e-12)$value +
      stats::integrate(f, 0.3, Inf, rel.tol = 1e-12)$value
    expect_equal(variance(d, th), num, tolerance = 1e-8)
  }
  expect_equal(mean(enet_distrib(), list(mu = 1.5, lambda = 1, alpha = 0.5)),
               1.5)
})

test_that("score, hessian and the response derivatives are exact", {
  skip_if_not_installed("numDeriv")
  d <- enet_distrib()
  yv <- c(-2.1, -0.4, 0.9, 3.3)
  # alpha lives in the open unit interval and numDeriv::hessian steps by
  # a tenth of the value by default, which walks straight out of it
  ma <- list(d = 1e-4, r = 4)
  for (lam in c(0.3, 2, 20)) for (al in c(0.1, 0.4, 0.8)) {
    th <- list(mu = 0.3, lambda = lam, alpha = al)
    v0 <- c(0.3, lam, al)
    ll <- function(v) sum(distrib_pdf(d, yv, list(mu = v[1], lambda = v[2],
                                                  alpha = v[3]), log = TRUE))
    g <- distrib_gradient(d, yv, th)
    expect_equal(c(sum(g$mu), sum(g$lambda), sum(g$alpha)),
                 numDeriv::grad(ll, v0, method.args = ma), tolerance = 1e-6,
                 info = sprintf("lambda %g alpha %g", lam, al))
    H <- distrib_hessian(d, yv, th)
    Ha <- matrix(c(sum(H$mu_mu), sum(H$mu_lambda), sum(H$mu_alpha),
                   sum(H$mu_lambda), sum(H$lambda_lambda), sum(H$lambda_alpha),
                   sum(H$mu_alpha), sum(H$lambda_alpha), sum(H$alpha_alpha)),
                 3, 3)
    expect_equal(Ha, numDeriv::hessian(ll, v0, method.args = ma),
                 tolerance = 1e-3, info = sprintf("lambda %g alpha %g", lam, al))
    cy <- distrib_cross_y(d, yv, th)
    ref <- numDeriv::jacobian(function(v)
      distrib_grad_y(d, yv, list(mu = v[1], lambda = v[2], alpha = v[3])),
      v0, method.args = ma)
    expect_equal(unname(cbind(cy$mu, cy$lambda, cy$alpha)), ref,
                 tolerance = 1e-6)
  }
  # the response derivatives are elementary and exact at every order
  th <- list(mu = 0.3, lambda = 2, alpha = 0.4)
  expect_equal(distrib_grad_y(d, yv, th),
               -2 * 0.4 * sign(yv - 0.3) - 2 * 0.6 * (yv - 0.3))
  expect_equal(distrib_hess_y(d, yv, th), rep(-2 * 0.6, 4))
})

test_that("a gradient that is 5 percent wrong is still caught", {
  # the guard against a validator that has been blunted. The broken
  # gradient goes on a SUBCLASS: registering a method mutates the generic
  # in place, so overwriting the real one would break it for the rest of
  # the session.
  EnetBad <- S7::new_class("EnetBad", parent = S7::S7_class(enet_distrib()),
                           package = NULL)
  S7::method(distrib_gradient, EnetBad) <- function(
      distrib, y, theta, scale = c("parameter", "link"), ...) {
    g <- distrib_gradient(enet_distrib(), y, theta)
    g$lambda <- g$lambda * 1.05
    g
  }
  good <- enet_distrib()
  bad <- EnetBad(
    distrib_name = "enet bad", dimension = "univariate",
    bounds = c(-Inf, Inf),
    params = good@params, params_interpretation = good@params_interpretation,
    n_params = 3, params_bounds = good@params_bounds,
    link_params = good@link_params, params_smooth = good@params_smooth
  )
  set.seed(5)
  res <- check_distrib(bad, theta = list(mu = 0, lambda = 2, alpha = 0.6),
                       orders = 1:2, verbose = FALSE)
  expect_true(any(grepl("gradient", res$check[res$status != "OK"])))
})
