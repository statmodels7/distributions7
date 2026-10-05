# The expected information of a family with a non-smooth parameter, through
# the generic routes. The observed l_ij has a point mass at the kink, so the
# order-2 "integrate" route takes the Bartlett form -E[l_i l_j], and the
# quadrature places a knot at every kink the family declares. References are
# mpmath quadratures at 34 digits of -E[l_i l_j], split at the kink.

user_laplace <- function(with_kink) {
  K <- S7::new_class(if (with_kink) "UserLapKink" else "UserLapNoKink",
                     parent = continuous_distrib, package = NULL)
  S7::method(distrib_pdf, K) <- function(distrib, y, theta, log = FALSE) {
    ld <- -log(2 * theta[[2]]) - abs(y - theta[[1]]) / theta[[2]]
    if (log) ld else exp(ld)
  }
  S7::method(distrib_gradient, K) <- function(distrib, y, theta,
                                              scale = c("parameter", "link"), ...) {
    r <- y - theta[[1]]
    list(mu = sign(r) / theta[[2]], sigma = (abs(r) / theta[[2]] - 1) / theta[[2]])
  }
  S7::method(distrib_hessian, K) <- function(distrib, y, theta,
                                             scale = c("parameter", "link"), ...) {
    r <- y - theta[[1]]
    s <- theta[[2]]
    list(mu_mu = 0 * r, sigma_sigma = (1 - 2 * abs(r) / s) / s^2,
         mu_sigma = -sign(r) / s^2)
  }
  if (with_kink) {
    S7::method(kink_decomposition, K) <- function(distrib) {
      kink_spec(phi = "abs", v = function(y, theta) y - theta$mu,
                dv = function(y, theta) list(mu = -1, sigma = 0),
                coef = function(theta) -1 / theta$sigma)
    }
  }
  K(distrib_name = "user laplace", dimension = "univariate", bounds = c(-Inf, Inf),
    params = c("mu", "sigma"), params_interpretation = c(mu = "loc", sigma = "scale"),
    n_params = 2, params_bounds = list(mu = c(-Inf, Inf), sigma = c(0, Inf)),
    link_params = list(mu = linkfunctions7::identity_link(),
                       sigma = linkfunctions7::log_link()),
    params_smooth = c(mu = FALSE, sigma = TRUE))
}

folded_laplace_ref <- list(
  list(th = list(mu = 0.5, sigma = 1.2),
       r = c(mu_mu = -0.27369344995195455456, sigma_sigma = -0.61090525806349373675,
             mu_sigma = -0.17531291437187078745)),
  list(th = list(mu = -1, sigma = 2),
       r = c(mu_mu = -0.11552928931500243963, sigma_sigma = -0.21050546477193256427,
             mu_sigma = 0.067235355342498780187))
)

test_that("kink_knots() locates the kink of a family, a fold and a truncation", {
  th <- list(mu = 0.5, sigma = 2)
  expect_equal(kink_knots(laplace_distrib(), th, -10, 10), 0.5)
  expect_equal(kink_knots(folded(laplace_distrib()), list(mu = -0.5, sigma = 2), 0, 10), 0.5)
  expect_equal(kink_knots(truncated(laplace_distrib(), lower = 0), th, 0, 10), 0.5)
  expect_length(kink_knots(gaussian1_distrib(), list(mu = 0, sigma = 1), -10, 10), 0L)
  # a v with two roots gives both
  hub <- kink_spec(phi = "hinge", v = function(y, theta) abs(y - theta$mu) - theta$k,
                   dv = function(y, theta) list(mu = -sign(y - theta$mu), k = -1),
                   coef = function(theta) 1)
  expect_equal(.kink_roots(hub, list(mu = 1, k = 0.3), -20, 20), c(0.7, 1.3))
})

test_that("the order-2 'integrate' route of a folded user family meets 30-digit values", {
  # Before the change the fallback integrated the observed l_mumu, which is
  # 0 almost everywhere, and returned +0.42 at the first point.
  fu <- folded(user_laplace(TRUE))
  expect_null(distrib_scalar_route(fu))
  for (x in folded_laplace_ref) {
    e <- distrib_expected_hessian(fu, 1, x$th)
    for (k in names(x$r)) expect_equal(e[[k]], x$r[[k]], tolerance = 1e-12, label = k)
  }
  # without a declared kink the form is the same and no knot is placed
  fn <- folded(user_laplace(FALSE))
  for (x in folded_laplace_ref) {
    e <- distrib_expected_hessian(fn, 1, x$th)
    for (k in names(x$r)) expect_equal(e[[k]], x$r[[k]], tolerance = 1e-7, label = k)
  }
})

test_that("the generic routes place a knot at a kink near a quantile", {
  # The median of this folded laplace2 is 1.50044, 4e-4 from the kink at 1.5;
  # without the knot mu_lambda came out -0.007947, 14 per cent off.
  d <- folded(laplace2_distrib())
  th <- list(mu = -1.5, lambda = 2.3)
  r <- c(mu_mu = -5.279348364723493958, lambda_lambda = -0.16917083628598816907,
         mu_lambda = -0.0069467186585908969318)
  for (a in c("bartlett", "integrate")) {
    e <- expected_derivative(d, 1, th, 2L, a)
    for (k in names(r)) expect_equal(e[[k]], r[[k]], tolerance = 1e-12, label = paste(a, k))
  }
  d <- folded(enet_distrib())
  th <- list(mu = -1.4, lambda = 2.3, alpha = 0.55)
  r <- c(mu_mu = -3.2381019848271992738, lambda_lambda = -0.11935186185505059692,
         alpha_alpha = -0.13021920914801414025, mu_lambda = -0.015020058765778763658,
         mu_alpha = 0.0124581790951593554, lambda_alpha = -0.080067545752604242056)
  e <- expected_derivative(d, 1, th, 2L, "integrate")
  for (k in names(r)) expect_equal(e[[k]], r[[k]], tolerance = 1e-12, label = k)
})
