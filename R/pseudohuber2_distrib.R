#' @include distrib.R generics.R numerical_functions.R expected_loc_scale.R expected_derivatives.R y_higher.R theta2_families.R cdf_derivatives.R cdf_derivatives_families.R cross_derivatives_families.R moments.R cdf_mapped_higher.R cdf_compiled.R
NULL

#' @title Pseudo-Huber Distribution Class, Standard-Deviation Parametrization
#' @name PseudoHuber2Distrib
#'
#' @description
#' The S7 class of the pseudo-Huber family parametrized by its location
#' \eqn{\mu}, its standard deviation \eqn{\sigma} and its shape \eqn{\nu}. It
#' is the family of [pseudohuber_distrib()] with the scale \eqn{\sigma_1} of
#' that family replaced by \eqn{\sigma = \sigma_1 \sqrt{R(\nu)}}, where
#' \eqn{R(\nu) = \sqrt{\nu}\, K_2(\sqrt{\nu}) / K_1(\sqrt{\nu})}. It inherits
#' from `continuous_distrib`.
#'
#' Build one with [pseudohuber2_distrib()], which supplies the three link
#' functions and fills the properties in. This page documents the raw S7
#' constructor, which takes the parent's properties and validates none of the
#' relationships between them.
#'
#' @inheritParams distrib
#'
#' @return An S7 object of class `PseudoHuber2Distrib`, inheriting from
#'   `continuous_distrib` and from `distrib`. For an object built by
#'   [pseudohuber2_distrib()] the properties hold `"pseudo huber2"`,
#'   `"univariate"`, `c(-Inf, Inf)`, `c("mu", "sigma", "nu")`, the
#'   interpretations `c(mu = "location", sigma = "standard deviation",
#'   nu = "shape")`, `3`, and the domains \eqn{(-\infty, \infty)},
#'   \eqn{(0, \infty)}, \eqn{(0, \infty)}.
#'
#' @section Methods:
#' Registered on this class in this file: [distrib_pdf()], [distrib_cdf()],
#' [distrib_quantile()], [distrib_rng()], [distrib_gradient()],
#' [distrib_hessian()], [distrib_deriv3()], [distrib_deriv4()],
#' [distrib_deriv5()], [distrib_grad_y()], [distrib_hess_y()],
#' [distrib_deriv3_y()], [distrib_deriv4_y()], [distrib_cross_y()],
#' [distrib_cross2_y()], [distrib_grad_y_hess()], [distrib_hess_y_hess()],
#' [distrib_grad_cdf()], [distrib_hess_cdf()], [distrib_deriv3_cdf()],
#' [distrib_deriv4_cdf()], the expected information with
#' its two derivatives through [register_loc_scale_expected()], and the four
#' moments [mean()], [variance()], [skewness()] and [kurtosis()].
#'
#' @seealso [pseudohuber2_distrib()] to build one; [PseudoHuberDistrib] for
#'   the scale parametrization.
#'
#' @examples
#' d <- pseudohuber2_distrib()
#' S7::S7_inherits(d, continuous_distrib)
#' d@params_interpretation
PseudoHuber2Distrib <- S7::new_class("PseudoHuber2Distrib",
                                     parent = continuous_distrib)


#' @title Pseudo-Huber Density, Standard-Deviation Parametrization
#' @name distrib_pdf.PseudoHuber2Distrib
#'
#' @description
#' Computes the density
#' \deqn{f(y; \mu, \sigma, \nu) = \dfrac{\sqrt{R(\nu)}}{2 \sigma \sqrt{\nu}\,
#'       K_1(\sqrt{\nu})} \exp\left(-\sqrt{\nu + R(\nu)
#'       \left(\dfrac{y-\mu}{\sigma}\right)^2}\right), \qquad
#'       R(\nu) = \dfrac{\sqrt{\nu}\, K_2(\sqrt{\nu})}{K_1(\sqrt{\nu})},}
#' whose variance is \eqn{\sigma^2}. The log-density is
#' \eqn{-D - \log 2 - \log\sigma + h(\nu)} with
#' \eqn{D = \sqrt{\nu + R(\nu)(y - \mu)^2/\sigma^2}} and
#' \eqn{h(\nu) = -\tfrac14 \log\nu + \tfrac12 \log K_2(\sqrt\nu) -
#' \tfrac32 \log K_1(\sqrt\nu)}. \eqn{R} and \eqn{h} are formed with the
#' exponentially scaled Bessel functions, so they stay finite where
#' \eqn{K_1(\sqrt\nu)} underflows.
#'
#' @param distrib A `PseudoHuber2Distrib` object, from
#'   [pseudohuber2_distrib()].
#' @param y A numeric vector of observations.
#' @param theta A named list with components `mu`, `sigma` and `nu`, each of
#'   length 1 or of the length of `y`.
#' @param log Logical of length 1. When `TRUE` the log-density is returned.
#' @param ... Unused.
#'
#' @return A numeric vector of densities, one per observation.
#'
#' @seealso [pseudohuber2_distrib()]; [distrib_pdf.PseudoHuberDistrib()] for
#'   the scale parametrization.
#'
#' @examples
#' d <- pseudohuber2_distrib()
#' th <- list(mu = 0.4, sigma = 1.5, nu = 2)
#' integrate(function(v) distrib_pdf(d, v, th), -Inf, Inf)$value
#'
#' # the same density as pseudohuber_distrib() at sigma1 = sigma / sqrt(R)
#' R <- sqrt(2) * besselK(sqrt(2), 2) / besselK(sqrt(2), 1)
#' y <- c(-2, 0, 3)
#' all.equal(distrib_pdf(d, y, th),
#'           distrib_pdf(pseudohuber_distrib(), y,
#'                       list(mu = 0.4, sigma = 1.5 / sqrt(R), nu = 2)))
S7::method(distrib_pdf, PseudoHuber2Distrib) <- function(distrib, y, theta,
                                                         log = FALSE, ...) {
  mu <- theta[[1]]
  sigma <- theta[[2]]
  nu <- theta[[3]]
  nt <- pseudohuber2_nu_terms_cpp(as.numeric(nu))
  D <- sqrt(nu + nt$R * ((y - mu) / sigma)^2)
  log_val <- -D - log(2) - log(sigma) + nt$h
  if (log) log_val else exp(log_val)
}


#' @title Pseudo-Huber Distribution Function, Standard-Deviation
#'   Parametrization
#' @name distrib_cdf.PseudoHuber2Distrib
#'
#' @description
#' Computes \eqn{P(Y \le q)} by numerical integration of the density, taken
#' by the compiled rule of [compiled_cdf()] over the tail on the side of
#' \eqn{q} away from \eqn{\mu}.
#'
#' @param distrib A `PseudoHuber2Distrib` object.
#' @param q A numeric vector of quantiles.
#' @param theta A named list with components `mu`, `sigma` and `nu`.
#' @param lower.tail Logical; if `FALSE`, \eqn{P(Y > q)} is returned.
#' @param log.p Logical; if `TRUE`, the logarithm is returned.
#' @param ... Unused.
#'
#' @return A numeric vector of probabilities.
#'
#' @seealso [distrib_quantile.PseudoHuber2Distrib()], which inverts it.
#'
#' @examples
#' d <- pseudohuber2_distrib()
#' distrib_cdf(d, c(-1, 0.4, 2), list(mu = 0.4, sigma = 1.5, nu = 2))
S7::method(distrib_cdf, PseudoHuber2Distrib) <- compiled_cdf


#' @title Pseudo-Huber Quantile Function, Standard-Deviation Parametrization
#' @name distrib_quantile.PseudoHuber2Distrib
#'
#' @description
#' Inverts [distrib_cdf.PseudoHuber2Distrib()] by root finding on the lower
#' half, in a bracket that starts ten standard deviations below \eqn{\mu} and
#' is doubled until it contains the quantile; an upper quantile is the
#' reflection of the lower one about \eqn{\mu}.
#'
#' @param distrib A `PseudoHuber2Distrib` object.
#' @param p A numeric vector of probabilities.
#' @param theta A named list with components `mu`, `sigma` and `nu`.
#' @param lower.tail Logical; if `FALSE`, `p` is an upper-tail probability.
#' @param log.p Logical; if `TRUE`, `p` is given on the log scale.
#' @param ... Unused.
#'
#' @return A numeric vector of quantiles; `NaN` for a probability outside
#'   \eqn{[0, 1]}.
#'
#' @seealso [distrib_cdf.PseudoHuber2Distrib()].
#'
#' @examples
#' d <- pseudohuber2_distrib()
#' th <- list(mu = 0.4, sigma = 1.5, nu = 2)
#' q <- distrib_quantile(d, c(0.1, 0.5, 0.9), th)
#' distrib_cdf(d, q, th)
S7::method(distrib_quantile, PseudoHuber2Distrib) <- function(distrib, p, theta,
                                                              lower.tail = TRUE,
                                                              log.p = FALSE,
                                                              ...) {
  if (log.p) p <- exp(p)
  if (!lower.tail) p <- 1 - p
  all_params <- expand_params(c(list(.p = p), theta))
  rows <- transpose_params(all_params)
  vapply(rows, function(r) {
    r <- as.list(r)
    th <- r[distrib@params]
    pi <- r$.p
    if (is.na(pi) || pi < 0 || pi > 1) return(NaN)
    if (pi == 0) return(-Inf)
    if (pi == 1) return(Inf)
    m <- th[[1]]
    if (pi == 0.5) return(m)
    reflect <- pi > 0.5
    p_low <- if (reflect) 1 - pi else pi
    lo <- m - 10 * th[[2]]
    it <- 0L
    while (distrib_cdf(distrib, lo, th) > p_low && it < 100L) {
      lo <- m - 2 * (m - lo)
      it <- it + 1L
    }
    q_low <- stats::uniroot(
      function(q) distrib_cdf(distrib, q, th) - p_low,
      lower = lo, upper = m, tol = .Machine$double.eps^0.5
    )$root
    if (reflect) 2 * m - q_low else q_low
  }, numeric(1))
}


#' @title Pseudo-Huber Random Generation, Standard-Deviation Parametrization
#' @name distrib_rng.PseudoHuber2Distrib
#'
#' @description
#' Draws `n` variates as the normal variance mixture of
#' [distrib_rng.PseudoHuberDistrib()], at the scale
#' \eqn{\sigma_1 = \sigma / \sqrt{R(\nu)}} of that family.
#'
#' @param distrib A `PseudoHuber2Distrib` object.
#' @param n The number of draws.
#' @param theta A named list with components `mu`, `sigma` and `nu`.
#' @param ... Unused.
#'
#' @return A numeric vector of `n` draws.
#'
#' @seealso [distrib_rng.PseudoHuberDistrib()].
#'
#' @examples
#' set.seed(1)
#' x <- distrib_rng(pseudohuber2_distrib(), 200,
#'                  list(mu = 0, sigma = 2, nu = 3))
#' sd(x)
S7::method(distrib_rng, PseudoHuber2Distrib) <- function(distrib, n, theta,
                                                         ...) {
  nu <- as.numeric(theta[[3]])
  s1 <- as.numeric(theta[[2]]) / sqrt(pseudohuber2_nu_terms_cpp(nu)$R)
  pseudohuber_rng_cpp(as.integer(n), as.numeric(theta[[1]]), s1, nu)
}


#' @title Pseudo-Huber Derivatives in the Parameters, Standard-Deviation
#'   Parametrization
#' @name distrib_gradient.PseudoHuber2Distrib
#' @aliases distrib_hessian.PseudoHuber2Distrib
#'   distrib_deriv3.PseudoHuber2Distrib distrib_deriv4.PseudoHuber2Distrib
#'   distrib_deriv5.PseudoHuber2Distrib
#'
#' @description
#' Return the derivatives of the log-density in \eqn{(\mu, \sigma, \nu)} of
#' orders one to five, in closed form, each order by its own compiled
#' kernel.
#'
#' @details
#' With \eqn{r = y - \mu} the log-density is
#' \eqn{\ell = -D - \log 2 - \log\sigma + h(\nu)}, \eqn{D = \sqrt{Q}},
#' \eqn{Q = \nu + R(\nu) r^2/\sigma^2}. \eqn{Q} is a product of one-variable
#' functions, so each of its partial derivatives is closed:
#' \deqn{\partial_\mu^i \partial_\sigma^j \partial_\nu^k Q =
#'       [i = j = 0,\, k = 1] + R^{(k)}(\nu)\, B_i\, C_j,}
#' with \eqn{B_0 = r^2}, \eqn{B_1 = -2r}, \eqn{B_2 = 2}, \eqn{B_i = 0} for
#' \eqn{i \ge 3}, and \eqn{C_j = (-1)^j (j+1)!\, \sigma^{-(j+2)}}. A component
#' of \eqn{D} is Faa di Bruno's sum over the set partitions of its indices,
#' with \eqn{(d/dQ)^m \sqrt{Q} = c_m Q^{1/2 - m}}. The derivatives of
#' \eqn{R} and \eqn{h} are read off \eqn{\log K_n(\sqrt\nu)}, whose
#' derivatives in \eqn{t = \sqrt\nu} follow from the ratios
#' \eqn{K_n^{(j)}(t)/K_n(t) = (-1/2)^j \sum_i \binom{j}{i} K_{|n-j+2i|}(t)/
#' K_n(t)}, formed with exponentially scaled Bessel functions. A kernel of
#' order \eqn{m} computes \eqn{R} and \eqn{h} to order \eqn{m} and the
#' partitions of order \eqn{m} only.
#'
#' The expected derivatives of orders three and four are those of
#' [expected_derivative()]; the expected information is registered through
#' [register_loc_scale_expected()].
#'
#' @param distrib A `PseudoHuber2Distrib` object.
#' @param y A numeric vector of observations.
#' @param theta A named list with components `mu`, `sigma` and `nu`.
#' @param scale `"parameter"` or `"link"`; the link scale is applied by the
#'   generic.
#' @param expected Logical; for orders three and four, whether the expected
#'   derivative is returned.
#' @param approx,nsim The approximation and Monte Carlo size of
#'   [expected_derivative()].
#' @param ... Unused.
#'
#' @return A named list with one numeric vector per component, named by the
#'   parameters joined with `_`: 3 components at order one, 6 at order two
#'   (`mu_mu`, `sigma_sigma`, `nu_nu`, `mu_sigma`, `mu_nu`, `sigma_nu`), and
#'   10, 15 and 21 at orders three to five in lexicographic order.
#'
#' @seealso [distrib_pdf.PseudoHuber2Distrib()].
#'
#' @examples
#' d <- pseudohuber2_distrib()
#' th <- list(mu = 0.4, sigma = 1.5, nu = 2)
#' y <- c(-1, 0.5, 3)
#' g <- distrib_gradient(d, y, th)
#' # against a central difference in sigma
#' h <- 1e-6
#' up <- distrib_pdf(d, y, list(mu = 0.4, sigma = 1.5 + h, nu = 2), log = TRUE)
#' dn <- distrib_pdf(d, y, list(mu = 0.4, sigma = 1.5 - h, nu = 2), log = TRUE)
#' all.equal(g$sigma, (up - dn) / (2 * h), tolerance = 1e-6)
#' names(distrib_deriv5(d, y, th))
S7::method(distrib_gradient, PseudoHuber2Distrib) <- function(
    distrib, y, theta, scale = c("parameter", "link"), ...) {
  pseudohuber2_gradient_cpp(y, theta[[1]], theta[[2]], theta[[3]])
}

S7::method(distrib_hessian, PseudoHuber2Distrib) <- function(
    distrib, y, theta, scale = c("parameter", "link"), ...) {
  pseudohuber2_hessian_cpp(y, theta[[1]], theta[[2]], theta[[3]])
}

S7::method(distrib_deriv3, PseudoHuber2Distrib) <- function(
    distrib, y, theta, expected = FALSE, scale = c("parameter", "link"),
    approx = c("integrate", "bartlett", "mc", "opg"), nsim = 10000, ...) {
  if (expected) {
    expected_derivative(distrib, y, theta, order = 3L,
                        approx = match.arg(approx), nsim = nsim)
  } else {
    pseudohuber2_deriv3_cpp(y, theta[[1]], theta[[2]], theta[[3]])
  }
}

S7::method(distrib_deriv4, PseudoHuber2Distrib) <- function(
    distrib, y, theta, expected = FALSE, scale = c("parameter", "link"),
    approx = c("integrate", "bartlett", "mc", "opg"), nsim = 10000, ...) {
  if (expected) {
    expected_derivative(distrib, y, theta, order = 4L,
                        approx = match.arg(approx), nsim = nsim)
  } else {
    pseudohuber2_deriv4_cpp(y, theta[[1]], theta[[2]], theta[[3]])
  }
}

S7::method(distrib_deriv5, PseudoHuber2Distrib) <- function(
    distrib, y, theta, scale = c("parameter", "link"), ...) {
  deriv5_scale(distrib, y, theta, pseudohuber2_deriv5_cpp(y, theta[[1]], theta[[2]], theta[[3]]),
               match.arg(scale))
}

register_loc_scale_expected(PseudoHuber2Distrib)


#' @title Pseudo-Huber Derivatives in the Response, Standard-Deviation
#'   Parametrization
#' @name distrib_grad_y.PseudoHuber2Distrib
#' @aliases distrib_hess_y.PseudoHuber2Distrib
#'   distrib_cross_y.PseudoHuber2Distrib
#'   distrib_deriv3_y.PseudoHuber2Distrib distrib_deriv4_y.PseudoHuber2Distrib
#'   distrib_cross2_y.PseudoHuber2Distrib
#'   distrib_grad_y_hess.PseudoHuber2Distrib
#'   distrib_hess_y_hess.PseudoHuber2Distrib
#'
#' @description
#' The response enters the log-density only through \eqn{r = y - \mu}, so
#' \eqn{\partial_y = -\partial_\mu}. `distrib_grad_y()` returns
#' \eqn{-R(\nu)\, r/(\sigma^2 D)} and `distrib_hess_y()` returns
#' \eqn{-R(\nu)\, \nu/(\sigma^2 D^3)}; `distrib_cross_y()` returns the
#' \eqn{\mu} row of the Hessian with its sign changed; the third and fourth
#' response derivatives are the pure-\eqn{\mu} components of the third and
#' fourth parameter derivatives, with the sign \eqn{(-1)^k}. The derivatives
#' of the response derivatives in the parameters are the location-scale ones
#' of [partial_loc_scale_grad_y_hess()] and its siblings.
#'
#' @param distrib A `PseudoHuber2Distrib` object.
#' @param y A numeric vector of observations.
#' @param theta A named list with components `mu`, `sigma` and `nu`.
#' @param scale `"parameter"` or `"link"`, for `distrib_cross_y()`.
#' @param ... Unused.
#'
#' @return A numeric vector for `distrib_grad_y()` and `distrib_hess_y()`; a
#'   named list, one component per parameter, for `distrib_cross_y()`.
#'
#' @seealso [distrib_gradient.PseudoHuber2Distrib()].
#'
#' @examples
#' d <- pseudohuber2_distrib()
#' th <- list(mu = 0.4, sigma = 1.5, nu = 2)
#' y <- c(-1, 0.5, 3)
#' all.equal(distrib_grad_y(d, y, th), -distrib_gradient(d, y, th)$mu)
S7::method(distrib_grad_y, PseudoHuber2Distrib) <- function(distrib, y, theta,
                                                            ...) {
  mu <- theta[[1]]; sigma <- theta[[2]]; nu <- theta[[3]]
  R <- pseudohuber2_nu_terms_cpp(as.numeric(nu))$R
  r <- y - mu
  D <- sqrt(nu + R * r^2 / sigma^2)
  -R * r / (sigma^2 * D)
}

S7::method(distrib_hess_y, PseudoHuber2Distrib) <- function(distrib, y, theta,
                                                            ...) {
  mu <- theta[[1]]; sigma <- theta[[2]]; nu <- theta[[3]]
  R <- pseudohuber2_nu_terms_cpp(as.numeric(nu))$R
  r <- y - mu
  D <- sqrt(nu + R * r^2 / sigma^2)
  -R * nu / (sigma^2 * D^3)
}

S7::method(distrib_cross_y, PseudoHuber2Distrib) <- function(
    distrib, y, theta, scale = c("parameter", "link"), ...) {
  h <- pseudohuber2_hessian_cpp(y, theta[[1]], theta[[2]], theta[[3]])
  n <- length(y)
  stats::setNames(list(rep_len(-h$mu_mu, n), rep_len(-h$mu_sigma, n),
                       rep_len(-h$mu_nu, n)), distrib@params)
}

S7::method(distrib_deriv3_y, PseudoHuber2Distrib) <- loc_deriv_y_k(3L)
S7::method(distrib_deriv4_y, PseudoHuber2Distrib) <- loc_deriv_y_k(4L)
S7::method(distrib_cross2_y, PseudoHuber2Distrib) <- partial_loc_scale_cross2_y
S7::method(distrib_grad_y_hess, PseudoHuber2Distrib) <-
  partial_loc_scale_grad_y_hess
S7::method(distrib_hess_y_hess, PseudoHuber2Distrib) <-
  partial_loc_scale_hess_y_hess


#' @title Pseudo-Huber Distribution-Function Derivatives,
#'   Standard-Deviation Parametrization
#' @name distrib_grad_cdf.PseudoHuber2Distrib
#' @aliases distrib_hess_cdf.PseudoHuber2Distrib
#'   distrib_deriv3_cdf.PseudoHuber2Distrib
#'   distrib_deriv4_cdf.PseudoHuber2Distrib
#'
#' @description
#' The derivatives of the distribution function in the parameters, by the
#' location-scale identities of [partial_loc_scale_grad_cdf()],
#' [partial_loc_scale_hess_cdf()] and, at orders three and four,
#' [partial_loc_scale_deriv_cdf_k()]: the location and scale components are
#' closed in the density, and the components in \eqn{\nu} are differences of
#' the quadrature.
#'
#' @param distrib A `PseudoHuber2Distrib` object.
#' @param q A numeric vector of quantiles.
#' @param theta A named list with components `mu`, `sigma` and `nu`.
#' @param lower.tail,log As in [distrib_grad_cdf()].
#' @param ... Unused.
#'
#' @return A named list, one component per parameter (first order) or per
#'   pair (second order).
#'
#' @seealso [distrib_cdf.PseudoHuber2Distrib()].
#'
#' @examples
#' d <- pseudohuber2_distrib()
#' th <- list(mu = 0.4, sigma = 1.5, nu = 2)
#' all.equal(distrib_grad_cdf(d, 1, th, log = FALSE)$mu,
#'           -distrib_pdf(d, 1, th))
S7::method(distrib_grad_cdf, PseudoHuber2Distrib) <- compiled_grad_cdf
S7::method(distrib_hess_cdf, PseudoHuber2Distrib) <- compiled_hess_cdf
S7::method(distrib_deriv3_cdf, PseudoHuber2Distrib) <-
  partial_loc_scale_deriv_cdf_k(3L)
S7::method(distrib_deriv4_cdf, PseudoHuber2Distrib) <-
  partial_loc_scale_deriv_cdf_k(4L)


#' @title Pseudo-Huber Moments, Standard-Deviation Parametrization
#' @name mean.PseudoHuber2Distrib
#' @aliases variance.PseudoHuber2Distrib skewness.PseudoHuber2Distrib
#'   kurtosis.PseudoHuber2Distrib
#'
#' @description
#' The mean is \eqn{\mu}, the variance \eqn{\sigma^2} and the skewness zero.
#' The excess kurtosis depends on \eqn{\nu} alone,
#' \eqn{3 K_3(\sqrt\nu) K_1(\sqrt\nu) / K_2(\sqrt\nu)^2 - 3}, which runs from
#' 3 at \eqn{\nu \to 0} to 0 at \eqn{\nu \to \infty}.
#'
#' @param x A `PseudoHuber2Distrib` object.
#' @param theta A named list with components `mu`, `sigma` and `nu`.
#' @param ... Unused.
#'
#' @return A numeric vector, one value per parameter row.
#'
#' @seealso [pseudohuber2_distrib()].
#'
#' @examples
#' d <- pseudohuber2_distrib()
#' th <- list(mu = 0.4, sigma = 1.5, nu = 2)
#' c(mean(d, th), variance(d, th), skewness(d, th), kurtosis(d, th))
S7::method(mean, PseudoHuber2Distrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  rep(theta[[1]], length.out = max(lengths(theta[seq_len(3)])))
}

S7::method(variance, PseudoHuber2Distrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  theta[[2]]^2 + moment_const(theta, 3L, 0)
}

S7::method(skewness, PseudoHuber2Distrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  rep(0, length.out = max(lengths(theta[seq_len(3)])))
}

S7::method(kurtosis, PseudoHuber2Distrib) <- function(x, theta, ...) {
  theta <- align_theta(x, theta)
  sq_nu <- sqrt(theta[[3]])
  k1 <- besselK(sq_nu, 1, expon.scaled = TRUE)
  k2 <- besselK(sq_nu, 2, expon.scaled = TRUE)
  k3 <- besselK(sq_nu, 3, expon.scaled = TRUE)
  3 * (k3 * k1) / (k2^2) - 3 + moment_const(theta, 3L, 0)
}


#' Pseudo-Huber Distribution, Standard-Deviation Parametrization
#'
#' @description
#' Creates a pseudo-Huber distribution object parametrized by the location
#' \eqn{\mu}, the standard deviation \eqn{\sigma} and the shape \eqn{\nu}.
#'
#' @details
#' The family is that of [pseudohuber_distrib()], whose scale \eqn{\sigma_1}
#' is replaced by the standard deviation
#' \eqn{\sigma = \sigma_1 \sqrt{R(\nu)}},
#' \eqn{R(\nu) = \sqrt\nu\, K_2(\sqrt\nu)/K_1(\sqrt\nu)}. As
#' \eqn{\nu \to \infty} the distribution tends to the gaussian with standard
#' deviation \eqn{\sigma}, so the gaussian limit is reached along the single
#' coordinate \eqn{\nu} at a fixed \eqn{\sigma}. In the scale parametrization
#' the same limit needs \eqn{\sigma_1 \to 0} together with
#' \eqn{\nu \to \infty}, and the correlation of
#' \eqn{(\log\sigma_1, \log\nu)} in the inverse expected information is
#' -0.93 at \eqn{\nu = 1} and -0.9997 at \eqn{\nu = 1000}; in this
#' parametrization the correlation of \eqn{(\log\sigma, \log\nu)} is -0.37
#' and -0.05.
#'
#' Every derivative of the log-density in the parameters to order five is
#' in closed form (see [distrib_gradient.PseudoHuber2Distrib()]), the
#' expected information is computed by one quadrature per distinct shape,
#' and the distribution function by quadrature of the density.
#'
#' @param link_mu A link for the location, by default
#'   [linkfunctions7::identity_link()].
#' @param link_sigma A link for the standard deviation, by default
#'   [linkfunctions7::log_link()].
#' @param link_nu A link for the shape, by default
#'   [linkfunctions7::log_link()].
#'
#' @return An S7 object of class `PseudoHuber2Distrib`, inheriting from
#'   `continuous_distrib`, with `distrib_name` `"pseudo huber2"`, `params`
#'   `c("mu", "sigma", "nu")` and `link_params` the three links given here.
#'
#' @references
#' Barndorff-Nielsen, O. (1978). Hyperbolic distributions and distributions on
#' hyperbolae. *Scandinavian Journal of Statistics*, **5**(3), 151-157.
#'
#' @importFrom linkfunctions7 identity_link log_link
#'
#' @examples
#' d <- pseudohuber2_distrib()
#' th <- list(mu = 0.4, sigma = 1.5, nu = 2)
#' # sigma is the standard deviation
#' c(variance = variance(d, th), sigma2 = 1.5^2)
#' # a large shape is the gaussian with standard deviation sigma
#' yy <- c(-1, 0.5, 2)
#' all.equal(distrib_pdf(d, yy, list(mu = 0.4, sigma = 1.5, nu = 1e12)),
#'           dnorm(yy, 0.4, 1.5), tolerance = 1e-6)
#'
#' @seealso [pseudohuber_distrib()] for the scale parametrization;
#'   [student_t2_distrib()] for the Student t parametrized by its standard
#'   deviation; [PseudoHuber2Distrib] for the class.
#' @export
pseudohuber2_distrib <- function(link_mu = identity_link(),
                                 link_sigma = log_link(),
                                 link_nu = log_link()) {
  PseudoHuber2Distrib(
    distrib_name = "pseudo huber2",
    dimension = "univariate",
    bounds = c(-Inf, Inf),
    params = c("mu", "sigma", "nu"),
    params_interpretation = c(mu = "location", sigma = "standard deviation",
                              nu = "shape"),
    n_params = 3,
    params_bounds = list(mu = c(-Inf, Inf), sigma = c(0, Inf),
                         nu = c(0, Inf)),
    link_params = list(mu = link_mu, sigma = link_sigma, nu = link_nu)
  )
}
