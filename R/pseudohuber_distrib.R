#' @include distrib.R generics.R numerical_functions.R
NULL

#' @title Pseudo-Huber Distribution Class
#' @name PseudoHuberDistrib
#'
#' @description
#' The S7 class of the pseudo-Huber family, a location-scale family on the
#' whole real line whose log-density is the negative pseudo-Huber loss
#' \eqn{-\sqrt{\nu + z^2}}, with a shape \eqn{\nu > 0} interpolating between a
#' Laplace at \eqn{\nu \to 0} and a Gaussian at \eqn{\nu \to \infty}. It is the
#' symmetric hyperbolic distribution, a special case of the generalized
#' hyperbolic. It inherits from `continuous_distrib`; the eleven methods listed
#' below are registered on it in this file.
#'
#' Build one with [pseudohuber_distrib()], which supplies the three link
#' functions and fills the properties in. This page documents the raw S7
#' constructor, which takes the parent's properties and validates none of the
#' relationships between them.
#'
#' @inheritParams distrib
#'
#' @return An S7 object of class `PseudoHuberDistrib`, inheriting from
#'   `continuous_distrib` and from `distrib`. Its properties are the parent's:
#'   `distrib_name`, `dimension`, `bounds`, `params`, `params_interpretation`,
#'   `n_params`, `params_bounds`, `link_params` and `params_smooth`. For an
#'   object built by [pseudohuber_distrib()] they hold `"pseudo huber"`,
#'   `"univariate"`, `c(-Inf, Inf)`, `c("mu", "sigma", "nu")`, the
#'   interpretations `c(mu = "location", sigma = "scale", nu = "shape")`, `3`,
#'   and the domains \eqn{(-\infty, \infty)}, \eqn{(0, \infty)},
#'   \eqn{(0, \infty)}.
#'
#' @seealso [pseudohuber_distrib()] to build one;
#'   [laplace_distrib()] and [gaussian1_distrib()] for the two limits;
#'   [student_t1_distrib()] for the other robust three-parameter family here;
#'   [distrib_pdf.PseudoHuberDistrib()] for the density.
#'
#' @section Methods:
#' Registered on this class in this file:
#'   [`distrib_cdf()`][distrib_cdf.PseudoHuberDistrib],
#'   [`distrib_deriv3()`][distrib_deriv3.PseudoHuberDistrib],
#'   [`distrib_deriv4()`][distrib_deriv4.PseudoHuberDistrib],
#'   [`distrib_expected_hessian()`][distrib_expected_hessian.PseudoHuberDistrib],
#'   [`distrib_grad_y()`][distrib_grad_y.PseudoHuberDistrib],
#'   [`distrib_gradient()`][distrib_gradient.PseudoHuberDistrib],
#'   [`distrib_hess_y()`][distrib_hess_y.PseudoHuberDistrib],
#'   [`distrib_hessian()`][distrib_hessian.PseudoHuberDistrib],
#'   [`distrib_pdf()`][distrib_pdf.PseudoHuberDistrib],
#'   [`distrib_quantile()`][distrib_quantile.PseudoHuberDistrib],
#'   [`distrib_rng()`][distrib_rng.PseudoHuberDistrib], and the predicate
#'   [`expected_hessian_exact()`][expected_hessian_exact.PseudoHuberDistrib],
#'   which answers `FALSE` here.
#'
#' Registered from other files: the mixed derivative
#'   [`distrib_cross_y()`][distrib_cross_y.PseudoHuberDistrib] and the
#'   distribution-function derivatives
#'   [`distrib_grad_cdf()`][distrib_grad_cdf.PseudoHuberDistrib] and
#'   [`distrib_hess_cdf()`][distrib_hess_cdf.PseudoHuberDistrib], plus the four
#'   moments [`mean()`][mean.PseudoHuberDistrib],
#'   [`variance()`][variance.PseudoHuberDistrib],
#'   [`skewness()`][skewness.PseudoHuberDistrib] and
#'   [`kurtosis()`][kurtosis.PseudoHuberDistrib] in `moments.R`.
#'
#' Everything else is inherited from [continuous_distrib()].
#'
#' @examples
#' d <- pseudohuber_distrib()
#' S7::S7_inherits(d, continuous_distrib)
#'
#' # The properties a consumer reads to drive the family without knowing it.
#' d@params
#' d@params_interpretation
#'
#' # This is the one family here whose expected information is not written
#' # out, and the predicate says so.
#' expected_hessian_exact(d)
PseudoHuberDistrib <- S7::new_class("PseudoHuberDistrib", parent = continuous_distrib)

# --- S7 METHODS IMPLEMENTATION ---

#' @title Pseudo-Huber Probability Density Function
#' @name distrib_pdf.PseudoHuberDistrib
#' @description
#' Computes the pseudo-Huber density
#' \deqn{f(y; \mu, \sigma, \nu) = \dfrac{1}{2 \sigma \sqrt{\nu}\, K_1(\sqrt{\nu})}
#'       \exp\left(-\sqrt{\nu + \left(\dfrac{y-\mu}{\sigma}\right)^2}\right),}
#' with \eqn{K_1} the modified Bessel function of the second kind. The
#' exponent is the negative pseudo-Huber loss, quadratic in the residual near
#' the location and linear far from it. That exponent is the pseudo-Huber
#' loss, and this density is its exponential.
#'
#' The normalizing constant is formed on the log scale through the
#' **exponentially scaled** Bessel function,
#' \eqn{\log K_1(x) = \log\{e^x K_1(x)\} - x}. The Bessel terms are
#' degree-homogeneous, so the scaled form is exact and stays finite where
#' \eqn{K_1(\sqrt{\nu})} itself underflows, verified to \eqn{\nu = 2000}.
#'
#' @param distrib A `PseudoHuberDistrib` object, from
#'   [pseudohuber_distrib()].
#' @param y A numeric vector of observations. Every real value is in the
#'   support, so no value is rejected.
#' @param theta A named list with components `mu`, `sigma` and `nu`, each a
#'   numeric vector of length 1 or of the length of `y`. A component of length
#'   1 is recycled. `sigma` and `nu` must be strictly positive.
#' @param log Logical of length 1. When `TRUE` the log-density is returned.
#'   Defaults to `FALSE`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of densities, of length
#'   `max(length(y), length(mu), length(sigma), length(nu))`, one value per
#'   observation.
#'
#' @section Notation:
#' \eqn{\mu} is the location, \eqn{\sigma > 0} the scale, \eqn{\nu > 0} the
#' shape and \eqn{K_1} the modified Bessel function of the second kind of order
#' one, `besselK(x, 1)` in R.
#'
#' @seealso [distrib_cdf.PseudoHuberDistrib()] for the distribution function,
#'   [distrib_gradient.PseudoHuberDistrib()] for the derivatives of the
#'   log-density, [laplace_distrib()] and [gaussian1_distrib()] for the two
#'   limits, and [distrib_pdf()] for the generic.
#'
#' @examples
#' d <- pseudohuber_distrib()
#' y <- c(-2.5, 0.3, 1.8)
#' th <- list(mu = 0.4, sigma = 1.2, nu = 2)
#'
#' # The formula written out, with the Bessel constant.
#' D <- sqrt(2 + ((y - 0.4) / 1.2)^2)
#' all.equal(distrib_pdf(d, y, th),
#'           exp(-D) / (2 * 1.2 * sqrt(2) * besselK(sqrt(2), 1)))
#'
#' # It integrates to one.
#' integrate(function(v) distrib_pdf(d, v, th), -Inf, Inf)$value
#'
#' # A small shape is a Laplace of scale sigma; a large one a Gaussian of
#' # standard deviation sigma * nu^(1/4).
#' yy <- c(0.5, 1, 2, 4)
#' max(abs(distrib_pdf(d, yy, list(mu = 0, sigma = 1, nu = 1e-8)) -
#'         0.5 * exp(-abs(yy))))
#' max(abs(distrib_pdf(d, yy, list(mu = 0, sigma = 1, nu = 1e8)) -
#'         dnorm(yy, 0, 100)))
S7::method(distrib_pdf, PseudoHuberDistrib) <- function(distrib, y, theta, log = FALSE, ...) {
  mu <- theta[[1]]
  sigma <- theta[[2]]
  nu <- theta[[3]]

  sq_nu <- sqrt(nu)
  D <- sqrt(nu + ((y - mu) / sigma)^2)

  # log K_1(x) = log(besselK(x, 1, scaled)) - x  (scaled version avoids underflow)
  log_norm <- log(2) + log(sigma) + 0.5 * log(nu) + log(besselK(sq_nu, 1, expon.scaled = TRUE)) - sq_nu
  log_val <- -D - log_norm

  if (log) log_val else exp(log_val)
}

#' @title Pseudo-Huber Cumulative Distribution Function
#' @name distrib_cdf.PseudoHuberDistrib
#' @description
#' Computes \eqn{F(q) = P(Y \le q)} by numerical integration of the density.
#' The family has no elementary distribution function, so there is nothing
#' closed form to call.
#'
#' Two devices keep the quadrature honest. The law is symmetric about
#' \eqn{\mu}, so a quantile above the location is **reflected**,
#' \eqn{F(q) = 1 - F(2\mu - q)}, and only the lower tail is ever integrated,
#' where the integrand decays away from a finite endpoint. And every quantile
#' is one **row** of a single batched quadrature through [quad_rows()], so a
#' vector of `q` is integrated in a single call.
#'
#' A row that fails to reach the requested accuracy signals an error naming the
#' positions, instead of returning a plausible number.
#'
#' @param distrib A `PseudoHuberDistrib` object, from
#'   [pseudohuber_distrib()].
#' @param q A numeric vector of quantiles.
#' @param theta A named list with components `mu`, `sigma` and `nu`, each a
#'   numeric vector of length 1 or of the length of `q`. A component of length
#'   1 is recycled; a vector gives one integration per parameter setting.
#'   `sigma` and `nu` must be strictly positive.
#' @param lower.tail Logical of length 1. When `TRUE`, the default,
#'   probabilities are \eqn{P(Y \le q)}; when `FALSE` they are \eqn{P(Y > q)},
#'   formed as \eqn{1 - F}.
#' @param log.p Logical of length 1. When `TRUE` the logarithm of the
#'   probability is returned, taken after the quadrature, so it carries the
#'   quadrature's own accuracy rather than improving on it in the far tail.
#'   Defaults to `FALSE`.
#'
#' @return A numeric vector of probabilities in \eqn{[0, 1]}, of length
#'   `max(length(q), length(mu), length(sigma), length(nu))`, clamped to that
#'   range.
#'
#' @seealso [distrib_quantile.PseudoHuberDistrib()] for the inverse,
#'   [distrib_pdf.PseudoHuberDistrib()] for the integrand,
#'   [distrib_grad_cdf.PseudoHuberDistrib()] for the derivatives of this
#'   function, which are closed form in the location and the scale, and
#'   [distrib_cdf()] for the generic.
#'
#' @examples
#' d <- pseudohuber_distrib()
#' th <- list(mu = 0.4, sigma = 1.2, nu = 2)
#'
#' # The quadrature, and the symmetry the method exploits.
#' distrib_cdf(d, c(-1, 0.4, 2), th)
#' distrib_cdf(d, 0.4 - 1.5, th) + distrib_cdf(d, 0.4 + 1.5, th)
#'
#' # At the location it is one half exactly, the law being symmetric.
#' distrib_cdf(d, 0.4, th)
#'
#' # It agrees with a direct integration of the density.
#' c(method = distrib_cdf(d, 2, th),
#'   integral = integrate(function(v) distrib_pdf(d, v, th), -Inf, 2)$value)
S7::method(distrib_cdf, PseudoHuberDistrib) <- function(distrib, q, theta, lower.tail = TRUE, log.p = FALSE) {
  all_params <- expand_params(c(list(.q = q), theta))
  qv <- all_params$.q
  th_cols <- all_params[distrib@params]

  # By symmetry F(q) = 1 - F(2 mu - q): integrate only over the lower tail,
  # where the integrand is anchored at the finite endpoint. Every quantile is
  # one row of a single batched quadrature.
  mu <- th_cols[[1L]]
  left <- qv <= mu
  up <- ifelse(left, qv, 2 * mu - qv)

  integrand <- function(x, i) {
    xv <- as.numeric(x)
    idx <- rep(i, times = ncol(x))
    distrib_pdf(distrib, xv, lapply(th_cols, function(v) v[idx]))
  }
  vals <- quad_rows(integrand, -Inf, up)
  if (anyNA(vals)) {
    stop(sprintf(
      "The cdf quadrature did not reach the requested accuracy at quantile(s) %s.",
      paste(which(is.na(vals)), collapse = ", ")
    ), call. = FALSE)
  }
  res <- ifelse(left, vals, 1 - vals)

  res <- pmin(pmax(res, 0), 1)
  if (!lower.tail) res <- 1 - res
  if (log.p) log(res) else res
}

#' @title Pseudo-Huber Quantile Function
#' @name distrib_quantile.PseudoHuberDistrib
#' @description
#' Computes \eqn{Q(p)} by root-finding on the numerical distribution function,
#' the family having no elementary quantile. Symmetry about \eqn{\mu} does most
#' of the work: \eqn{Q(1/2) = \mu} exactly and without a search, and a
#' probability above one half is reflected through
#' \eqn{Q(p) = 2\mu - Q(1-p)}, so the search always runs in the lower half.
#'
#' The bracket starts at ten standard deviations below the location and doubles
#' its width until it contains the root, at most 100 times.
#' [stats::uniroot()] then closes it to `sqrt(.Machine$double.eps)`. Each
#' evaluation of the objective is a quadrature, so this is by far the dearest
#' method of the family.
#'
#' @param distrib A `PseudoHuberDistrib` object, from
#'   [pseudohuber_distrib()].
#' @param p A numeric vector of probabilities in \eqn{[0, 1]}, or of their
#'   logarithms when `log.p = TRUE`. `NaN` is returned for `NA` and for a value
#'   outside the range; 0 gives `-Inf` and 1 gives `Inf`.
#' @param theta A named list with components `mu`, `sigma` and `nu`, each a
#'   numeric vector of length 1 or of the length of `p`. A component of length
#'   1 is recycled. `sigma` and `nu` must be strictly positive.
#' @param lower.tail Logical of length 1. When `TRUE`, the default, `p` is
#'   \eqn{P(Y \le q)}; when `FALSE` it is \eqn{P(Y > q)}.
#' @param log.p Logical of length 1. When `TRUE`, `p` is read as a logarithm.
#'   Defaults to `FALSE`.
#'
#' @return A numeric vector of quantiles on \eqn{[-\infty, \infty]}, of length
#'   `max(length(p), length(mu), length(sigma), length(nu))`.
#'
#' @seealso [distrib_cdf.PseudoHuberDistrib()] for the function inverted here,
#'   [distrib_rng.PseudoHuberDistrib()], which draws by inverting it at uniform
#'   variates, and [distrib_quantile()] for the generic.
#'
#' @examples
#' d <- pseudohuber_distrib()
#' th <- list(mu = 0.4, sigma = 1.2, nu = 2)
#'
#' # The median is the location, returned without a search.
#' distrib_quantile(d, 0.5, th)
#'
#' # The quartiles, and the round trip back through the distribution function.
#' q <- distrib_quantile(d, c(0.25, 0.75), th)
#' q
#' distrib_cdf(d, q, th)
#'
#' # Symmetry: the two quartiles are equidistant from the location.
#' q - 0.4
S7::method(distrib_quantile, PseudoHuberDistrib) <- function(distrib, p, theta, lower.tail = TRUE, log.p = FALSE) {
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

    # Reflect to the lower half where the CDF is computed directly
    reflect <- pi > 0.5
    p_low <- if (reflect) 1 - pi else pi

    # Bracket [lo, m] expanding by multiples of the standard deviation
    sd_th <- sqrt(variance(distrib, th))
    lo <- m - 10 * sd_th
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

#' @title Pseudo-Huber Random Number Generator
#' @name distrib_rng.PseudoHuberDistrib
#' @description
#' Draws `n` independent variates by inverse transform: uniform variates from
#' [stats::runif()] passed through
#' [distrib_quantile.PseudoHuberDistrib()]. Each draw therefore costs a
#' root-find over a quadrature, which makes this the slowest generator in the
#' package; a sample of a few thousand is comfortable, a sample of a million is
#' not. The draws depend on `.Random.seed` in the usual way.
#'
#' @param distrib A `PseudoHuberDistrib` object, from
#'   [pseudohuber_distrib()].
#' @param n A single positive integer, the number of draws.
#' @param theta A named list with components `mu`, `sigma` and `nu`, each a
#'   numeric vector of length 1 or of length `n`. A component of length 1 is
#'   recycled, so a vector of length `n` draws one variate per parameter
#'   setting. `sigma` and `nu` must be strictly positive.
#'
#' @return A numeric vector of `n` draws.
#'
#' @seealso [distrib_quantile.PseudoHuberDistrib()] for the inversion,
#'   [fit_distrib()] to estimate the parameters back from a sample, and
#'   [distrib_rng()] for the generic.
#'
#' @examples
#' d <- pseudohuber_distrib()
#' th <- list(mu = 0.4, sigma = 1.2, nu = 2)
#'
#' # The sample moments recover the location and the Bessel-ratio variance.
#' set.seed(6)
#' z <- distrib_rng(d, 2000, th)
#' rbind(sample = c(mean(z), var(z)),
#'       theoretical = c(mean(d, th), variance(d, th)))
S7::method(distrib_rng, PseudoHuberDistrib) <- function(distrib, n, theta) {
  distrib_quantile(distrib, stats::runif(n), theta)
}

#' @title Pseudo-Huber Score
#' @name distrib_gradient.PseudoHuberDistrib
#' @description
#' Computes the first derivatives of the pseudo-Huber log-density with respect
#' to \eqn{\mu}, \eqn{\sigma} and \eqn{\nu}, one value per observation, in
#' closed form. With \eqn{r = y - \mu} and
#' \eqn{D = \sqrt{\nu + (r/\sigma)^2}},
#' \deqn{\dfrac{\partial \ell}{\partial \mu} = \dfrac{r}{\sigma^2 D}, \qquad
#'       \dfrac{\partial \ell}{\partial \sigma}
#'         = \dfrac{1}{\sigma}\left(\dfrac{r^2}{\sigma^2 D} - 1\right),}
#' \deqn{\dfrac{\partial \ell}{\partial \nu} = -\dfrac{1}{2}\left[
#'   \dfrac{1}{\nu} + \dfrac{1}{D}
#'   + \dfrac{K_1'(\sqrt{\nu})}{\sqrt{\nu}\, K_1(\sqrt{\nu})}\right].}
#'
#' The location component **redescends**: it grows like \eqn{r/\sigma^2} near
#' the location and tends to \eqn{\pm 1/\sigma} far from it, so a gross outlier
#' contributes a bounded amount to the estimating equation. That bounded
#' influence is the whole reason the pseudo-Huber loss exists, and here it is
#' the score of a genuine likelihood.
#'
#' With `scale = "link"` the generic applies the chain rule for the links the
#' family carries. This method always returns the parameter scale. The
#' arithmetic runs in a compiled kernel; the Bessel ratio enters only the
#' \eqn{\nu} component and is formed from the exponentially scaled functions.
#'
#' @param distrib A `PseudoHuberDistrib` object, from
#'   [pseudohuber_distrib()].
#' @param y A numeric vector of observations.
#' @param theta A named list with components `mu`, `sigma` and `nu`, each a
#'   numeric vector of length 1 or of the length of `y`. A component of length
#'   1 is recycled. `sigma` and `nu` must be strictly positive.
#' @param scale One of `"parameter"` (the default) or `"link"`, matched by
#'   [base::match.arg()]. Read by the generic, not by this method.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A named list of three numeric vectors, `mu`, `sigma` and `nu`, each
#'   of length `max(length(y), length(mu), length(sigma), length(nu))`.
#'
#' @section Notation:
#' \eqn{\ell} is the log-density of one observation, \eqn{\mu} the location,
#' \eqn{\sigma > 0} the scale, \eqn{\nu > 0} the shape, \eqn{r = y - \mu} and
#' \eqn{K_1} the modified Bessel function of the second kind of order one.
#'
#' @seealso [distrib_hessian.PseudoHuberDistrib()] for the second derivatives,
#'   [distrib_expected_hessian.PseudoHuberDistrib()] for their expectation,
#'   [distrib_grad_y.PseudoHuberDistrib()] for the derivative in the response,
#'   and [distrib_gradient()] for the generic.
#'
#' @examples
#' d <- pseudohuber_distrib()
#' y <- c(-2.5, 0.3, 1.8)
#' th <- list(mu = 0.4, sigma = 1.2, nu = 2)
#' g <- distrib_gradient(d, y, th)
#'
#' # The location and scale components, written out.
#' r <- y - 0.4; D <- sqrt(2 + (r / 1.2)^2)
#' all.equal(g$mu, r / (1.2^2 * D))
#' all.equal(g$sigma, (r^2 / (1.2^2 * D) - 1) / 1.2)
#'
#' # numDeriv on the summed log-density reproduces the summed score.
#' fn <- function(p)
#'   sum(distrib_pdf(d, y, list(mu = p[1], sigma = p[2], nu = p[3]), log = TRUE))
#' rbind(numeric = numDeriv::grad(fn, c(0.4, 1.2, 2)),
#'       closed = vapply(g, sum, numeric(1)))
#'
#' # The location score is bounded by 1 / sigma, where a Gaussian's grows
#' # without bound.
#' rr <- c(1, 4, 16, 64)
#' rbind(residual = rr,
#'       pseudohuber = rr / (1.2^2 * sqrt(2 + (rr / 1.2)^2)),
#'       gaussian = rr / 1.2^2,
#'       bound = rep(1 / 1.2, 4))
S7::method(distrib_gradient, PseudoHuberDistrib) <- function(distrib, y, theta, scale = c("parameter", "link"), ...) {
  pseudohuber_gradient_cpp(y, theta[[1]], theta[[2]], theta[[3]])
}

#' @title Pseudo-Huber Observed Hessian
#' @name distrib_hessian.PseudoHuberDistrib
#' @description
#' Computes the six distinct second derivatives of the pseudo-Huber log-density
#' with respect to \eqn{\mu}, \eqn{\sigma} and \eqn{\nu}, one value per
#' observation, in closed form. With \eqn{r = y - \mu},
#' \eqn{D = \sqrt{\nu + (r/\sigma)^2}},
#' \eqn{R_1 = K_1'(\sqrt{\nu})/K_1(\sqrt{\nu})} and
#' \eqn{R_2 = K_1''(\sqrt{\nu})/K_1(\sqrt{\nu})},
#' \deqn{\dfrac{\partial^2 \ell}{\partial \mu^2} = -\dfrac{\nu}{\sigma^2 D^3},
#'       \qquad
#'       \dfrac{\partial^2 \ell}{\partial \sigma^2}
#'         = \dfrac{\sigma^4 - 3\sigma^2 r^2 D^{-1} + r^4 D^{-3}}{\sigma^6},}
#' \deqn{\dfrac{\partial^2 \ell}{\partial \nu^2} = \dfrac{1}{4D^3}
#'   + \dfrac{1}{2\nu^2}
#'   + \dfrac{1}{4}\left(\dfrac{R_1}{\nu^{3/2}} + \dfrac{R_1^2}{\nu}
#'     - \dfrac{R_2}{\nu}\right),}
#' \deqn{\dfrac{\partial^2 \ell}{\partial \mu \, \partial \sigma}
#'         = \dfrac{-2\nu\sigma^2 r - r^3}{\sigma^2(\nu\sigma^2 + r^2)^{3/2}},
#'       \qquad
#'       \dfrac{\partial^2 \ell}{\partial \mu \, \partial \nu}
#'         = -\dfrac{r}{2\sigma^2 D^3}, \qquad
#'       \dfrac{\partial^2 \ell}{\partial \sigma \, \partial \nu}
#'         = -\dfrac{r^2}{2\sigma^3 D^3}.}
#'
#' The curvature in \eqn{\mu} is negative at every observation, unlike a
#' Student t's, so the log-density is concave in the location however far out
#' the residual is. What redescends is the score, not the curvature's sign.
#'
#' @param distrib A `PseudoHuberDistrib` object, from
#'   [pseudohuber_distrib()].
#' @param y A numeric vector of observations.
#' @param theta A named list with components `mu`, `sigma` and `nu`, each a
#'   numeric vector of length 1 or of the length of `y`. A component of length
#'   1 is recycled. `sigma` and `nu` must be strictly positive.
#' @param scale One of `"parameter"` (the default) or `"link"`, matched by
#'   [base::match.arg()]. Read by the generic, not by this method.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A named list of six numeric vectors, `mu_mu`, `sigma_sigma`,
#'   `nu_nu`, `mu_sigma`, `mu_nu` and `sigma_nu`, each of length
#'   `max(length(y), length(mu), length(sigma), length(nu))`. The six name the
#'   distinct entries of a symmetric \eqn{3 \times 3} matrix per observation.
#'
#' @section Notation:
#' \eqn{\ell} is the log-density of one observation, \eqn{\mu} the location,
#' \eqn{\sigma > 0} the scale, \eqn{\nu > 0} the shape, \eqn{r = y - \mu} and
#' \eqn{K_1} the modified Bessel function of the second kind of order one.
#'
#' @seealso [distrib_gradient.PseudoHuberDistrib()] for the score,
#'   [distrib_expected_hessian.PseudoHuberDistrib()] for the expectation of
#'   this quantity, [distrib_deriv3.PseudoHuberDistrib()] for the order above,
#'   and [distrib_hessian()] for the generic.
#'
#' @examples
#' d <- pseudohuber_distrib()
#' y <- c(-2.5, 0.3, 1.8)
#' th <- list(mu = 0.4, sigma = 1.2, nu = 2)
#' h <- distrib_hessian(d, y, th)
#' names(h)
#'
#' # The location entry, written out, and its sign.
#' r <- y - 0.4; D <- sqrt(2 + (r / 1.2)^2)
#' all.equal(h$mu_mu, -2 / (1.2^2 * D^3))
#'
#' # numDeriv on the summed log-density reproduces the summed matrix.
#' fn <- function(p)
#'   sum(distrib_pdf(d, y, list(mu = p[1], sigma = p[2], nu = p[3]), log = TRUE))
#' H <- numDeriv::hessian(fn, c(0.4, 1.2, 2))
#' rbind(numeric = c(H[1, 1], H[2, 2], H[3, 3], H[1, 2], H[1, 3], H[2, 3]),
#'       closed = vapply(h, sum, numeric(1)))
#'
#' # The curvature in the location stays negative however far out y is.
#' distrib_hessian(d, 0.4 + c(1, 10, 100), th)$mu_mu
S7::method(distrib_hessian, PseudoHuberDistrib) <- function(distrib, y, theta, scale = c("parameter", "link"), ...) {
  pseudohuber_hessian_cpp(y, theta[[1]], theta[[2]], theta[[3]])
}

#' @title Pseudo-Huber Expected Hessian
#' @name distrib_expected_hessian.PseudoHuberDistrib
#' @description
#' Returns the expectation of the observed Hessian under the model. **There is
#' no closed form**, so the four components that do not vanish are obtained by
#' the strategy `approx` names, normally a numerical integration of the
#' observed Hessian against the density through [expectation()]. The two
#' components containing \eqn{\mu} an odd number of times are then **replaced
#' by exact zeros**: the law is symmetric about \eqn{\mu}, so
#' \eqn{\mathbb{E}[r] = \mathbb{E}[r^3] = 0} and the \eqn{\mu\sigma} and
#' \eqn{\mu\nu} entries vanish. The location is therefore orthogonal to both
#' other parameters, and \eqn{\hat\mu} is asymptotically independent of them.
#'
#' The method **improves** the approximation rather than replacing it, which is
#' why [expected_hessian_exact.PseudoHuberDistrib()] answers `FALSE`. Reading
#' the method's owning class would say the family writes its information out;
#' it does not, and the cost says so: measured at 100 observations this takes
#' about 11 seconds, where the families that do write it out answer in a median
#' of 0.183 milliseconds.
#'
#' @param distrib A `PseudoHuberDistrib` object, from
#'   [pseudohuber_distrib()].
#' @param y A numeric vector of observations. Its length sets the length of
#'   each returned component.
#' @param theta A named list with components `mu`, `sigma` and `nu`, each a
#'   numeric vector of length 1 or of the length of `y`. `sigma` and `nu` must
#'   be strictly positive.
#' @param scale One of `"parameter"` (the default) or `"link"`, matched by
#'   [base::match.arg()]. Read by the generic, not by this method.
#' @param approx One of `"bartlett"` (the default), `"integrate"`, `"mc"` or
#'   `"opg"`, the strategy [expected_derivative()] uses. **Read here**, unlike
#'   on the families that write their information out.
#' @param nsim A single positive integer, the sample size when
#'   `approx = "mc"`. Defaults to `10000`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A named list of six numeric vectors, `mu_mu`, `sigma_sigma`,
#'   `nu_nu`, `mu_sigma`, `mu_nu` and `sigma_nu`, each of length `length(y)`.
#'   `mu_sigma` and `mu_nu` are exactly zero.
#'
#' @seealso [distrib_hessian.PseudoHuberDistrib()] for the quantity this is the
#'   expectation of, [expected_hessian_exact.PseudoHuberDistrib()] for the
#'   predicate that reports this is not a closed form,
#'   [fisher_scoring()], which reads that predicate, and
#'   [distrib_expected_hessian()] for the generic.
#'
#' @examples
#' d <- pseudohuber_distrib()
#' y <- c(-2.5, 0.3, 1.8)
#' th <- list(mu = 0.4, sigma = 1.2, nu = 2)
#' eh <- distrib_expected_hessian(d, y, th)
#' vapply(eh, function(v) v[1], numeric(1))
#'
#' # The two entries odd in the residual are exactly zero by symmetry, so the
#' # location is orthogonal to the scale and the shape.
#' c(eh$mu_sigma[1], eh$mu_nu[1])
#'
#' # Unlike the families that write their information out, this one reads
#' # `approx`: a Monte Carlo strategy gives a different, noisier answer.
#' set.seed(1)
#' vapply(distrib_expected_hessian(d, y, th, approx = "mc", nsim = 2000),
#'        function(v) v[1], numeric(1))
S7::method(distrib_expected_hessian, PseudoHuberDistrib) <- function(distrib, y, theta, scale = c("parameter", "link"), approx = c("bartlett", "integrate", "mc", "opg"), nsim = 10000, ...) {
  n <- length(y)
  out <- expected_derivative(distrib, y, theta, order = 2L,
                             approx = match.arg(approx), nsim = nsim)
  # exact by symmetry: E[r] = E[r^3] = 0
  out$mu_sigma <- rep(0, n)
  out$mu_nu <- rep(0, n)
  out
}

#' @title The Pseudo-Huber Does Not Write Its Expected Information Out
#' @name expected_hessian_exact.PseudoHuberDistrib
#' @description
#' Answers `FALSE`, declaring that this family's expected information is a
#' numerical approximation and not a formula, so that callers who branch on the
#' distinction branch correctly.
#'
#' @details
#' The predicate's default reads the class a method is registered on, which
#' here would answer `TRUE` and be wrong.
#' [distrib_expected_hessian.PseudoHuberDistrib()] is registered on this class,
#' but what it does is call [expected_derivative()] and then replace the two
#' components that vanish by symmetry: it **improves** the approximation rather
#' than replacing it.
#'
#' The cost is the discriminator. Measured at 100 observations the method takes
#' about 11 seconds, where the families that do write their information out
#' answer in a median of 0.183 milliseconds. Two consequences were live before
#' the declaration: [fit_distrib()] rejected a legitimate
#' `fisher_scoring(approx = )` here with a message saying the family computes
#' its expected information in closed form, which is untrue; and its
#' standard-error branch entered a multi-second quadrature believing it a
#' formula.
#'
#' @param x A `PseudoHuberDistrib` object, from [pseudohuber_distrib()].
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return `FALSE`, a logical of length 1.
#'
#' @seealso [expected_hessian_exact()] for the generic and its default,
#'   [distrib_expected_hessian.PseudoHuberDistrib()] for the method this
#'   describes, and [fisher_scoring()] for the consumer.
#' @keywords internal
S7::method(expected_hessian_exact, PseudoHuberDistrib) <- function(x, ...) {
  FALSE
}

#' @title Pseudo-Huber Third-Order Derivatives
#' @name distrib_deriv3.PseudoHuberDistrib
#' @description
#' Computes the ten distinct third derivatives of the pseudo-Huber log-density
#' in \eqn{\mu}, \eqn{\sigma} and \eqn{\nu}, **in closed form**, in a compiled
#' kernel. Every component but the pure-\eqn{\nu} one is a rational function of
#' \eqn{r = y - \mu}, \eqn{\sigma} and \eqn{D = \sqrt{\nu + (r/\sigma)^2}};
#' Bessel functions enter through \eqn{\nu} alone, and the exponentially scaled
#' forms are used so that a large \eqn{\nu} does not overflow.
#'
#' **The expected third derivatives have no closed form.** With
#' `expected = TRUE` the method calls [expected_derivative()], which integrates
#' the observed derivatives against the density by the strategy `approx` names.
#' That is the one place on this page where `approx` and `nsim` are read.
#'
#' @param distrib A `PseudoHuberDistrib` object, from
#'   [pseudohuber_distrib()].
#' @param y A numeric vector of observations. With `expected = TRUE` only its
#'   length is read.
#' @param theta A named list with components `mu`, `sigma` and `nu`, each a
#'   numeric vector of length 1 or of the length of `y`. A component of length
#'   1 is recycled. `sigma` and `nu` must be strictly positive.
#' @param expected Logical of length 1. When `TRUE` the expectation under the
#'   model is returned in place of the value at the data, computed numerically.
#'   Defaults to `FALSE`.
#' @param scale One of `"parameter"` (the default) or `"link"`, matched by
#'   [base::match.arg()]. Read by the generic, not by this method.
#' @param approx One of `"integrate"` (the default here), `"bartlett"`, `"mc"`
#'   or `"opg"`. Read only when `expected = TRUE`.
#' @param nsim A single positive integer, the sample size when
#'   `approx = "mc"`. Read only when `expected = TRUE`. Defaults to `10000`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A named list of ten numeric vectors, `mu_mu_mu` through
#'   `nu_nu_nu`, each of length
#'   `max(length(y), length(mu), length(sigma), length(nu))`.
#'
#' @seealso [distrib_hessian.PseudoHuberDistrib()] for the order below,
#'   [distrib_deriv4.PseudoHuberDistrib()] for the order above,
#'   [expected_derivative()] for the numerical expectation, and
#'   [distrib_deriv3()] for the generic.
#'
#' @examples
#' d <- pseudohuber_distrib()
#' y <- c(-2.5, 0.3, 1.8)
#' th <- list(mu = 0.4, sigma = 1.2, nu = 2)
#' d3 <- distrib_deriv3(d, y, th)
#' names(d3)
#'
#' # A central difference of the Hessian reproduces the pure-location
#' # component, which is what says the closed form is the right one.
#' eps <- 1e-5
#' up <- distrib_hessian(d, y, list(mu = 0.4 + eps, sigma = 1.2, nu = 2))$mu_mu
#' dn <- distrib_hessian(d, y, list(mu = 0.4 - eps, sigma = 1.2, nu = 2))$mu_mu
#' all.equal((up - dn) / (2 * eps), d3$mu_mu_mu, tolerance = 1e-6)
#'
#' # The expected branch is a quadrature and takes a strategy; the components
#' # odd in the residual come back at zero.
#' vapply(distrib_deriv3(d, 0.4, th, expected = TRUE),
#'        function(v) v[1], numeric(1))
S7::method(distrib_deriv3, PseudoHuberDistrib) <- function(distrib, y, theta, expected = FALSE, scale = c("parameter", "link"), approx = c("integrate", "bartlett", "mc", "opg"), nsim = 10000, ...) {
  if (expected) {
    expected_derivative(distrib, y, theta, order = 3L,
                        approx = match.arg(approx), nsim = nsim)
  } else {
    pseudohuber_deriv3_cpp(y, theta[[1]], theta[[2]], theta[[3]])
  }
}

#' @title Pseudo-Huber Fourth-Order Derivatives
#' @name distrib_deriv4.PseudoHuberDistrib
#' @description
#' Computes the fifteen distinct fourth derivatives of the pseudo-Huber
#' log-density in \eqn{\mu}, \eqn{\sigma} and \eqn{\nu}, **in closed form**, in
#' a compiled kernel. The Bessel handling is
#' [distrib_deriv3.PseudoHuberDistrib()]'s: the functions enter through
#' \eqn{\nu} alone, in their exponentially scaled forms.
#'
#' **The expected fourth derivatives have no closed form.** With
#' `expected = TRUE` the method calls [expected_derivative()], which is the one
#' place on this page where `approx` and `nsim` are read.
#'
#' @param distrib A `PseudoHuberDistrib` object, from
#'   [pseudohuber_distrib()].
#' @param y A numeric vector of observations. With `expected = TRUE` only its
#'   length is read.
#' @param theta A named list with components `mu`, `sigma` and `nu`, each a
#'   numeric vector of length 1 or of the length of `y`. A component of length
#'   1 is recycled. `sigma` and `nu` must be strictly positive.
#' @param expected Logical of length 1. When `TRUE` the expectation under the
#'   model is returned in place of the value at the data, computed numerically.
#'   Defaults to `FALSE`.
#' @param scale One of `"parameter"` (the default) or `"link"`, matched by
#'   [base::match.arg()]. Read by the generic, not by this method.
#' @param approx One of `"integrate"` (the default here), `"bartlett"`, `"mc"`
#'   or `"opg"`. Read only when `expected = TRUE`.
#' @param nsim A single positive integer, the sample size when
#'   `approx = "mc"`. Read only when `expected = TRUE`. Defaults to `10000`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A named list of fifteen numeric vectors named for the multi-index
#'   they carry, from `mu_mu_mu_mu` to `nu_nu_nu_nu`, each of length
#'   `max(length(y), length(mu), length(sigma), length(nu))`.
#'
#' @seealso [distrib_deriv3.PseudoHuberDistrib()] for the order below and the
#'   Bessel handling, [distrib_hessian.PseudoHuberDistrib()] for the second
#'   order, [expected_derivative()] for the numerical expectation, and
#'   [distrib_deriv4()] for the generic.
#'
#' @examples
#' d <- pseudohuber_distrib()
#' y <- c(-2.5, 0.3, 1.8)
#' th <- list(mu = 0.4, sigma = 1.2, nu = 2)
#' d4 <- distrib_deriv4(d, y, th)
#' length(d4)
#' names(d4)[1:4]
#'
#' # A central difference of the third order reproduces the pure-location
#' # component.
#' eps <- 1e-4
#' up <- distrib_deriv3(d, y, list(mu = 0.4 + eps, sigma = 1.2, nu = 2))$mu_mu_mu
#' dn <- distrib_deriv3(d, y, list(mu = 0.4 - eps, sigma = 1.2, nu = 2))$mu_mu_mu
#' all.equal((up - dn) / (2 * eps), d4$mu_mu_mu_mu, tolerance = 1e-5)
S7::method(distrib_deriv4, PseudoHuberDistrib) <- function(distrib, y, theta, expected = FALSE, scale = c("parameter", "link"), approx = c("integrate", "bartlett", "mc", "opg"), nsim = 10000, ...) {
  if (expected) {
    expected_derivative(distrib, y, theta, order = 4L,
                        approx = match.arg(approx), nsim = nsim)
  } else {
    pseudohuber_deriv4_cpp(y, theta[[1]], theta[[2]], theta[[3]])
  }
}

#' @title Pseudo-Huber First Derivative in the Response
#' @name distrib_grad_y.PseudoHuberDistrib
#' @description
#' Computes \eqn{\partial \ell / \partial y} in closed form. With
#' \eqn{r = y - \mu} and \eqn{D = \sqrt{\nu + (r/\sigma)^2}},
#' \deqn{\dfrac{\partial \ell}{\partial y} = -\dfrac{r}{\sigma^2 D}.}
#' The family is a location family in \eqn{\mu}, so this is exactly the
#' negative of the location score
#' [distrib_gradient.PseudoHuberDistrib()]`$mu`, and like it, it is
#' **bounded**: it tends to \eqn{\mp 1/\sigma} as the residual runs away. This
#' quantity is what a quantile residual's delta-method standard error and a
#' change of variable in the response both need.
#'
#' @param distrib A `PseudoHuberDistrib` object, from
#'   [pseudohuber_distrib()].
#' @param y A numeric vector of observations.
#' @param theta A named list with components `mu`, `sigma` and `nu`, each a
#'   numeric vector of length 1 or of the length of `y`. A component of length
#'   1 is recycled. `sigma` and `nu` must be strictly positive.
#'
#' @return A numeric vector of length
#'   `max(length(y), length(mu), length(sigma), length(nu))`, one value per
#'   observation.
#'
#' @seealso [distrib_hess_y.PseudoHuberDistrib()] for the second derivative in
#'   the response, [distrib_cross_y.PseudoHuberDistrib()] for the mixed
#'   derivative, [distrib_gradient.PseudoHuberDistrib()] for the derivatives in
#'   the parameters, and [distrib_grad_y()] for the generic.
#'
#' @examples
#' d <- pseudohuber_distrib()
#' y <- c(-2.5, 0.3, 1.8)
#' th <- list(mu = 0.4, sigma = 1.2, nu = 2)
#'
#' # The closed form, written out.
#' r <- y - 0.4; D <- sqrt(2 + (r / 1.2)^2)
#' all.equal(distrib_grad_y(d, y, th), -r / (1.2^2 * D))
#'
#' # A location family, so this is minus the location score.
#' all.equal(distrib_grad_y(d, y, th), -distrib_gradient(d, y, th)$mu)
#'
#' # A central difference of the log-density in y reproduces it.
#' eps <- 1e-6
#' all.equal((distrib_pdf(d, y + eps, th, log = TRUE) -
#'            distrib_pdf(d, y - eps, th, log = TRUE)) / (2 * eps),
#'           distrib_grad_y(d, y, th), tolerance = 1e-6)
#'
#' # Bounded by 1 / sigma however far out the response is.
#' c(distrib_grad_y(d, 0.4 + c(10, 1000, 1e6), th), bound = -1 / 1.2)
S7::method(distrib_grad_y, PseudoHuberDistrib) <- function(distrib, y, theta) {
  mu <- theta[[1]]; sigma <- theta[[2]]; nu <- theta[[3]]
  r <- y - mu
  sigma2 <- sigma^2
  D <- sqrt(nu + r^2 / sigma2)
  -r / (sigma2 * D)
}

#' @title Pseudo-Huber Second Derivative in the Response
#' @name distrib_hess_y.PseudoHuberDistrib
#' @description
#' Computes \eqn{\partial^2 \ell / \partial y^2} in closed form. With
#' \eqn{r = y - \mu} and \eqn{D = \sqrt{\nu + (r/\sigma)^2}},
#' \deqn{\dfrac{\partial^2 \ell}{\partial y^2} = -\dfrac{\nu}{\sigma^2 D^3}.}
#' The family is a location family, so this equals the pure-location entry of
#' [distrib_hessian.PseudoHuberDistrib()], with no sign change; two derivatives
#' in \eqn{\mu} carry two factors of \eqn{-1}. It is **negative everywhere**,
#' so the log-density is concave in the response, and it decays like
#' \eqn{|r|^{-3}} in the tails.
#'
#' @param distrib A `PseudoHuberDistrib` object, from
#'   [pseudohuber_distrib()].
#' @param y A numeric vector of observations.
#' @param theta A named list with components `mu`, `sigma` and `nu`, each a
#'   numeric vector of length 1 or of the length of `y`. A component of length
#'   1 is recycled. `sigma` and `nu` must be strictly positive.
#'
#' @return A numeric vector of length
#'   `max(length(y), length(mu), length(sigma), length(nu))`, one value per
#'   observation, every entry negative.
#'
#' @seealso [distrib_grad_y.PseudoHuberDistrib()] for the first derivative in
#'   the response, [distrib_hessian.PseudoHuberDistrib()] for the second
#'   derivatives in the parameters, and [distrib_hess_y()] for the generic.
#'
#' @examples
#' d <- pseudohuber_distrib()
#' y <- c(-2.5, 0.3, 1.8)
#' th <- list(mu = 0.4, sigma = 1.2, nu = 2)
#'
#' # The closed form, written out.
#' r <- y - 0.4; D <- sqrt(2 + (r / 1.2)^2)
#' all.equal(distrib_hess_y(d, y, th), -2 / (1.2^2 * D^3))
#'
#' # A location family, so this is the pure-location entry of the Hessian.
#' all.equal(distrib_hess_y(d, y, th), distrib_hessian(d, y, th)$mu_mu)
#'
#' # A central difference of the first derivative reproduces it.
#' eps <- 1e-5
#' all.equal((distrib_grad_y(d, y + eps, th) -
#'            distrib_grad_y(d, y - eps, th)) / (2 * eps),
#'           distrib_hess_y(d, y, th), tolerance = 1e-6)
#'
#' # Negative everywhere, decaying like |r|^-3 in the tails.
#' distrib_hess_y(d, 0.4 + c(1, 10, 100), th)
S7::method(distrib_hess_y, PseudoHuberDistrib) <- function(distrib, y, theta) {
  mu <- theta[[1]]; sigma <- theta[[2]]; nu <- theta[[3]]
  r <- y - mu
  sigma2 <- sigma^2
  D <- sqrt(nu + r^2 / sigma2)
  -nu / (sigma2 * D^3)
}

# --- CONSTRUCTOR WRAPPER ---

#' Pseudo-Huber Distribution, Location, Scale and Shape
#'
#' @description
#' Builds the distribution object for the pseudo-Huber family, whose
#' log-density is the negative pseudo-Huber loss
#' \eqn{-\sqrt{\nu + \{(y-\mu)/\sigma\}^2}}. It is the likelihood counterpart
#' of that loss: the score is bounded, so a gross outlier contributes a limited
#' amount to the estimating equation, and the shape \eqn{\nu > 0} interpolates
#' between a Laplace at \eqn{\nu \to 0} and a Gaussian at \eqn{\nu \to \infty}.
#' The family is the symmetric hyperbolic distribution, a special case of the
#' generalized hyperbolic.
#'
#' The returned object carries closed-form derivatives of the log-density to
#' fourth order in the parameters and closed first and second derivatives in
#' the response. The distribution function, the quantile and the expected
#' information are numerical; that is the price of the family.
#'
#' @param link_mu A `link` object from `linkfunctions7` for the location
#'   \eqn{\mu}. Defaults to [linkfunctions7::identity_link()], the location
#'   ranging over the whole line already.
#' @param link_sigma A `link` object from `linkfunctions7` for the scale
#'   \eqn{\sigma}. Defaults to [linkfunctions7::log_link()], which maps
#'   \eqn{(0, \infty)} onto the line and so keeps every fitted value positive.
#' @param link_nu A `link` object from `linkfunctions7` for the shape
#'   \eqn{\nu}. Defaults to [linkfunctions7::log_link()], for the same reason.
#'
#' @details
#' # The parametrization
#'
#' The density on \eqn{y \in (-\infty, \infty)} is
#' \deqn{f(y; \mu, \sigma, \nu) = \dfrac{1}{2 \sigma \sqrt{\nu}\, K_1(\sqrt{\nu})}
#'       \exp\left(-\sqrt{\nu + \left(\dfrac{y-\mu}{\sigma}\right)^2}\right),}
#' with \eqn{K_1} the modified Bessel function of the second kind. The exponent
#' is quadratic in the residual near the location, where
#' \eqn{\sqrt{\nu + z^2} \approx \sqrt{\nu} + z^2/(2\sqrt{\nu})}, and linear
#' far from it, where it is \eqn{|z|}. That is the pseudo-Huber loss, and this
#' family is its exponential.
#'
#' The Bessel terms are degree-homogeneous, so the exponentially scaled
#' `besselK(x, nu, expon.scaled = TRUE)` is exact and avoids the overflow the
#' unscaled form meets; the constant stays finite to \eqn{\nu = 2000}.
#'
#' # The two limits, with their rates
#'
#' At \eqn{\nu \to 0} the exponent is \eqn{|y-\mu|/\sigma} and the family is
#' the **Laplace** with scale \eqn{\sigma}: measured at \eqn{\mu = 0},
#' \eqn{\sigma = 1}, the largest gap over \eqn{y \in \{0.5, 1, 2, 4\}} is
#' `4.9e-05` at \eqn{\nu = 10^{-4}} and `1.9e-12` at \eqn{10^{-12}}, so the
#' approach is \eqn{O(\nu)}.
#'
#' At \eqn{\nu \to \infty} the exponent is
#' \eqn{\sqrt{\nu} + (y-\mu)^2/(2\sqrt{\nu}\sigma^2)} and the family is the
#' **Gaussian** with standard deviation \eqn{\sigma\nu^{1/4}}: the largest gap
#' is `1.5e-04` at \eqn{\nu = 10^4} and `1.5e-07` at \eqn{10^8}. Note that the
#' spread grows with \eqn{\nu}, so \eqn{\sigma} alone is not the standard
#' deviation.
#'
#' # Moments
#'
#' Every moment is a ratio of Bessel functions at \eqn{\sqrt{\nu}}:
#' \deqn{\mathbb{E}(Y) = \mu, \qquad
#'       \operatorname{Var}(Y) = \sigma^2 \sqrt{\nu}\,
#'         \dfrac{K_2(\sqrt{\nu})}{K_1(\sqrt{\nu})},}
#' the skewness is 0 by symmetry, and the excess kurtosis is
#' \eqn{3K_3(\sqrt{\nu})K_1(\sqrt{\nu})/K_2(\sqrt{\nu})^2 - 3}. The kurtosis
#' falls as \eqn{\nu} grows, from 2.93 at \eqn{\nu = 0.01} through 1.86 at 1 to
#' 0.030 at \eqn{10^4}, tending to the Gaussian's 0.
#'
#' # Derivatives
#'
#' Write \eqn{r = y - \mu} and \eqn{D = \sqrt{\nu + (r/\sigma)^2}}. The score is
#' \deqn{\dfrac{\partial \ell}{\partial \mu} = \dfrac{r}{\sigma^2 D}, \qquad
#'       \dfrac{\partial \ell}{\partial \sigma}
#'         = \dfrac{1}{\sigma}\left(\dfrac{r^2}{\sigma^2 D} - 1\right),}
#' \deqn{\dfrac{\partial \ell}{\partial \nu} = -\dfrac{1}{2}\left[
#'   \dfrac{1}{\nu} + \dfrac{1}{D}
#'   + \dfrac{K_1'(\sqrt{\nu})}{\sqrt{\nu}\, K_1(\sqrt{\nu})}\right].}
#' The location component is **bounded by** \eqn{1/\sigma}, which is the
#' family's whole point; the curvature in \eqn{\mu}, unlike a Student t's, is
#' negative at every observation, so the log-density stays concave in the
#' location. Orders three and four are closed form too, in compiled kernels, as
#' are the derivatives in the response and the mixed derivative
#' [distrib_cross_y.PseudoHuberDistrib()], whose \eqn{\nu} component is
#' \eqn{r/(2\sigma^2 D^3)}.
#'
#' # What is numerical, and what it costs
#'
#' The distribution function has no elementary form and is a quadrature,
#' batched over the quantiles and reflected about \eqn{\mu} so that only the
#' lower tail is integrated. The quantile inverts it by root-finding, and the
#' generator inverts it at uniform variates, so a sample costs one root-find
#' per draw.
#'
#' The **expected information has no closed form either**. Its four non-zero
#' components come from [expected_derivative()], and the two containing
#' \eqn{\mu} an odd number of times are replaced by exact zeros, the law being
#' symmetric. That makes the location orthogonal to the scale and the shape.
#' [expected_hessian_exact.PseudoHuberDistrib()] declares the approximation, so
#' a caller who branches on the distinction branches correctly; at 100
#' observations it costs about 11 seconds against a median of 0.183
#' milliseconds for a family that writes its information out.
#'
#' # Estimation
#'
#' [fit_distrib()] maximizes the log-likelihood on the link scale. Given the
#' cost of the expected information, `method = optimizers7::newton()` on the
#' observed Hessian is the cheaper route here, and it is closed form.
#'
#' @section Notation:
#' \eqn{\ell} is the log-density of one observation, \eqn{\mu} the location,
#' \eqn{\sigma > 0} the scale, \eqn{\nu > 0} the shape, \eqn{r = y - \mu},
#' \eqn{D = \sqrt{\nu + (r/\sigma)^2}} and \eqn{K_m} the modified Bessel
#' function of the second kind of order \eqn{m}, `besselK(x, m)` in R.
#' \eqn{\eta} is a parameter on the unconstrained scale of its link, with
#' \eqn{\theta = g^{-1}(\eta)}.
#'
#' @return An S7 object of class `PseudoHuberDistrib`, inheriting from
#'   `continuous_distrib`, with `distrib_name` `"pseudo huber"`, `dimension`
#'   `"univariate"`, `bounds` `c(-Inf, Inf)`, `params`
#'   `c("mu", "sigma", "nu")`, `n_params` `3`, `params_bounds` the domains
#'   \eqn{(-\infty, \infty)}, \eqn{(0, \infty)} and \eqn{(0, \infty)}, and
#'   `link_params` the three links given here.
#'
#' @references
#' Barndorff-Nielsen, O. (1978). Hyperbolic distributions and distributions on
#' hyperbolae. *Scandinavian Journal of Statistics*, **5**(3), 151-157.
#'
#' Charbonnier, P., Blanc-Feraud, L., Aubert, G. and Barlaud, M. (1997).
#' Deterministic edge-preserving regularization in computed imaging.
#' *IEEE Transactions on Image Processing*, **6**(2), 298-311.
#'
#' @importFrom linkfunctions7 identity_link log_link
#' @importFrom stats runif uniroot
#'
#' @examples
#' d <- pseudohuber_distrib()
#' d
#'
#' # The density integrates to one.
#' th <- list(mu = 0.4, sigma = 1.2, nu = 2)
#' integrate(function(v) distrib_pdf(d, v, th), -Inf, Inf)$value
#'
#' # Every moment is a Bessel ratio; the skewness is zero by symmetry.
#' c(mean = mean(d, th), var = variance(d, th),
#'   skew = skewness(d, th), kurt = kurtosis(d, th))
#'
#' # A small shape is a Laplace of scale sigma; a large one a Gaussian of
#' # standard deviation sigma * nu^(1/4).
#' yy <- c(0.5, 1, 2, 4)
#' c(laplace = max(abs(distrib_pdf(d, yy, list(mu = 0, sigma = 1, nu = 1e-8)) -
#'                     0.5 * exp(-abs(yy)))),
#'   gaussian = max(abs(distrib_pdf(d, yy, list(mu = 0, sigma = 1, nu = 1e8)) -
#'                      dnorm(yy, 0, 100))))
#'
#' # The location score is bounded, where a Gaussian's grows without bound.
#' rr <- c(1, 4, 16, 64)
#' rbind(residual = rr,
#'       pseudohuber = distrib_gradient(d, 0.4 + rr, th)$mu,
#'       gaussian = rr / 1.2^2)
#'
#' # The expected information is a quadrature here, and the family says so.
#' expected_hessian_exact(d)
#'
#' # Fitting recovers the location and the scale; the shape is the hardest of
#' # the three to pin down from a moderate sample.
#' set.seed(5)
#' z <- distrib_rng(d, 600, list(mu = 1, sigma = 2, nu = 2))
#' coef(fit_distrib(d, z))
#'
#' @seealso
#' [laplace_distrib()] and [gaussian1_distrib()] for the two limits;
#' [student_t1_distrib()] for the other robust three-parameter family here,
#' whose score redescends to zero where this one flattens to a bound;
#' [fit_distrib()] to estimate the parameters; [check_distrib()] to validate a
#' family of your own against the same battery this one passes;
#' [PseudoHuberDistrib] for the class.
#' @export
pseudohuber_distrib <- function(link_mu = identity_link(), link_sigma = log_link(), link_nu = log_link()) {

  PseudoHuberDistrib(
    distrib_name = "pseudo huber",
    dimension = "univariate",
    bounds = c(-Inf, Inf),

    params = c("mu", "sigma", "nu"),
    params_interpretation = c(mu = "location", sigma = "scale", nu = "shape"),
    n_params = 3,

    params_bounds = list(
      mu = c(-Inf, Inf),
      sigma = c(0, Inf),
      nu = c(0, Inf)
    ),

    link_params = list(
      mu = link_mu,
      sigma = link_sigma,
      nu = link_nu
    )
  )

}
