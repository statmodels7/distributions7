#' @include distrib.R generics.R numerical_functions.R
NULL

#' @title S7 Class for Pseudo-Huber Distribution
#' @name PseudoHuberDistrib
#'
#' @description A subclass of \code{continuous_distrib} representing the Pseudo-Huber distribution,
#' whose density is defined by the Pseudo-Huber loss kernel (a special case of the
#' Generalized Hyperbolic distribution).
#' @inheritParams distrib
#' @return An object of class \code{PseudoHuberDistrib}.
#' @seealso \code{\link{pseudohuber_distrib}}
#'
#' @section Methods:
#' Methods implemented for this class:
#'   \code{\link[=distrib_cdf.PseudoHuberDistrib]{distrib_cdf()}},
#'   \code{\link[=distrib_deriv3.PseudoHuberDistrib]{distrib_deriv3()}},
#'   \code{\link[=distrib_deriv4.PseudoHuberDistrib]{distrib_deriv4()}},
#'   \code{\link[=distrib_expected_hessian.PseudoHuberDistrib]{distrib_expected_hessian()}},
#'   \code{\link[=distrib_grad_y.PseudoHuberDistrib]{distrib_grad_y()}},
#'   \code{\link[=distrib_gradient.PseudoHuberDistrib]{distrib_gradient()}},
#'   \code{\link[=distrib_hess_y.PseudoHuberDistrib]{distrib_hess_y()}},
#'   \code{\link[=distrib_hessian.PseudoHuberDistrib]{distrib_hessian()}},
#'   \code{\link[=distrib_pdf.PseudoHuberDistrib]{distrib_pdf()}},
#'   \code{\link[=distrib_quantile.PseudoHuberDistrib]{distrib_quantile()}},
#'   \code{\link[=distrib_rng.PseudoHuberDistrib]{distrib_rng()}},
#'   \code{\link[=kurtosis]{kurtosis()}},
#'   \code{\link[=mean.distrib]{mean()}},
#'   \code{\link[=skewness]{skewness()}},
#'   \code{\link[=variance]{variance()}}
#'
#' Everything else is inherited from \code{\link{continuous_distrib}}.
PseudoHuberDistrib <- S7::new_class("PseudoHuberDistrib", parent = continuous_distrib)

# --- S7 METHODS IMPLEMENTATION ---

#' @title Pseudo-Huber Probability Density Function
#' @name distrib_pdf.PseudoHuberDistrib
#' @description
#' Computes the probability density function for the Pseudo-Huber distribution:
#' \deqn{f(y; \mu, \sigma, \nu) = \dfrac{1}{2 \sigma \sqrt{\nu} K_1(\sqrt{\nu})} \exp\left( - \sqrt{\nu + \left(\dfrac{y-\mu}{\sigma}\right)^2} \right)}
#' where \eqn{K_1} is the modified Bessel function of the second kind.
#'
#' @param distrib A \code{PseudoHuberDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameters \code{mu}, \code{sigma} and \code{nu}.
#' @param log Logical; if \code{TRUE}, returns the log-density.
#' @return A numeric vector of density values.
#' @seealso \code{\link{pseudohuber_distrib}}
S7::method(distrib_pdf, PseudoHuberDistrib) <- function(distrib, y, theta, log = FALSE) {
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
#' Computes the cumulative distribution function for the Pseudo-Huber distribution
#' by numerical integration of the density. The distribution is symmetric around
#' \eqn{\mu}, so the upper tail is computed by reflection.
#'
#' @param distrib A \code{PseudoHuberDistrib} object.
#' @param q A numeric vector of quantiles.
#' @param theta A list containing the parameters \code{mu}, \code{sigma} and \code{nu}.
#' @param lower.tail Logical; if \code{TRUE} (default), probabilities are \eqn{P(Y \le q)}, otherwise \eqn{P(Y > q)}.
#' @param log.p Logical; if \code{TRUE}, probabilities \eqn{p} are given as \eqn{\log(p)}.
#' @return A numeric vector of cumulative probabilities.
#' @seealso \code{\link{pseudohuber_distrib}}
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
#' Computes the quantile function for the Pseudo-Huber distribution by root-finding
#' on the numerical CDF. Symmetry around \eqn{\mu} is exploited: \eqn{Q(1/2) = \mu} and
#' \eqn{Q(p) = 2\mu - Q(1-p)} for \eqn{p > 1/2}.
#'
#' @param distrib A \code{PseudoHuberDistrib} object.
#' @param p A numeric vector of probabilities.
#' @param theta A list containing the parameters \code{mu}, \code{sigma} and \code{nu}.
#' @param lower.tail Logical; if \code{TRUE} (default), probabilities are \eqn{P(Y \le p)}, otherwise \eqn{P(Y > p)}.
#' @param log.p Logical; if \code{TRUE}, probabilities \eqn{p} are given as \eqn{\log(p)}.
#' @return A numeric vector of quantiles.
#' @seealso \code{\link{pseudohuber_distrib}}
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
#' Generates random numbers from the Pseudo-Huber distribution via inverse transform
#' sampling (numerical quantile function applied to uniform draws).
#'
#' @param distrib A \code{PseudoHuberDistrib} object.
#' @param n Number of observations to generate.
#' @param theta A list containing the parameters \code{mu}, \code{sigma} and \code{nu}.
#' @return A numeric vector of random draws.
#' @seealso \code{\link{pseudohuber_distrib}}
S7::method(distrib_rng, PseudoHuberDistrib) <- function(distrib, n, theta) {
  distrib_quantile(distrib, stats::runif(n), theta)
}

#' @title Pseudo-Huber Analytical Gradient
#' @name distrib_gradient.PseudoHuberDistrib
#' @description
#' Computes the analytical gradient of the Pseudo-Huber log-density.
#' Let \eqn{r = y - \mu} and \eqn{D = \sqrt{\nu + (r/\sigma)^2}}:
#'
#' \deqn{\dfrac{\partial \ell}{\partial \mu} = \dfrac{r}{\sigma^2 D}}
#' \deqn{\dfrac{\partial \ell}{\partial \sigma} = \dfrac{1}{\sigma} \left( \dfrac{r^2}{\sigma^2 D} - 1 \right)}
#' \deqn{\dfrac{\partial \ell}{\partial \nu} = -\dfrac{1}{2} \left[ \dfrac{1}{\nu} + \dfrac{1}{D} + \dfrac{K_1'(\sqrt{\nu})}{\sqrt{\nu}\, K_1(\sqrt{\nu})} \right]}
#'
#' @param distrib A \code{PseudoHuberDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameters \code{mu}, \code{sigma} and \code{nu}.
#' @return A list containing the vectors of first derivatives.
#' @seealso \code{\link{pseudohuber_distrib}}
S7::method(distrib_gradient, PseudoHuberDistrib) <- function(distrib, y, theta, scale = c("parameter", "link"), ...) {
  pseudohuber_gradient_cpp(y, theta[[1]], theta[[2]], theta[[3]])
}

#' @title Pseudo-Huber Analytical Observed Hessian
#' @name distrib_hessian.PseudoHuberDistrib
#' @description
#' Computes the analytical observed Hessian of the Pseudo-Huber log-density.
#' Let \eqn{r = y - \mu}, \eqn{D = \sqrt{\nu + (r/\sigma)^2}},
#' \eqn{R_1 = K_1'(\sqrt{\nu})/K_1(\sqrt{\nu})} and \eqn{R_2 = K_1''(\sqrt{\nu})/K_1(\sqrt{\nu})}:
#'
#' \deqn{\dfrac{\partial^2 \ell}{\partial \mu^2} = -\dfrac{\nu}{\sigma^2 D^3}}
#' \deqn{\dfrac{\partial^2 \ell}{\partial \sigma^2} = \dfrac{\sigma^4 - 3\sigma^2 r^2 D^{-1} + r^4 D^{-3}}{\sigma^6}}
#' \deqn{\dfrac{\partial^2 \ell}{\partial \nu^2} = \dfrac{1}{4D^3} + \dfrac{1}{2\nu^2} + \dfrac{1}{4}\left(\dfrac{R_1}{\nu^{3/2}} + \dfrac{R_1^2}{\nu} - \dfrac{R_2}{\nu}\right)}
#' \deqn{\dfrac{\partial^2 \ell}{\partial \mu \partial \sigma} = \dfrac{-2\nu\sigma^2 r - r^3}{\sigma^2(\nu\sigma^2 + r^2)^{3/2}}}
#' \deqn{\dfrac{\partial^2 \ell}{\partial \mu \partial \nu} = -\dfrac{r}{2\sigma^2 D^3}}
#' \deqn{\dfrac{\partial^2 \ell}{\partial \sigma \partial \nu} = -\dfrac{r^2}{2\sigma^3 D^3}}
#'
#' @param distrib A \code{PseudoHuberDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameters \code{mu}, \code{sigma} and \code{nu}.
#' @return A list containing the vectors of second derivatives.
#' @seealso \code{\link{pseudohuber_distrib}}
S7::method(distrib_hessian, PseudoHuberDistrib) <- function(distrib, y, theta, scale = c("parameter", "link"), ...) {
  pseudohuber_hessian_cpp(y, theta[[1]], theta[[2]], theta[[3]])
}

#' @title Pseudo-Huber Analytical Expected Hessian
#' @name distrib_expected_hessian.PseudoHuberDistrib
#' @description
#' Computes the expected Hessian of the Pseudo-Huber log-density.
#' No closed form exists, so the non-zero components are evaluated by numerical
#' integration of the observed Hessian against the density (via \code{\link{expectation}}).
#' By symmetry \eqn{\mathbb{E}[r] = \mathbb{E}[r^3] = 0}, hence the \eqn{\mu\sigma} and
#' \eqn{\mu\nu} components are exactly 0.
#'
#' @param distrib A \code{PseudoHuberDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameters \code{mu}, \code{sigma} and \code{nu}.
#' @return A list containing the vectors of expected second derivatives.
#' @seealso \code{\link{pseudohuber_distrib}}
S7::method(distrib_expected_hessian, PseudoHuberDistrib) <- function(distrib, y, theta, scale = c("parameter", "link"), approx = c("bartlett", "integrate", "mc", "opg"), nsim = 10000, ...) {
  n <- length(y)
  out <- expected_derivative(distrib, y, theta, order = 2L,
                             approx = match.arg(approx), nsim = nsim)
  # exact by symmetry: E[r] = E[r^3] = 0
  out$mu_sigma <- rep(0, n)
  out$mu_nu <- rep(0, n)
  out
}

#' @title Pseudo-Huber Analytical Third-Order Derivatives
#' @name distrib_deriv3.PseudoHuberDistrib
#' @description
#' Closed-form observed third-order derivatives of the Pseudo-Huber log-density.
#' Bessel functions enter only through the pure-\eqn{\nu} component; the
#' exponentially scaled forms are used so that large \eqn{\nu} does not overflow.
#' The expected third derivatives have no closed form, so \code{expected = TRUE}
#' is handled by the strategies in \code{\link{expected_derivative_methods}}.
#' @param distrib A \code{PseudoHuberDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameters \code{mu}, \code{sigma} and \code{nu}.
#' @param expected Logical; if \code{TRUE}, returns the (approximated) expected derivatives.
#' @param approx,nsim Passed to \code{\link{expected_derivative_methods}} when \code{expected = TRUE}.
#' @return A named list of third-derivative component vectors.
#' @seealso \code{\link{pseudohuber_distrib}}
S7::method(distrib_deriv3, PseudoHuberDistrib) <- function(distrib, y, theta, expected = FALSE, scale = c("parameter", "link"), approx = c("integrate", "bartlett", "mc", "opg"), nsim = 10000, ...) {
  if (expected) {
    expected_derivative(distrib, y, theta, order = 3L,
                        approx = match.arg(approx), nsim = nsim)
  } else {
    pseudohuber_deriv3_cpp(y, theta[[1]], theta[[2]], theta[[3]])
  }
}

#' @title Pseudo-Huber Analytical Fourth-Order Derivatives
#' @name distrib_deriv4.PseudoHuberDistrib
#' @description
#' Closed-form observed fourth-order derivatives of the Pseudo-Huber log-density
#' (see \code{\link{distrib_deriv3.PseudoHuberDistrib}} for the Bessel handling).
#' The expected fourth derivatives have no closed form and are handled by the
#' strategies in \code{\link{expected_derivative_methods}}.
#' @param distrib A \code{PseudoHuberDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameters \code{mu}, \code{sigma} and \code{nu}.
#' @param expected Logical; if \code{TRUE}, returns the (approximated) expected derivatives.
#' @param approx,nsim Passed to \code{\link{expected_derivative_methods}} when \code{expected = TRUE}.
#' @return A named list of fourth-derivative component vectors.
#' @seealso \code{\link{pseudohuber_distrib}}
S7::method(distrib_deriv4, PseudoHuberDistrib) <- function(distrib, y, theta, expected = FALSE, scale = c("parameter", "link"), approx = c("integrate", "bartlett", "mc", "opg"), nsim = 10000, ...) {
  if (expected) {
    expected_derivative(distrib, y, theta, order = 4L,
                        approx = match.arg(approx), nsim = nsim)
  } else {
    pseudohuber_deriv4_cpp(y, theta[[1]], theta[[2]], theta[[3]])
  }
}

#' @title Pseudo-Huber Response Derivatives
#' @name distrib_grad_y.PseudoHuberDistrib
#' @description
#' Closed-form derivatives of the Pseudo-Huber log-density with respect to the
#' response. Let \eqn{r = y - \mu} and \eqn{D = \sqrt{\nu + (r/\sigma)^2}}:
#' \eqn{\partial \ell / \partial y = -r/(\sigma^2 D)} and
#' \eqn{\partial^2 \ell / \partial y^2 = -\nu/(\sigma^2 D^3)}.
#' @param distrib A \code{PseudoHuberDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameters \code{mu}, \code{sigma} and \code{nu}.
#' @return A numeric vector.
#' @seealso \code{\link{pseudohuber_distrib}}
S7::method(distrib_grad_y, PseudoHuberDistrib) <- function(distrib, y, theta) {
  mu <- theta[[1]]; sigma <- theta[[2]]; nu <- theta[[3]]
  r <- y - mu
  sigma2 <- sigma^2
  D <- sqrt(nu + r^2 / sigma2)
  -r / (sigma2 * D)
}

#' @title Pseudo-Huber Response Second Derivative
#' @name distrib_hess_y.PseudoHuberDistrib
#' @description Closed-form \eqn{\partial^2 \ell / \partial y^2 = -\nu/(\sigma^2 D^3)}, \eqn{D = \sqrt{\nu + ((y-\mu)/\sigma)^2}}.
#' @param distrib A \code{PseudoHuberDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameters \code{mu}, \code{sigma} and \code{nu}.
#' @return A numeric vector.
#' @seealso \code{\link{pseudohuber_distrib}}
S7::method(distrib_hess_y, PseudoHuberDistrib) <- function(distrib, y, theta) {
  mu <- theta[[1]]; sigma <- theta[[2]]; nu <- theta[[3]]
  r <- y - mu
  sigma2 <- sigma^2
  D <- sqrt(nu + r^2 / sigma2)
  -nu / (sigma2 * D^3)
}

# --- CONSTRUCTOR WRAPPER ---

#' Pseudo-Huber Distribution Object (Location-Scale Parameterization)
#'
#' @description
#' Creates a distribution object for the Pseudo-Huber distribution, whose density
#' corresponds to the exponential of the negative Pseudo-Huber loss. It is a special
#' case of the Generalized Hyperbolic distribution.
#'
#' @param link_mu A link function object for the location parameter \eqn{\mu}.
#'   Defaults to \code{\link[linkfunctions7]{identity_link}}.
#' @param link_sigma A link function object for the scale parameter \eqn{\sigma}.
#'   Defaults to \code{\link[linkfunctions7]{log_link}} to ensure positivity.
#' @param link_nu A link function object for the shape parameter \eqn{\nu}.
#'   Defaults to \code{\link[linkfunctions7]{log_link}} to ensure positivity.
#'
#' @details
#' The probability density function is:
#' \deqn{f(y; \mu, \sigma, \nu) = \dfrac{1}{2 \sigma \sqrt{\nu} K_1(\sqrt{\nu})} \exp\left( - \sqrt{\nu + \left(\dfrac{y-\mu}{\sigma}\right)^2} \right)}
#' where \eqn{K_1} is the modified Bessel function of the second kind.
#'
#' \strong{Moments:}
#' \itemize{
#'   \item Expected value: \eqn{\mathbb{E}(Y) = \mu}
#'   \item Variance: \eqn{\mathbb{V}(Y) = \sigma^2 \sqrt{\nu}\, \dfrac{K_2(\sqrt{\nu})}{K_1(\sqrt{\nu})}}
#'   \item Skewness: 0 (symmetric)
#'   \item Excess kurtosis: \eqn{3 \dfrac{K_3(\sqrt{\nu}) K_1(\sqrt{\nu})}{K_2(\sqrt{\nu})^2} - 3}
#' }
#'
#' \strong{Score} (with \eqn{r = y-\mu} and \eqn{D = \sqrt{\nu + (r/\sigma)^2}};
#' \eqn{K_1'} the derivative of \eqn{K_1}):
#' \deqn{\dfrac{\partial \ell}{\partial \mu} = \dfrac{r}{\sigma^2 D}, \qquad
#'       \dfrac{\partial \ell}{\partial \sigma} = \dfrac{1}{\sigma}\left(\dfrac{r^2}{\sigma^2 D} - 1\right)}
#' \deqn{\dfrac{\partial \ell}{\partial \nu} = -\dfrac{1}{2}\left[\dfrac{1}{\nu} + \dfrac{1}{D} + \dfrac{K_1'(\sqrt{\nu})}{\sqrt{\nu}\,K_1(\sqrt{\nu})}\right]}
#' The observed Hessian is available in closed form via
#' \code{\link{distrib_hessian.PseudoHuberDistrib}}; the expected Hessian has no
#' closed form and is obtained by numerical integration.
#'
#' \strong{Parameter Domains:}
#' \itemize{
#'   \item \eqn{\mu \in (-\infty, +\infty)}
#'   \item \eqn{\sigma \in (0, +\infty)}
#'   \item \eqn{\nu \in (0, +\infty)}
#' }
#'
#' \strong{Note:} The CDF, quantile function and RNG have no closed form and rely on
#' numerical integration / root-finding, so they are slower than for the other
#' distributions in the package. Response derivatives (\code{\link{distrib_grad_y}},
#' \code{\link{distrib_hess_y}}) are available in closed form; third- and
#' fourth-order parameter derivatives use the numerical fallback.
#'
#' @seealso
#' \itemize{
#'   \item \code{\link{distrib_pdf.PseudoHuberDistrib}} for the density function.
#'   \item \code{\link{distrib_cdf.PseudoHuberDistrib}} for the cumulative distribution function.
#'   \item \code{\link{distrib_quantile.PseudoHuberDistrib}} for the quantile function.
#'   \item \code{\link{distrib_rng.PseudoHuberDistrib}} for random number generation.
#'   \item \code{\link{distrib_gradient.PseudoHuberDistrib}} for the analytical gradient.
#'   \item \code{\link{distrib_hessian.PseudoHuberDistrib}} for the analytical observed Hessian.
#'   \item \code{\link{distrib_expected_hessian.PseudoHuberDistrib}} for the expected Hessian.
#' }
#'
#' @return An S7 object of class \code{PseudoHuberDistrib} (inheriting from \code{continuous_distrib}) representing the Pseudo-Huber distribution.
#'
#' @importFrom linkfunctions7 identity_link log_link
#' @importFrom stats runif uniroot
#' @examples
#' d <- pseudohuber_distrib()
#' d@params
#'
#' theta <- list(mu = 0, sigma = 1, nu = 1)
#' distrib_pdf(d, c(-1, 0, 1), theta)
#' distrib_gradient(d, c(-1, 0, 1), theta)
#'
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
