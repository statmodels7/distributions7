#' @include distrib.R generics.R
NULL

#' @title S7 Class for Laplace Distribution
#' @name LaplaceDistrib
#'
#' @description A subclass of \code{continuous_distrib} representing the Laplace
#' (double-exponential) distribution. It is the canonical example of a distribution
#' whose log-likelihood is \strong{not differentiable} in its location parameter.
#' @inheritParams distrib
#' @return An object of class \code{LaplaceDistrib}.
#' @seealso \code{\link{laplace_distrib}}
#'
#' @section Methods:
#' Methods implemented for this class:
#'   \code{\link[=distrib_cdf.LaplaceDistrib]{distrib_cdf()}},
#'   \code{\link[=distrib_expected_hessian.LaplaceDistrib]{distrib_expected_hessian()}},
#'   \code{\link[=distrib_grad_y.LaplaceDistrib]{distrib_grad_y()}},
#'   \code{\link[=distrib_gradient.LaplaceDistrib]{distrib_gradient()}},
#'   \code{\link[=distrib_hess_y.LaplaceDistrib]{distrib_hess_y()}},
#'   \code{\link[=distrib_hessian.LaplaceDistrib]{distrib_hessian()}},
#'   \code{\link[=distrib_deriv3.LaplaceDistrib]{distrib_deriv3()}},
#'   \code{\link[=distrib_deriv4.LaplaceDistrib]{distrib_deriv4()}},
#'   \code{\link[=distrib_pdf.LaplaceDistrib]{distrib_pdf()}},
#'   \code{\link[=distrib_quantile.LaplaceDistrib]{distrib_quantile()}},
#'   \code{\link[=distrib_rng.LaplaceDistrib]{distrib_rng()}},
#'   \code{\link[=kurtosis]{kurtosis()}},
#'   \code{\link[=mean.distrib]{mean()}},
#'   \code{\link[=skewness]{skewness()}},
#'   \code{\link[=variance]{variance()}}
#'
#' Everything else is inherited from \code{\link{continuous_distrib}}.
LaplaceDistrib <- S7::new_class("LaplaceDistrib", parent = continuous_distrib)

# --- S7 METHODS IMPLEMENTATION ---

#' @title Laplace Probability Density Function
#' @name distrib_pdf.LaplaceDistrib
#' @description
#' Computes the probability density function for the Laplace distribution:
#' \deqn{f(y; \mu, \sigma) = \dfrac{1}{2\sigma} \exp\left(-\dfrac{|y-\mu|}{\sigma}\right)}
#' @param distrib A \code{LaplaceDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameters \code{mu} and \code{sigma}.
#' @param log Logical; if \code{TRUE}, returns the log-density.
#' @return A numeric vector of density values.
#' @seealso \code{\link{laplace_distrib}}
S7::method(distrib_pdf, LaplaceDistrib) <- function(distrib, y, theta, log = FALSE) {
  mu <- theta[[1]]
  b <- theta[[2]]
  log_d <- -log(2 * b) - abs(y - mu) / b
  if (log) log_d else exp(log_d)
}

#' @title Laplace Cumulative Distribution Function
#' @name distrib_cdf.LaplaceDistrib
#' @description
#' Computes the cumulative distribution function for the Laplace distribution:
#' \deqn{F(q; \mu, \sigma) = \begin{cases} \dfrac{1}{2}\exp\left(\dfrac{q-\mu}{\sigma}\right) & q < \mu \\ 1 - \dfrac{1}{2}\exp\left(-\dfrac{q-\mu}{\sigma}\right) & q \ge \mu \end{cases}}
#' @param distrib A \code{LaplaceDistrib} object.
#' @param q A numeric vector of quantiles.
#' @param theta A list containing the parameters \code{mu} and \code{sigma}.
#' @param lower.tail Logical; if \code{TRUE} (default), probabilities are \eqn{P(Y \le q)}, otherwise \eqn{P(Y > q)}.
#' @param log.p Logical; if \code{TRUE}, probabilities \eqn{p} are given as \eqn{\log(p)}.
#' @return A numeric vector of cumulative probabilities.
#' @seealso \code{\link{laplace_distrib}}
S7::method(distrib_cdf, LaplaceDistrib) <- function(distrib, q, theta, lower.tail = TRUE, log.p = FALSE) {
  mu <- theta[[1]]
  b <- theta[[2]]
  res <- ifelse(q < mu, 0.5 * exp((q - mu) / b), 1 - 0.5 * exp(-(q - mu) / b))
  if (!lower.tail) res <- 1 - res
  if (log.p) log(res) else res
}

#' @title Laplace Quantile Function
#' @name distrib_quantile.LaplaceDistrib
#' @description
#' Computes the quantile function (inverse CDF) for the Laplace distribution:
#' \deqn{Q(p; \mu, \sigma) = \mu - \sigma\,\mathrm{sign}(p - \tfrac{1}{2})\,\log\left(1 - 2\left|p - \tfrac{1}{2}\right|\right)}
#' @param distrib A \code{LaplaceDistrib} object.
#' @param p A numeric vector of probabilities.
#' @param theta A list containing the parameters \code{mu} and \code{sigma}.
#' @param lower.tail Logical; if \code{TRUE} (default), probabilities are \eqn{P(Y \le p)}, otherwise \eqn{P(Y > p)}.
#' @param log.p Logical; if \code{TRUE}, probabilities \eqn{p} are given as \eqn{\log(p)}.
#' @return A numeric vector of quantiles.
#' @seealso \code{\link{laplace_distrib}}
S7::method(distrib_quantile, LaplaceDistrib) <- function(distrib, p, theta, lower.tail = TRUE, log.p = FALSE) {
  mu <- theta[[1]]
  b <- theta[[2]]
  if (log.p) p <- exp(p)
  if (!lower.tail) p <- 1 - p
  mu - b * sign(p - 0.5) * log(1 - 2 * abs(p - 0.5))
}

#' @title Laplace Random Number Generator
#' @name distrib_rng.LaplaceDistrib
#' @description
#' Generates random numbers from the Laplace distribution by inverse-transform
#' sampling.
#' @param distrib A \code{LaplaceDistrib} object.
#' @param n Number of observations to generate.
#' @param theta A list containing the parameters \code{mu} and \code{sigma}.
#' @return A numeric vector of random draws.
#' @seealso \code{\link{laplace_distrib}}
S7::method(distrib_rng, LaplaceDistrib) <- function(distrib, n, theta) {
  distrib_quantile(distrib, stats::runif(n), theta)
}

#' @title Laplace Analytical Gradient
#' @name distrib_gradient.LaplaceDistrib
#' @description
#' Computes the analytical gradient of the Laplace log-density. Note that the
#' derivative with respect to \eqn{\mu} exists only almost everywhere (there is a
#' kink at \eqn{y = \mu}, a set of probability zero); the subgradient value 0 is
#' returned there.
#'
#' \deqn{\dfrac{\partial \ell}{\partial \mu} = \dfrac{\mathrm{sign}(y-\mu)}{\sigma}}
#' \deqn{\dfrac{\partial \ell}{\partial \sigma} = \dfrac{1}{\sigma}\left(\dfrac{|y-\mu|}{\sigma} - 1\right)}
#'
#' @param distrib A \code{LaplaceDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameters \code{mu} and \code{sigma}.
#' @return A list containing the vectors of first derivatives.
#' @seealso \code{\link{laplace_distrib}}
S7::method(distrib_gradient, LaplaceDistrib) <- function(distrib, y, theta, scale = c("parameter", "link"), ...) {
  mu <- theta[[1]]
  b <- theta[[2]]
  r <- y - mu
  list(
    mu = sign(r) / b,
    sigma = (abs(r) / b - 1) / b
  )
}

#' @title Laplace Analytical Observed Hessian
#' @name distrib_hessian.LaplaceDistrib
#' @description
#' Computes the analytical observed Hessian of the Laplace log-density. Because the
#' log-density is piecewise linear in \eqn{\mu}, the second derivative with respect
#' to \eqn{\mu} is \strong{zero} almost everywhere:
#'
#' \deqn{\dfrac{\partial^2 \ell}{\partial \mu^2} = 0, \qquad
#'       \dfrac{\partial^2 \ell}{\partial \mu \partial \sigma} = -\dfrac{\mathrm{sign}(y-\mu)}{\sigma^2}, \qquad
#'       \dfrac{\partial^2 \ell}{\partial \sigma^2} = \dfrac{\sigma - 2|y-\mu|}{\sigma^3}}
#'
#' The degeneracy of \eqn{\partial^2 \ell / \partial \mu^2} means Newton-Raphson
#' cannot update \eqn{\mu}; use \code{\link{distrib_expected_hessian}} (Fisher
#' scoring), which supplies the correct information \eqn{1/\sigma^2} for \eqn{\mu}.
#'
#' @param distrib A \code{LaplaceDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameters \code{mu} and \code{sigma}.
#' @return A list containing the vectors of second derivatives.
#' @seealso \code{\link{laplace_distrib}}
S7::method(distrib_hessian, LaplaceDistrib) <- function(distrib, y, theta, scale = c("parameter", "link"), ...) {
  mu <- theta[[1]]
  b <- theta[[2]]
  r <- y - mu
  n <- length(y)
  list(
    mu_mu = rep(0, n),
    sigma_sigma = (b - 2 * abs(r)) / b^3,
    mu_sigma = -sign(r) / b^2
  )
}

#' @title Laplace Analytical Expected Hessian (Fisher Information)
#' @name distrib_expected_hessian.LaplaceDistrib
#' @description
#' Computes the expected Hessian (negative Fisher information) of the Laplace
#' log-density. Because the log-likelihood is not differentiable in \eqn{\mu}, the
#' second Bartlett identity fails: \eqn{\mathbb{E}[\partial^2 \ell / \partial \mu^2] = 0},
#' yet the Fisher information for \eqn{\mu} is \eqn{1/\sigma^2}. The expected Hessian is
#' therefore defined here from the variance of the score, which is what the Fisher
#' information \emph{is} whenever the score exists, whether or not the identity
#' relating it to \eqn{-\mathbb{E}[H]} holds:
#'
#' \deqn{\mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial \mu^2}\right] = -\dfrac{1}{\sigma^2}, \qquad
#'       \mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial \mu \partial \sigma}\right] = 0, \qquad
#'       \mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial \sigma^2}\right] = -\dfrac{1}{\sigma^2}}
#'
#' @param distrib A \code{LaplaceDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameters \code{mu} and \code{sigma}.
#' @return A list containing the vectors of expected second derivatives.
#' @seealso \code{\link{laplace_distrib}}
S7::method(distrib_expected_hessian, LaplaceDistrib) <- function(distrib, y, theta, scale = c("parameter", "link"), approx = c("bartlett", "integrate", "mc", "opg"), nsim = 10000, ...) {
  b <- theta[[2]]
  n <- length(y)
  list(
    mu_mu = rep(-1 / b^2, length.out = n),
    sigma_sigma = rep(-1 / b^2, length.out = n),
    mu_sigma = rep(0, n)
  )
}

#' @title Laplace Analytical Third-Order Derivatives
#' @name distrib_deriv3.LaplaceDistrib
#' @description
#' Closed-form third-order derivatives of the Laplace log-density, almost
#' everywhere (observed, or expected when \code{expected = TRUE}). With
#' \eqn{s = \mathrm{sign}(y-\mu)} and \eqn{a = \lvert y-\mu \rvert}, the only
#' non-zero components are
#' \eqn{\ell^{(\mu\sigma\sigma)} = 2s/\sigma^3} and
#' \eqn{\ell^{(\sigma\sigma\sigma)} = -2/\sigma^3 + 6a/\sigma^4}; the kink at \eqn{y = \mu} is the same
#' one the gradient carries, and \code{params_smooth} records it.
#' @param distrib A \code{LaplaceDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameters \code{mu} and \code{sigma}.
#' @param expected Logical; if \code{TRUE}, returns the expected third derivatives.
#' @return A named list of third-derivative component vectors.
#' @seealso \code{\link{laplace_distrib}}
S7::method(distrib_deriv3, LaplaceDistrib) <- function(distrib, y, theta, expected = FALSE, scale = c("parameter", "link"), approx = c("integrate", "bartlett", "mc", "opg"), nsim = 10000, ...) {
  if (expected) laplace_deriv3_expected_cpp(y, theta[[1]], theta[[2]])
  else laplace_deriv3_cpp(y, theta[[1]], theta[[2]])
}

#' @title Laplace Analytical Fourth-Order Derivatives
#' @name distrib_deriv4.LaplaceDistrib
#' @description
#' Closed-form fourth-order derivatives of the Laplace log-density, almost
#' everywhere (observed, or expected when \code{expected = TRUE}), in the
#' notation of \code{\link{distrib_deriv3.LaplaceDistrib}}: the non-zero
#' components are \eqn{\ell^{(\mu\sigma\sigma\sigma)} = -6s/\sigma^4} and
#' \eqn{\ell^{(\sigma\sigma\sigma\sigma)} = 6/\sigma^4 - 24a/\sigma^5}.
#' @param distrib A \code{LaplaceDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameters \code{mu} and \code{sigma}.
#' @param expected Logical; if \code{TRUE}, returns the expected fourth derivatives.
#' @return A named list of fourth-derivative component vectors.
#' @seealso \code{\link{laplace_distrib}}
S7::method(distrib_deriv4, LaplaceDistrib) <- function(distrib, y, theta, expected = FALSE, scale = c("parameter", "link"), approx = c("integrate", "bartlett", "mc", "opg"), nsim = 10000, ...) {
  if (expected) laplace_deriv4_expected_cpp(y, theta[[1]], theta[[2]])
  else laplace_deriv4_cpp(y, theta[[1]], theta[[2]])
}

#' @title Laplace Response Derivatives
#' @name distrib_grad_y.LaplaceDistrib
#' @description
#' Closed-form derivative of the Laplace log-density with respect to the response,
#' \eqn{\partial \ell / \partial y = -\mathrm{sign}(y-\mu)/\sigma} (the second derivative
#' is 0 almost everywhere). The analytic form is provided because finite differences
#' would be inaccurate across the kink at \eqn{y = \mu}.
#' @param distrib A \code{LaplaceDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameters \code{mu} and \code{sigma}.
#' @return A numeric vector.
#' @seealso \code{\link{laplace_distrib}}
S7::method(distrib_grad_y, LaplaceDistrib) <- function(distrib, y, theta) {
  -sign(y - theta[[1]]) / theta[[2]]
}

#' @title Laplace Response Second Derivative
#' @name distrib_hess_y.LaplaceDistrib
#' @description Closed-form \eqn{\partial^2 \ell / \partial y^2 = 0} (almost everywhere).
#' @param distrib A \code{LaplaceDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing the parameters \code{mu} and \code{sigma}.
#' @return A numeric vector of zeros.
#' @seealso \code{\link{laplace_distrib}}
S7::method(distrib_hess_y, LaplaceDistrib) <- function(distrib, y, theta) {
  rep(0, length.out = length(y))
}

# --- CONSTRUCTOR WRAPPER ---

#' Laplace Distribution Object
#'
#' @description
#' Creates a distribution object for the Laplace (double-exponential) distribution,
#' parameterized by location (\eqn{\mu}) and scale (\eqn{\sigma}). This is the reference
#' example of a distribution whose log-likelihood is not differentiable in a
#' parameter (\eqn{\mu}); see \strong{Details} for how the package handles this.
#'
#' @param link_mu A link function object for the location parameter \eqn{\mu}.
#'   Defaults to \code{\link[linkfunctions7]{identity_link}}.
#' @param link_sigma A link function object for the scale parameter \eqn{\sigma}.
#'   Defaults to \code{\link[linkfunctions7]{log_link}} to ensure positivity.
#'
#' @details
#' The Laplace (double-exponential) distribution has location \eqn{\mu} and scale
#' \eqn{\sigma}.
#'
#' \strong{Probability density function:}
#' \deqn{f(y; \mu, \sigma) = \dfrac{1}{2\sigma} \exp\left(-\dfrac{|y-\mu|}{\sigma}\right)}
#'
#' \strong{Cumulative distribution function:}
#' \deqn{F(q; \mu, \sigma) = \begin{cases} \tfrac{1}{2}\exp\left(\tfrac{q-\mu}{\sigma}\right) & q < \mu \\ 1 - \tfrac{1}{2}\exp\left(-\tfrac{q-\mu}{\sigma}\right) & q \ge \mu \end{cases}}
#'
#' \strong{Quantile function:}
#' \deqn{Q(p; \mu, \sigma) = \mu - \sigma\,\mathrm{sign}(p - \tfrac{1}{2})\,\log\left(1 - 2\left|p - \tfrac{1}{2}\right|\right)}
#'
#' \strong{Score} (defined almost everywhere):
#' \deqn{\dfrac{\partial \ell}{\partial \mu} = \dfrac{\mathrm{sign}(y-\mu)}{\sigma}, \qquad
#'       \dfrac{\partial \ell}{\partial \sigma} = \dfrac{1}{\sigma}\left(\dfrac{|y-\mu|}{\sigma} - 1\right)}
#'
#' \strong{Observed Hessian:}
#' \deqn{\dfrac{\partial^2 \ell}{\partial \mu^2} = 0, \quad
#'       \dfrac{\partial^2 \ell}{\partial \mu\,\partial \sigma} = -\dfrac{\mathrm{sign}(y-\mu)}{\sigma^2}, \quad
#'       \dfrac{\partial^2 \ell}{\partial \sigma^2} = \dfrac{\sigma - 2|y-\mu|}{\sigma^3}}
#'
#' \strong{Expected Hessian} (Fisher information from the score variance; see below):
#' \deqn{\mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial \mu^2}\right] = -\dfrac{1}{\sigma^2}, \quad
#'       \mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial \sigma^2}\right] = -\dfrac{1}{\sigma^2}, \quad
#'       \mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial \mu\,\partial \sigma}\right] = 0}
#'
#' \strong{Moments:} mean \eqn{\mu}, variance \eqn{2\sigma^2}, skewness 0, excess kurtosis 3.
#'
#' \strong{Non-differentiability.}
#' The density has a kink at \eqn{y = \mu}, so the log-likelihood is not
#' differentiable in \eqn{\mu}. The package marks this via
#' \code{params_smooth = c(mu = FALSE, sigma = TRUE)} and handles it as follows:
#' \itemize{
#'   \item the \strong{score} (\code{\link{distrib_gradient}}) exists almost
#'     everywhere, equal to \eqn{\mathrm{sign}(y-\mu)/\sigma};
#'   \item the \strong{observed Hessian} (\code{\link{distrib_hessian}}) has
#'     \eqn{\partial^2 \ell / \partial \mu^2 = 0}, so Newton-Raphson cannot update
#'     \eqn{\mu};
#'   \item the \strong{expected Hessian} (\code{\link{distrib_expected_hessian}}) is
#'     implemented in closed form from the variance of the score, giving the correct
#'     Fisher information \eqn{1/\sigma^2} for \eqn{\mu} and making Fisher scoring the
#'     appropriate estimation method. Because the closed form exists, the
#'     \code{approx} argument is ignored for this distribution.
#' }
#'
#' \strong{Parameter Domains:}
#' \itemize{
#'   \item \eqn{\mu \in (-\infty, +\infty)}
#'   \item \eqn{\sigma \in (0, +\infty)}
#' }
#'
#' @return An S7 object of class \code{LaplaceDistrib} (inheriting from \code{continuous_distrib}).
#'
#' @importFrom linkfunctions7 identity_link log_link
#' @importFrom stats runif
#' @examples
#' d <- laplace_distrib()
#' d@params
#'
#' theta <- list(mu = 0, sigma = 1)
#' distrib_pdf(d, c(-1, 0, 1), theta)
#' distrib_gradient(d, c(-1, 0, 1), theta)
#'
#' @export
laplace_distrib <- function(link_mu = identity_link(), link_sigma = log_link()) {
  LaplaceDistrib(
    distrib_name = "laplace",
    dimension = "univariate",
    bounds = c(-Inf, Inf),

    params = c("mu", "sigma"),
    params_interpretation = c(mu = "location", sigma = "scale"),
    n_params = 2,

    params_bounds = list(
      mu = c(-Inf, Inf),
      sigma = c(0, Inf)
    ),

    link_params = list(
      mu = link_mu,
      sigma = link_sigma
    ),

    params_smooth = c(mu = FALSE, sigma = TRUE)
  )
}
