#' @include distrib.R generics.R
NULL

#' @title S7 Class for the Generalized Pareto Distribution
#' @name GPDDistrib
#'
#' @description A subclass of \code{continuous_distrib} representing the
#'   generalized Pareto distribution in its scale and shape.
#' @inheritParams distrib
#' @return An object of class \code{GPDDistrib}.
#' @seealso \code{\link{gpd_distrib}}
#'
#' @section Methods:
#' Methods implemented for this class:
#'   \code{\link[=distrib_cdf.GPDDistrib]{distrib_cdf()}},
#'   \code{\link[=distrib_expected_hessian.GPDDistrib]{distrib_expected_hessian()}},
#'   \code{\link[=distrib_gradient.GPDDistrib]{distrib_gradient()}},
#'   \code{\link[=distrib_hessian.GPDDistrib]{distrib_hessian()}},
#'   \code{\link[=distrib_pdf.GPDDistrib]{distrib_pdf()}},
#'   \code{\link[=distrib_quantile.GPDDistrib]{distrib_quantile()}},
#'   \code{\link[=distrib_rng.GPDDistrib]{distrib_rng()}}
#'
#' Everything else is inherited from \code{\link{continuous_distrib}}.
GPDDistrib <- S7::new_class("GPDDistrib", parent = continuous_distrib)

#' The Upper Endpoint of a Generalized Pareto
#'
#' @description
#' \eqn{-\sigma/\xi} when \eqn{\xi < 0}, and infinity otherwise.
#'
#' @details
#' The endpoint depends on the parameters, which is the whole reason the
#' family needs care: for \eqn{\xi < 0} the support is bounded and moves with
#' \eqn{\sigma} and \eqn{\xi}, so the license to differentiate under the
#' integral sign is not automatic. See \code{\link{gpd_distrib}}.
#'
#' @param sigma The scale, a positive numeric vector.
#' @param xi The shape, a numeric vector.
#'
#' @return A numeric vector.
#'
#' @seealso \code{\link{gpd_distrib}}
#'
#' @keywords internal
gpd_endpoint <- function(sigma, xi) ifelse(xi < 0, -sigma / xi, Inf)

# --- S7 METHODS IMPLEMENTATION ---

#' @title Generalized Pareto Density
#' @name distrib_pdf.GPDDistrib
#' @description
#' \deqn{f(y) = \dfrac{1}{\sigma}
#'       \left(1 + \dfrac{\xi y}{\sigma}\right)^{-1/\xi - 1}}
#' with the exponential density as the limit at \eqn{\xi = 0}.
#' @param distrib A \code{GPDDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing \code{sigma} and \code{xi}.
#' @param log Logical; if \code{TRUE}, returns the log-density.
#' @return A numeric vector of density values.
#' @seealso \code{\link{gpd_distrib}}
S7::method(distrib_pdf, GPDDistrib) <- function(distrib, y, theta, log = FALSE, ..., threads = 1L) {
  out <- gpd_logpdf_cpp(y, theta[[1]], theta[[2]], threads)
  if (log) out else exp(out)
}

#' @title Generalized Pareto Distribution Function
#' @name distrib_cdf.GPDDistrib
#' @description
#' \deqn{F(q) = 1 - \left(1 + \dfrac{\xi q}{\sigma}\right)^{-1/\xi}}
#' @param distrib A \code{GPDDistrib} object.
#' @param q A numeric vector of quantiles.
#' @param theta A list containing \code{sigma} and \code{xi}.
#' @param lower.tail Logical; if \code{TRUE} (default), \eqn{P(Y \le q)}.
#' @param log.p Logical; if \code{TRUE}, probabilities are returned as logarithms.
#' @return A numeric vector of cumulative probabilities.
#' @seealso \code{\link{gpd_distrib}}
S7::method(distrib_cdf, GPDDistrib) <- function(distrib, q, theta,
                                                 lower.tail = TRUE,
                                                 log.p = FALSE) {
  # ifelse() returns a result the length of its TEST, so a scalar shape would
  # collapse a vector of quantiles to one number. The branch is taken with a
  # test recycled to the answer's length instead.
  s <- theta[[1]]
  x <- theta[[2]]
  z <- q / s
  n <- length(z)
  small <- rep_len(abs(x) < 1e-8, n)
  t <- 1 + rep_len(x, n) * z
  surv <- numeric(n)
  surv[small] <- exp(-z[small])
  surv[!small] <- t[!small]^(-1 / rep_len(x, n)[!small])
  surv[q < 0] <- 1
  surv[t <= 0] <- 0
  surv <- pmin(pmax(surv, 0), 1)
  p <- if (lower.tail) 1 - surv else surv
  if (log.p) log(p) else p
}

#' @title Generalized Pareto Quantile Function
#' @name distrib_quantile.GPDDistrib
#' @description
#' \deqn{Q(p) = \dfrac{\sigma}{\xi}\left((1-p)^{-\xi} - 1\right)}
#' with \eqn{-\sigma\log(1-p)} at \eqn{\xi = 0}.
#' @param distrib A \code{GPDDistrib} object.
#' @param p A numeric vector of probabilities.
#' @param theta A list containing \code{sigma} and \code{xi}.
#' @param lower.tail Logical; if \code{TRUE} (default), \code{p} is \eqn{P(Y \le q)}.
#' @param log.p Logical; if \code{TRUE}, \code{p} is given as its logarithm.
#' @return A numeric vector of quantiles.
#' @seealso \code{\link{gpd_distrib}}
S7::method(distrib_quantile, GPDDistrib) <- function(distrib, p, theta,
                                                      lower.tail = TRUE,
                                                      log.p = FALSE) {
  if (log.p) p <- exp(p)
  if (!lower.tail) p <- 1 - p
  s <- theta[[1]]
  x <- theta[[2]]
  n <- length(p)
  sv <- rep_len(s, n)
  xv <- rep_len(x, n)
  small <- abs(xv) < 1e-8
  out <- numeric(n)
  out[small] <- -sv[small] * log1p(-p[small])
  out[!small] <- sv[!small] * ((1 - p[!small])^(-xv[!small]) - 1) / xv[!small]
  out
}

#' @title Generalized Pareto Random Generation
#' @name distrib_rng.GPDDistrib
#' @description Inverse transform, the quantile function being elementary.
#' @param distrib A \code{GPDDistrib} object.
#' @param n The number of draws.
#' @param theta A list containing \code{sigma} and \code{xi}.
#' @return A numeric vector of length \code{n}.
#' @seealso \code{\link{gpd_distrib}}
S7::method(distrib_rng, GPDDistrib) <- function(distrib, n, theta) {
  distrib_quantile(distrib, stats::runif(n), theta)
}

#' @title Generalized Pareto Analytical Gradient
#' @name distrib_gradient.GPDDistrib
#' @description
#' With \eqn{z = y/\sigma}, \eqn{t = 1 + \xi z} and \eqn{u = z/t},
#' \deqn{\dfrac{\partial \ell}{\partial \sigma}
#'         = \dfrac{(\xi+1)u - 1}{\sigma}, \qquad
#'       \dfrac{\partial \ell}{\partial \xi}
#'         = \dfrac{\log t}{\xi^2} - \left(1 + \dfrac{1}{\xi}\right)u}
#' the second computed through a series near \eqn{\xi = 0}, where its limit is
#' \eqn{z^2/2 - z}.
#' @param distrib A \code{GPDDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing \code{sigma} and \code{xi}.
#' @param scale Either \code{"parameter"} or \code{"link"}; handled by the generic.
#' @param ... Unused.
#' @return A named list with the \code{sigma} and \code{xi} components.
#' @seealso \code{\link{gpd_distrib}}
S7::method(distrib_gradient, GPDDistrib) <- function(distrib, y, theta,
                                                      scale = c("parameter", "link"), ..., threads = 1L) {
  gpd_gradient_cpp(y, theta[[1]], theta[[2]], threads)
}

#' @title Generalized Pareto Analytical Observed Hessian
#' @name distrib_hessian.GPDDistrib
#' @description
#' The second derivatives of the same expressions, kept short by
#' \eqn{t - \xi z = 1}, which makes \eqn{\partial u/\partial\sigma} equal to
#' \eqn{-z/(\sigma t^2)}. The pure-\eqn{\xi} component goes through a series
#' near zero, where its two singular terms cancel.
#' @param distrib A \code{GPDDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing \code{sigma} and \code{xi}.
#' @param scale Either \code{"parameter"} or \code{"link"}; handled by the generic.
#' @param ... Unused.
#' @return A named list of second-derivative components.
#' @seealso \code{\link{gpd_distrib}}
S7::method(distrib_hessian, GPDDistrib) <- function(distrib, y, theta,
                                                     scale = c("parameter", "link"), ..., threads = 1L) {
  gpd_hessian_cpp(y, theta[[1]], theta[[2]], threads)
}

#' @title Generalized Pareto Analytical Expected Hessian
#' @name distrib_expected_hessian.GPDDistrib
#' @description
#' The closed form of Smith (1985), valid for \eqn{\xi > -1/2}:
#' \deqn{\mathbb{E}\left[\dfrac{\partial^2\ell}{\partial\sigma^2}\right]
#'         = \dfrac{-1}{(1+2\xi)\sigma^2}, \qquad
#'       \mathbb{E}\left[\dfrac{\partial^2\ell}{\partial\sigma\partial\xi}\right]
#'         = \dfrac{-1}{(1+2\xi)\sigma(1+\xi)}, \qquad
#'       \mathbb{E}\left[\dfrac{\partial^2\ell}{\partial\xi^2}\right]
#'         = \dfrac{-2}{(1+2\xi)(1+\xi)}}
#' At \eqn{\xi \le -1/2} the information does not exist and \code{NA} is
#' returned rather than a number; see \code{\link{gpd_distrib}}.
#' @param distrib A \code{GPDDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing \code{sigma} and \code{xi}.
#' @param scale Either \code{"parameter"} or \code{"link"}; handled by the generic.
#' @param approx Ignored; the expectation is closed form.
#' @param nsim Ignored.
#' @param ... Unused.
#' @return A named list of expected second-derivative components.
#' @seealso \code{\link{gpd_distrib}}
S7::method(distrib_expected_hessian, GPDDistrib) <- function(distrib, y, theta,
                                                              scale = c("parameter", "link"),
                                                              approx = c("bartlett", "integrate", "mc", "opg"),
                                                              nsim = 10000, ..., threads = 1L) {
  gpd_expected_hessian_cpp(y, theta[[1]], theta[[2]], threads)
}

# --- CONSTRUCTOR WRAPPER ---

#' Generalized Pareto Distribution Object
#'
#' @description
#' Creates a distribution object for the generalized Pareto distribution,
#' parametrized by a scale \eqn{\sigma} and a shape \eqn{\xi}.
#'
#' @param link_sigma A link function object for \eqn{\sigma}. Defaults to
#'   \code{\link[linkfunctions7]{log_link}} to ensure positivity.
#' @param link_xi A link function object for \eqn{\xi}. Defaults to
#'   \code{\link[linkfunctions7]{identity_link}}, the shape being free to take
#'   either sign.
#'
#' @details
#' The family of exceedances over a high threshold, and the natural companion
#' of \code{\link{gumbel_distrib}} in an analysis of extremes.
#'
#' \strong{Density and distribution function:}
#' \deqn{f(y) = \dfrac{1}{\sigma}\left(1 + \dfrac{\xi y}{\sigma}\right)^{-1/\xi-1},
#'       \qquad
#'       F(q) = 1 - \left(1 + \dfrac{\xi q}{\sigma}\right)^{-1/\xi}}
#' At \eqn{\xi = 0} both reduce to the exponential, which the implementation
#' reaches through a series rather than by a special case, so the parameter may
#' pass through zero during a fit.
#'
#' \strong{It is not parametrized by its mean}, unlike most families here. The
#' mean is \eqn{\sigma/(1-\xi)} and exists only for \eqn{\xi < 1}, so a mean
#' parametrization would leave the family undescribable exactly where it is
#' most used --- the heavy-tailed regime. This is the argument that keeps
#' \code{\link{mv_sigma}} and \code{\link{variance}} apart for the multivariate
#' \eqn{t}: a parametrization must not depend on a moment that need not exist.
#'
#' \strong{The support depends on the parameters when \eqn{\xi < 0}}, being
#' \eqn{[0, -\sigma/\xi]}, and this is the first family here of which that is
#' true. What it costs is the automatic license to differentiate under the
#' integral sign, on which the Bartlett identities rest. Two things survive and
#' one does not:
#' \itemize{
#'   \item the derivatives returned are correct as derivatives of the
#'     log-density at every admissible point, whatever the sign of \eqn{\xi};
#'   \item the expected information exists and is the closed form above for
#'     \eqn{\xi > -1/2}. The condition is exactly that the integrand be
#'     integrable: near the upper endpoint the second derivative grows like
#'     \eqn{(1-u)^{-2|\xi|}} in the probability scale, which is integrable if
#'     and only if \eqn{|\xi| < 1/2};
#'   \item below \eqn{\xi = -1/2} the information does not exist,
#'     \code{\link{distrib_expected_hessian}} returns \code{NA}, and the
#'     classical asymptotics of the maximum likelihood estimator do not hold
#'     (Smith, 1985).
#' }
#' The \code{bounds} of the object are \code{c(0, Inf)} because they are fixed
#' at construction while the true endpoint moves with the parameters; the
#' density is zero beyond it, so nothing computes a wrong number, but a caller
#' reading \code{bounds} learns less than usual.
#'
#' \strong{Moments:} mean \eqn{\sigma/(1-\xi)} for \eqn{\xi < 1}, variance
#' \eqn{\sigma^2/\{(1-\xi)^2(1-2\xi)\}} for \eqn{\xi < 1/2}.
#'
#' \strong{Parameter domains:}
#' \itemize{
#'   \item \eqn{\sigma \in (0, +\infty)}
#'   \item \eqn{\xi \in (-\infty, +\infty)}
#' }
#'
#' @return An S7 object of class \code{GPDDistrib}.
#'
#' @references
#' Smith, R. L. (1985). Maximum likelihood estimation in a class of nonregular
#' cases. \emph{Biometrika} 72, 67-90.
#'
#' Davison, A. C. and Smith, R. L. (1990). Models for exceedances over high
#' thresholds. \emph{Journal of the Royal Statistical Society B} 52, 393-442.
#'
#' @seealso \code{\link{gumbel_distrib}}, \code{\link{exponential_distrib}},
#'   \code{\link{weibull1_distrib}}
#'
#' @importFrom linkfunctions7 log_link identity_link
#' @importFrom stats runif
#' @examples
#' d <- gpd_distrib()
#' d@params
#'
#' theta <- list(sigma = 1.5, xi = 0.3)
#' distrib_pdf(d, c(0.2, 1, 4), theta)
#' distrib_gradient(d, c(0.2, 1, 4), theta)
#'
#' # at xi = 0 it is the exponential, reached by a series
#' max(abs(distrib_pdf(d, c(0.2, 1, 4), list(sigma = 1.5, xi = 0)) -
#'         dexp(c(0.2, 1, 4), 1 / 1.5)))
#'
#' # the information exists only above -1/2
#' distrib_expected_hessian(d, 0, list(sigma = 1.5, xi = -0.7))
#'
#' @export
gpd_distrib <- function(link_sigma = log_link(), link_xi = identity_link()) {
  GPDDistrib(
    distrib_name = "generalized pareto", dimension = "univariate",
    bounds = c(0, Inf),
    params = c("sigma", "xi"),
    params_interpretation = c(sigma = "scale", xi = "shape"),
    n_params = 2,
    params_bounds = list(sigma = c(0, Inf), xi = c(-Inf, Inf)),
    link_params = list(sigma = link_sigma, xi = link_xi)
  )
}
