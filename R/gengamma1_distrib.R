#' @include distrib.R generics.R
NULL

#' @title S7 Class for the Generalized Gamma Distribution
#' @name GenGamma1Distrib
#'
#' @description A subclass of \code{continuous_distrib} representing the
#'   generalized gamma distribution in Stacy's three-parameter form.
#' @inheritParams distrib
#' @return An object of class \code{GenGamma1Distrib}.
#' @seealso \code{\link{gengamma1_distrib}}
#'
#' @section Methods:
#' Methods implemented for this class:
#'   \code{\link[=distrib_cdf.GenGamma1Distrib]{distrib_cdf()}},
#'   \code{\link[=distrib_expected_hessian.GenGamma1Distrib]{distrib_expected_hessian()}},
#'   \code{\link[=distrib_gradient.GenGamma1Distrib]{distrib_gradient()}},
#'   \code{\link[=distrib_hessian.GenGamma1Distrib]{distrib_hessian()}},
#'   \code{\link[=distrib_pdf.GenGamma1Distrib]{distrib_pdf()}},
#'   \code{\link[=distrib_quantile.GenGamma1Distrib]{distrib_quantile()}},
#'   \code{\link[=distrib_rng.GenGamma1Distrib]{distrib_rng()}}
#'
#' Everything else is inherited from \code{\link{continuous_distrib}}.
GenGamma1Distrib <- S7::new_class("GenGamma1Distrib", parent = continuous_distrib)

# --- S7 METHODS IMPLEMENTATION ---

#' @title Generalized Gamma Density
#' @name distrib_pdf.GenGamma1Distrib
#' @description
#' \deqn{f(y) = \dfrac{p}{a^{d}\,\Gamma(d/p)}\, y^{d-1} e^{-(y/a)^{p}},
#'       \qquad y > 0}
#' @param distrib A \code{GenGamma1Distrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing \code{a}, \code{d} and \code{p}.
#' @param log Logical; if \code{TRUE}, returns the log-density.
#' @return A numeric vector of density values.
#' @seealso \code{\link{gengamma1_distrib}}
S7::method(distrib_pdf, GenGamma1Distrib) <- function(distrib, y, theta, log = FALSE, ..., threads = 1L) {
  out <- gengamma_logpdf_cpp(y, theta[[1]], theta[[2]], theta[[3]], threads)
  if (log) out else exp(out)
}

#' @title Generalized Gamma Distribution Function
#' @name distrib_cdf.GenGamma1Distrib
#' @description
#' \eqn{F(q) = P(d/p,\, (q/a)^{p})}, the regularized lower incomplete gamma
#' function, since \eqn{(Y/a)^{p}} is Gamma with shape \eqn{d/p} and unit
#' rate.
#' @param distrib A \code{GenGamma1Distrib} object.
#' @param q A numeric vector of quantiles.
#' @param theta A list containing \code{a}, \code{d} and \code{p}.
#' @param lower.tail Logical; if \code{TRUE} (default), \eqn{P(Y \le q)}.
#' @param log.p Logical; if \code{TRUE}, probabilities are returned as logarithms.
#' @return A numeric vector of cumulative probabilities.
#' @seealso \code{\link{gengamma1_distrib}}
S7::method(distrib_cdf, GenGamma1Distrib) <- function(distrib, q, theta,
                                                      lower.tail = TRUE,
                                                      log.p = FALSE) {
  a <- theta[[1]]; d <- theta[[2]]; p <- theta[[3]]
  w <- pmax(q, 0)^p / a^p
  stats::pgamma(w, shape = d / p, rate = 1,
                lower.tail = lower.tail, log.p = log.p)
}

#' @title Generalized Gamma Quantile Function
#' @name distrib_quantile.GenGamma1Distrib
#' @description
#' \eqn{Q(u) = a\,\{Q_{\Gamma}(u; d/p)\}^{1/p}}, inverting the same
#' representation the distribution function uses.
#' @param distrib A \code{GenGamma1Distrib} object.
#' @param p A numeric vector of probabilities.
#' @param theta A list containing \code{a}, \code{d} and \code{p}.
#' @param lower.tail Logical; if \code{TRUE} (default), \code{p} is \eqn{P(Y \le q)}.
#' @param log.p Logical; if \code{TRUE}, \code{p} is given as its logarithm.
#' @return A numeric vector of quantiles.
#' @seealso \code{\link{gengamma1_distrib}}
S7::method(distrib_quantile, GenGamma1Distrib) <- function(distrib, p, theta,
                                                           lower.tail = TRUE,
                                                           log.p = FALSE) {
  a <- theta[[1]]; d <- theta[[2]]; pw <- theta[[3]]
  g <- stats::qgamma(p, shape = d / pw, rate = 1,
                     lower.tail = lower.tail, log.p = log.p)
  a * g^(1 / pw)
}

#' @title Generalized Gamma Random Generation
#' @name distrib_rng.GenGamma1Distrib
#' @description
#' A Gamma draw raised to the power \eqn{1/p} and scaled, which is the
#' representation the family is defined by rather than an approximation of it.
#' @param distrib A \code{GenGamma1Distrib} object.
#' @param n The number of draws.
#' @param theta A list containing \code{a}, \code{d} and \code{p}.
#' @return A numeric vector of length \code{n}.
#' @seealso \code{\link{gengamma1_distrib}}
S7::method(distrib_rng, GenGamma1Distrib) <- function(distrib, n, theta) {
  a <- theta[[1]]; d <- theta[[2]]; p <- theta[[3]]
  a * stats::rgamma(n, shape = d / p, rate = 1)^(1 / p)
}

#' @title Generalized Gamma Analytical Gradient
#' @name distrib_gradient.GenGamma1Distrib
#' @description
#' With \eqn{w = (y/a)^{p}}, \eqn{L = \log(y/a)} and \eqn{k = d/p},
#' \deqn{\dfrac{\partial\ell}{\partial a} = \dfrac{pw - d}{a}, \qquad
#'       \dfrac{\partial\ell}{\partial d} = L - \dfrac{\psi(k)}{p}, \qquad
#'       \dfrac{\partial\ell}{\partial p} = \dfrac{1}{p}
#'         + \dfrac{d\,\psi(k)}{p^{2}} - wL}
#' @param distrib A \code{GenGamma1Distrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing \code{a}, \code{d} and \code{p}.
#' @param scale Either \code{"parameter"} or \code{"link"}; handled by the generic.
#' @param ... Unused.
#' @return A named list with the \code{a}, \code{d} and \code{p} components.
#' @seealso \code{\link{gengamma1_distrib}}
S7::method(distrib_gradient, GenGamma1Distrib) <- function(distrib, y, theta,
                                                           scale = c("parameter", "link"), ..., threads = 1L) {
  gengamma_gradient_cpp(y, theta[[1]], theta[[2]], theta[[3]], threads)
}

#' @title Generalized Gamma Analytical Observed Hessian
#' @name distrib_hessian.GenGamma1Distrib
#' @description
#' The second derivatives of the same expressions. The mixed
#' \eqn{a}--\eqn{d} component is \eqn{-1/a}, free of the data, the scale and
#' the first shape entering the log-density through \eqn{-d\log a} alone.
#' @param distrib A \code{GenGamma1Distrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing \code{a}, \code{d} and \code{p}.
#' @param scale Either \code{"parameter"} or \code{"link"}; handled by the generic.
#' @param ... Unused.
#' @return A named list of second-derivative components.
#' @seealso \code{\link{gengamma1_distrib}}
S7::method(distrib_hessian, GenGamma1Distrib) <- function(distrib, y, theta,
                                                          scale = c("parameter", "link"), ..., threads = 1L) {
  gengamma_hessian_cpp(y, theta[[1]], theta[[2]], theta[[3]], threads)
}

#' @title Generalized Gamma Analytical Expected Hessian
#' @name distrib_expected_hessian.GenGamma1Distrib
#' @description
#' Closed form. Every expectation the observed Hessian needs is a moment of
#' \eqn{u = (Y/a)^{p}}, which is Gamma with shape \eqn{k = d/p} and unit rate:
#' \eqn{\mathbb{E}[u] = k}, \eqn{\mathbb{E}[u\log u] = k\psi(k+1)} and
#' \eqn{\mathbb{E}[u(\log u)^{2}] = k\{\psi(k+1)^{2} + \psi'(k+1)\}}, so
#' \code{approx} is ignored.
#' @param distrib A \code{GenGamma1Distrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing \code{a}, \code{d} and \code{p}.
#' @param scale Either \code{"parameter"} or \code{"link"}; handled by the generic.
#' @param approx Ignored; the expectation is closed form.
#' @param nsim Ignored.
#' @param ... Unused.
#' @return A named list of expected second-derivative components.
#' @seealso \code{\link{gengamma1_distrib}}
S7::method(distrib_expected_hessian, GenGamma1Distrib) <- function(distrib, y, theta,
                                                                   scale = c("parameter", "link"),
                                                                   approx = c("bartlett", "integrate", "mc", "opg"),
                                                                   nsim = 10000, ..., threads = 1L) {
  gengamma_expected_hessian_cpp(y, theta[[1]], theta[[2]], theta[[3]], threads)
}

# --- CONSTRUCTOR WRAPPER ---

#' Generalized Gamma Distribution Object
#'
#' @description
#' Creates a distribution object for the generalized gamma distribution in
#' Stacy's form, with a scale \eqn{a} and two shapes \eqn{d} and \eqn{p}.
#'
#' @param link_a A link function object for \eqn{a}. Defaults to
#'   \code{\link[linkfunctions7]{log_link}}.
#' @param link_d A link function object for \eqn{d}. Defaults to
#'   \code{\link[linkfunctions7]{log_link}}.
#' @param link_p A link function object for \eqn{p}. Defaults to
#'   \code{\link[linkfunctions7]{log_link}}.
#'
#' @details
#' The flexible family for a positive response, and the one that makes a
#' choice between the gamma, the Weibull and the lognormal something to
#' estimate rather than to assume.
#'
#' \strong{Density:}
#' \deqn{f(y) = \dfrac{p}{a^{d}\,\Gamma(d/p)}\,y^{d-1}e^{-(y/a)^{p}}}
#'
#' \strong{What it nests}, which Stacy's parametrization is chosen to make
#' visible:
#' \itemize{
#'   \item \eqn{p = 1} is the \link[=gamma2_distrib]{gamma} with shape \eqn{d}
#'     and scale \eqn{a};
#'   \item \eqn{d = p} is the \link[=weibull1_distrib]{Weibull} with shape
#'     \eqn{p} and scale \eqn{a};
#'   \item \eqn{d = p = 1} is the \link[=exponential_distrib]{exponential};
#'   \item \eqn{p \to 0} with \eqn{d/p} held large approaches the
#'     \link[=lognormal1_distrib]{lognormal}.
#' }
#' The first three are exact and testable; the fourth is a limit and is not
#' reached at any admissible value.
#'
#' \strong{Score and information.} Writing \eqn{w = (y/a)^{p}},
#' \eqn{L = \log(y/a)} and \eqn{k = d/p},
#' \deqn{\dfrac{\partial\ell}{\partial a} = \dfrac{pw-d}{a}, \qquad
#'       \dfrac{\partial\ell}{\partial d} = L - \dfrac{\psi(k)}{p}, \qquad
#'       \dfrac{\partial\ell}{\partial p} = \dfrac{1}{p}
#'         + \dfrac{d\psi(k)}{p^{2}} - wL}
#' and the expected information is \strong{closed form}: \eqn{u = w} is Gamma
#' with shape \eqn{k} and unit rate, so every expectation the Hessian needs is
#' one of \eqn{\mathbb{E}[u] = k}, \eqn{\mathbb{E}[u\log u] = k\psi(k+1)} and
#' \eqn{\mathbb{E}[u(\log u)^{2}] = k\{\psi(k+1)^2 + \psi'(k+1)\}}. That is the
#' same device the Weibull and the Gumbel use, where the corresponding
#' variable is standard exponential.
#'
#' \strong{Moments:}
#' \eqn{\mathbb{E}[Y^{r}] = a^{r}\Gamma\{(d+r)/p\}/\Gamma(d/p)}, finite for
#' every \eqn{r > -d}.
#'
#' \strong{Parameter domains:} all three are positive.
#'
#' The three parameters are \strong{weakly identified together} on small
#' samples: \eqn{d} and \eqn{p} enter the density largely through their ratio,
#' and the profile likelihood in that direction is flat. A fit of all three
#' wants several hundred observations, and holding one with
#' \code{\link{fixed}} is often the better model.
#'
#' @return An S7 object of class \code{GenGamma1Distrib}.
#'
#' @references
#' Stacy, E. W. (1962). A generalization of the gamma distribution.
#' \emph{Annals of Mathematical Statistics} 33, 1187-1192.
#'
#' @seealso \code{\link{gamma2_distrib}}, \code{\link{weibull1_distrib}},
#'   \code{\link{lognormal1_distrib}}, \code{\link{exponential_distrib}}
#'
#' @importFrom linkfunctions7 log_link
#' @importFrom stats pgamma qgamma rgamma
#' @examples
#' d <- gengamma1_distrib()
#' d@params
#'
#' theta <- list(a = 2, d = 3, p = 1.5)
#' distrib_pdf(d, c(0.5, 2, 5), theta)
#'
#' # p = 1 is the gamma with shape d and scale a
#' max(abs(distrib_pdf(d, c(0.5, 2, 5), list(a = 2, d = 3, p = 1)) -
#'         dgamma(c(0.5, 2, 5), shape = 3, scale = 2)))
#'
#' # d = p is the Weibull with shape p and scale a
#' max(abs(distrib_pdf(d, c(0.5, 2, 5), list(a = 2, d = 1.5, p = 1.5)) -
#'         dweibull(c(0.5, 2, 5), shape = 1.5, scale = 2)))
#'
#' @export
gengamma1_distrib <- function(link_a = log_link(), link_d = log_link(),
                             link_p = log_link()) {
  GenGamma1Distrib(
    distrib_name = "gengamma1", dimension = "univariate",
    bounds = c(0, Inf),
    params = c("a", "d", "p"),
    params_interpretation = c(a = "scale", d = "shape", p = "power"),
    n_params = 3,
    params_bounds = list(a = c(0, Inf), d = c(0, Inf), p = c(0, Inf)),
    link_params = list(a = link_a, d = link_d, p = link_p)
  )
}
