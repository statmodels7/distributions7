#' @include distrib.R generics.R
NULL

#' @title S7 Class for the Beta-Binomial in Its Shapes
#' @name BetaBinom2Distrib
#'
#' @description A subclass of `discrete_distrib` for the beta-binomial in
#'   its canonical parametrization, the two beta shapes.
#' @inheritParams distrib
#' @param size The number of trials, a constant of the distribution.
#' @return An object of class `BetaBinom2Distrib`.
#' @seealso [betabinom2_distrib()], [betabinom1_distrib()]
#'
#' @section Methods:
#' Methods implemented for this class:
#'   [`distrib_deriv3()`][distrib_deriv3.BetaBinom2Distrib],
#'   [`distrib_deriv4()`][distrib_deriv4.BetaBinom2Distrib],
#'   [`distrib_expected_hessian()`][distrib_expected_hessian.BetaBinom2Distrib],
#'   [`distrib_gradient()`][distrib_gradient.BetaBinom2Distrib],
#'   [`distrib_hessian()`][distrib_hessian.BetaBinom2Distrib],
#'   [`distrib_pdf()`][distrib_pdf.BetaBinom2Distrib],
#'   [`distrib_rng()`][distrib_rng.BetaBinom2Distrib]
#'
#' Everything else is inherited from [discrete_distrib()].
BetaBinom2Distrib <- S7::new_class("BetaBinom2Distrib",
  parent = discrete_distrib,
  properties = list(size = S7::class_numeric)
)

#' Higher Derivatives of the Beta-Binomial in Its Shapes
#'
#' @description
#' The log-mass function is a sum of log-gamma terms, so its derivative of
#' order \eqn{k} is the same sum with \eqn{\psi^{(k-1)}} in place of
#' \eqn{\log\Gamma}. All four orders follow from one routine.
#'
#' @details
#' With \eqn{a} and \eqn{b} the shapes and \eqn{n} the size, the log-mass is
#' \eqn{\log\Gamma(y+a) + \log\Gamma(n-y+b) - \log\Gamma(n+a+b)
#'      - \log\Gamma(a) - \log\Gamma(b) + \log\Gamma(a+b)} up to a constant.
#' A derivative in \eqn{a} alone differentiates the first, third, fourth and
#' sixth terms; one in \eqn{b} alone the second, third, fifth and sixth; a
#' mixed one only the two terms carrying \eqn{a+b}.
#'
#' @param y A numeric vector of observations.
#' @param a,b The two shapes.
#' @param n The size.
#' @param k The polygamma order, `order - 1`.
#' @param i The number of \eqn{a} indices in the component.
#' @param j The number of \eqn{b} indices.
#'
#' @return A numeric vector.
#'
#' @seealso [betabinom2_distrib()]
#'
#' @keywords internal
betabinom2_component <- function(y, a, b, n, k, i, j) {
  mixed <- -psigamma(n + a + b, deriv = k) + psigamma(a + b, deriv = k)
  if (i > 0 && j > 0) return(rep(mixed, length.out = length(y)))
  if (j == 0) {
    return(psigamma(y + a, deriv = k) - psigamma(a, deriv = k) + mixed)
  }
  psigamma(n - y + b, deriv = k) - psigamma(b, deriv = k) + mixed
}

#' Every Component of a Beta-Binomial Derivative
#'
#' @description
#' Assembles the components of a derivative of the given order, named as
#' [deriv_names()] names them.
#'
#' @param y A numeric vector of observations.
#' @param a,b The two shapes.
#' @param n The size.
#' @param order The derivative order.
#' @param params The parameter names.
#'
#' @return A named list of component vectors.
#'
#' @seealso [betabinom2_distrib()]
#'
#' @keywords internal
betabinom2_derivs <- function(y, a, b, n, order, params) {
  idx <- deriv_indices(params, order)
  nm <- deriv_names(params, order)
  stats::setNames(lapply(seq_along(nm), function(m) {
    id <- idx[[m]]
    betabinom2_component(y, a, b, n, order - 1L, sum(id == 1L), sum(id == 2L))
  }), nm)
}

# --- S7 METHODS IMPLEMENTATION ---

#' @title Beta-Binomial Mass Function in Its Shapes
#' @name distrib_pdf.BetaBinom2Distrib
#' @description
#' \deqn{P(Y = y) = \binom{n}{y}\dfrac{B(y+\alpha, n-y+\beta)}{B(\alpha,\beta)}}
#' @param distrib A `BetaBinom2Distrib` object.
#' @param y A numeric vector of counts.
#' @param theta A list with `alpha` and `beta`.
#' @param log Logical; if `TRUE`, returns the log-probability.
#' @return A numeric vector.
#' @seealso [betabinom2_distrib()]
S7::method(distrib_pdf, BetaBinom2Distrib) <- function(distrib, y, theta, log = FALSE, ...) {
  a <- theta[[1]]
  b <- theta[[2]]
  n <- distrib@size
  # The support is tested BEFORE lchoose(), which warns on a non-integer count
  # rather than returning nothing: a mass function evaluated off its support is
  # zero, not a numerical complaint.
  ok <- y >= 0 & y <= n & y == round(y)
  out <- rep(-Inf, length(y))
  if (any(ok)) {
    out[ok] <- betabinom_log_mass(y[ok], a, b, n)
  }
  if (log) out else exp(out)
}

#' Log-Mass of the Beta-Binomial
#'
#' @description
#' The log-mass \eqn{\log\binom{n}{y} + \log B(y+\alpha, n-y+\beta) -
#' \log B(\alpha, \beta)} by whichever of two routes is accurate at the shapes
#' given.
#'
#' @details
#' The two beta functions are of magnitude \eqn{(\alpha+\beta)\log(\alpha+\beta)}
#' and their difference is of order one, so the ordinary route carries an
#' absolute error of \eqn{\varepsilon} times that magnitude and is used only
#' while this stays below `1e-8`. Beyond it the shifts are integers, so
#' each log-gamma difference is an exact sum of logarithms,
#' \deqn{\log\Gamma(\alpha+y) - \log\Gamma(\alpha) =
#'       \sum_{j=0}^{y-1}\log(\alpha+j),}
#' and the mass follows from three such sums without forming any quantity
#' larger than \eqn{n\log(\alpha+\beta)}. The sums also give the binomial limit
#' correctly as the shapes grow at a fixed ratio.
#'
#' @param y A numeric vector of counts, already known to lie on the support.
#' @param a,b The two shapes.
#' @param n The size.
#' @return A numeric vector of log-probabilities.
#' @keywords internal
betabinom_log_mass <- function(y, a, b, n) {
  s <- a + b
  # the two large terms are grouped so that their difference is taken before
  # the term of order one is added to it
  if (all(is.finite(s)) &&
      max(lgamma(s + n)) * .Machine$double.eps < 1e-8) {
    return(lchoose(n, y) + (lbeta(y + a, n - y + b) - lbeta(a, b)))
  }
  m <- length(y)
  a <- rep_len(a, m)
  b <- rep_len(b, m)
  s1 <- numeric(m)
  s2 <- numeric(m)
  s3 <- numeric(m)
  for (j in seq_len(n) - 1) {
    s1 <- s1 + (j < y) * log(a + j)
    s2 <- s2 + (j < n - y) * log(b + j)
    s3 <- s3 + log(a + b + j)
  }
  lchoose(n, y) + s1 + s2 - s3
}

#' @title Beta-Binomial Random Generation in Its Shapes
#' @name distrib_rng.BetaBinom2Distrib
#' @description A beta draw for the probability, then a binomial draw.
#' @param distrib A `BetaBinom2Distrib` object.
#' @param n The number of draws.
#' @param theta A list with `alpha` and `beta`.
#' @return A numeric vector of counts.
#' @seealso [betabinom2_distrib()]
S7::method(distrib_rng, BetaBinom2Distrib) <- function(distrib, n, theta) {
  p <- stats::rbeta(n, shape1 = theta[[1]], shape2 = theta[[2]])
  stats::rbinom(n, size = distrib@size, prob = p)
}

#' @title Beta-Binomial Analytical Gradient in Its Shapes
#' @name distrib_gradient.BetaBinom2Distrib
#' @description
#' \deqn{\dfrac{\partial\ell}{\partial\alpha}
#'         = \psi(y+\alpha) - \psi(n+\alpha+\beta)
#'           - \psi(\alpha) + \psi(\alpha+\beta)}
#' and the same with \eqn{n-y} and \eqn{\beta}.
#' @param distrib A `BetaBinom2Distrib` object.
#' @param y A numeric vector of counts.
#' @param theta A list with `alpha` and `beta`.
#' @param scale Either `"parameter"` or `"link"`; handled by the generic.
#' @param ... Unused.
#' @return A named list of first derivatives.
#' @seealso [betabinom2_distrib()]
S7::method(distrib_gradient, BetaBinom2Distrib) <- function(distrib, y, theta,
                                                             scale = c("parameter", "link"), ...) {
  a <- theta[[1]]
  b <- theta[[2]]
  n <- distrib@size
  ds <- digamma(a + b) - digamma(n + a + b)
  list(alpha = digamma(y + a) - digamma(a) + ds,
       beta = digamma(n - y + b) - digamma(b) + ds)
}

#' @title Beta-Binomial Analytical Observed Hessian in Its Shapes
#' @name distrib_hessian.BetaBinom2Distrib
#' @description The same sums with the trigamma function.
#' @param distrib A `BetaBinom2Distrib` object.
#' @param y A numeric vector of counts.
#' @param theta A list with `alpha` and `beta`.
#' @param scale Either `"parameter"` or `"link"`; handled by the generic.
#' @param ... Unused.
#' @return A named list of second derivatives.
#' @seealso [betabinom2_distrib()]
S7::method(distrib_hessian, BetaBinom2Distrib) <- function(distrib, y, theta,
                                                            scale = c("parameter", "link"), ...) {
  d <- betabinom2_derivs(y, theta[[1]], theta[[2]], distrib@size, 2L,
                         distrib@params)
  d[hess_names(distrib@params)]
}

#' @title Beta-Binomial Analytical Expected Hessian in Its Shapes
#' @name distrib_expected_hessian.BetaBinom2Distrib
#' @description
#' An **exact finite sum** over the support, the family being discrete on
#' \eqn{\{0, \dots, n\}}: the expectation is a weighted sum of at most
#' \eqn{n+1} terms rather than a quadrature or a sample.
#' @param distrib A `BetaBinom2Distrib` object.
#' @param y A numeric vector of counts.
#' @param theta A list with `alpha` and `beta`.
#' @param scale Either `"parameter"` or `"link"`; handled by the generic.
#' @param approx Ignored; the expectation is an exact sum.
#' @param nsim Ignored.
#' @param ... Unused.
#' @return A named list of expected second derivatives.
#' @seealso [betabinom2_distrib()]
S7::method(distrib_expected_hessian, BetaBinom2Distrib) <- function(distrib, y, theta,
                                                                     scale = c("parameter", "link"),
                                                                     approx = c("bartlett", "integrate", "mc", "opg"),
                                                                     nsim = 10000, ...) {
  betabinom2_expected(distrib, y, theta, 2L)
}

#' Expected Derivatives of the Beta-Binomial by Exact Summation
#'
#' @description
#' Averages a derivative over the support \eqn{\{0, \dots, n\}} weighted by the
#' mass function, which is exact because the support is finite.
#'
#' @param distrib A [BetaBinom2Distrib()] object.
#' @param y A numeric vector, used only for its length.
#' @param theta A list with `alpha` and `beta`.
#' @param order The derivative order.
#'
#' @return A named list of component vectors.
#'
#' @seealso [betabinom2_distrib()]
#'
#' @keywords internal
betabinom2_expected <- function(distrib, y, theta, order) {
  a <- theta[[1]]
  b <- theta[[2]]
  n <- distrib@size
  supp <- 0:n
  w <- distrib_pdf(distrib, supp, theta)
  d <- betabinom2_derivs(supp, a, b, n, order, distrib@params)
  nm <- if (order == 2L) hess_names(distrib@params) else names(d)
  stats::setNames(lapply(nm, function(k) {
    rep(sum(w * d[[k]]), length.out = length(y))
  }), nm)
}

#' @title Beta-Binomial Third-Order Derivatives in Its Shapes
#' @name distrib_deriv3.BetaBinom2Distrib
#' @description Closed form, with the expectation an exact sum over the
#'   support.
#' @param distrib A `BetaBinom2Distrib` object.
#' @param y A numeric vector of counts.
#' @param theta A list with `alpha` and `beta`.
#' @param expected Logical; if `TRUE`, returns the expected derivatives.
#' @param scale Either `"parameter"` or `"link"`; handled by the generic.
#' @param approx Ignored.
#' @param nsim Ignored.
#' @param ... Unused.
#' @return A named list of third-derivative components.
#' @seealso [betabinom2_distrib()]
S7::method(distrib_deriv3, BetaBinom2Distrib) <- function(distrib, y, theta, expected = FALSE,
                                                           scale = c("parameter", "link"),
                                                           approx = c("integrate", "bartlett", "mc", "opg"),
                                                           nsim = 10000, ...) {
  if (expected) {
    betabinom2_expected(distrib, y, theta, 3L)
  } else {
    betabinom2_derivs(y, theta[[1]], theta[[2]], distrib@size, 3L, distrib@params)
  }
}

#' @title Beta-Binomial Fourth-Order Derivatives in Its Shapes
#' @name distrib_deriv4.BetaBinom2Distrib
#' @description Closed form, with the expectation an exact sum over the
#'   support.
#' @param distrib A `BetaBinom2Distrib` object.
#' @param y A numeric vector of counts.
#' @param theta A list with `alpha` and `beta`.
#' @param expected Logical; if `TRUE`, returns the expected derivatives.
#' @param scale Either `"parameter"` or `"link"`; handled by the generic.
#' @param approx Ignored.
#' @param nsim Ignored.
#' @param ... Unused.
#' @return A named list of fourth-derivative components.
#' @seealso [betabinom2_distrib()]
S7::method(distrib_deriv4, BetaBinom2Distrib) <- function(distrib, y, theta, expected = FALSE,
                                                           scale = c("parameter", "link"),
                                                           approx = c("integrate", "bartlett", "mc", "opg"),
                                                           nsim = 10000, ...) {
  if (expected) {
    betabinom2_expected(distrib, y, theta, 4L)
  } else {
    betabinom2_derivs(y, theta[[1]], theta[[2]], distrib@size, 4L, distrib@params)
  }
}


#' Beta-Binomial Distribution in Its Shapes
#'
#' @description
#' Creates a beta-binomial distribution object in its canonical
#' parametrization, the two beta shapes.
#'
#' @details
#' The same law as [betabinom1_distrib()], which carries a mean
#' proportion and a dispersion: \eqn{\alpha = \mu/\sigma} and
#' \eqn{\beta = (1-\mu)/\sigma}.
#'
#' Every derivative is a sum of polygamma functions, since the log-mass is a
#' sum of log-gamma terms and differentiating it \eqn{k} times replaces each
#' by \eqn{\psi^{(k-1)}}. The expectations are **exact finite sums** over
#' \eqn{\{0, \dots, n\}} rather than quadratures.
#'
#' @section The distribution:
#' \deqn{P(Y=y) = \binom{n}{y}\frac{B(y+\alpha,\; n-y+\beta)}{B(\alpha, \beta)}}
#' on \eqn{y \in \{0, \dots, n\}}.
#'
#' \deqn{\mathbb{E}[Y] = \frac{n\alpha}{\alpha+\beta}, \qquad \operatorname{Var}(Y) = \frac{n\alpha\beta\,(\alpha+\beta+n)}{(\alpha+\beta)^{2}(\alpha+\beta+1)}}
#'
#' @param size The number of trials, a constant of the distribution.
#' @param link_alpha Link function for \eqn{\alpha}. Defaults to the log.
#' @param link_beta Link function for \eqn{\beta}. Defaults to the log.
#'
#' @return An S7 object of class [BetaBinom2Distrib()].
#'
#' @seealso [betabinom1_distrib()]
#'
#' @examples
#' d <- betabinom2_distrib(size = 10)
#' theta <- list(alpha = 2, beta = 3)
#' distrib_pdf(d, 0:10, theta)
#' c(mean = mean(d, theta), variance = variance(d, theta))
#'
#' @export
betabinom2_distrib <- function(size, link_alpha = log_link(),
                               link_beta = log_link()) {
  if (length(size) != 1L || !is.finite(size) || size < 1 || size != round(size)) {
    stop("'size' must be a single positive integer.", call. = FALSE)
  }
  BetaBinom2Distrib(
    distrib_name = paste0("betabinom2 [size=", size, "]"),
    dimension = "univariate",
    bounds = c(0, size),
    params = c("alpha", "beta"),
    params_interpretation = c(alpha = "shape", beta = "shape"),
    n_params = 2,
    params_bounds = list(alpha = c(0, Inf), beta = c(0, Inf)),
    link_params = list(alpha = link_alpha, beta = link_beta),
    size = size
  )
}
