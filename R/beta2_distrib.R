#' @include distrib.R generics.R
NULL

#' @title S7 Class for the Beta Distribution in Its Shapes
#' @name Beta2Distrib
#'
#' @description A subclass of `continuous_distrib` for the beta in its
#'   canonical parametrization, the two shapes.
#' @inheritParams distrib
#' @return An object of class `Beta2Distrib`.
#' @seealso [beta2_distrib()], [beta1_distrib()]
#'
#' @section Methods:
#' Methods implemented for this class:
#'   [`distrib_cdf()`][distrib_cdf.Beta2Distrib],
#'   [`distrib_deriv3()`][distrib_deriv3.Beta2Distrib],
#'   [`distrib_deriv4()`][distrib_deriv4.Beta2Distrib],
#'   [`distrib_expected_hessian()`][distrib_expected_hessian.Beta2Distrib],
#'   [`distrib_grad_y()`][distrib_grad_y.Beta2Distrib],
#'   [`distrib_hess_y()`][distrib_hess_y.Beta2Distrib],
#'   [`distrib_gradient()`][distrib_gradient.Beta2Distrib],
#'   [`distrib_hessian()`][distrib_hessian.Beta2Distrib],
#'   [`distrib_pdf()`][distrib_pdf.Beta2Distrib],
#'   [`distrib_quantile()`][distrib_quantile.Beta2Distrib],
#'   [`distrib_rng()`][distrib_rng.Beta2Distrib]
#'
#' Everything else is inherited from [continuous_distrib()].
Beta2Distrib <- S7::new_class("Beta2Distrib", parent = continuous_distrib)

# --- S7 METHODS IMPLEMENTATION ---

#' @title Beta Density in Its Shapes
#' @name distrib_pdf.Beta2Distrib
#' @description
#' \deqn{f(y) = \dfrac{y^{\alpha-1}(1-y)^{\beta-1}}{B(\alpha, \beta)}}
#' @param distrib A `Beta2Distrib` object.
#' @param y A numeric vector of observations in \eqn{(0, 1)}.
#' @param theta A list with `alpha` and `beta`.
#' @param log Logical; if `TRUE`, returns the log-density.
#' @return A numeric vector.
#' @seealso [beta2_distrib()]
S7::method(distrib_pdf, Beta2Distrib) <- function(distrib, y, theta, log = FALSE, ...) {
  stats::dbeta(y, shape1 = theta[[1]], shape2 = theta[[2]], log = log)
}

#' @title Beta Distribution Function in Its Shapes
#' @name distrib_cdf.Beta2Distrib
#' @description The incomplete beta function.
#' @param distrib A `Beta2Distrib` object.
#' @param q A numeric vector of quantiles.
#' @param theta A list with `alpha` and `beta`.
#' @param lower.tail Logical; if `TRUE`, probabilities are \eqn{P(Y \le q)}.
#' @param log.p Logical; if `TRUE`, returns log-probabilities.
#' @param ... Unused.
#' @return A numeric vector.
#' @seealso [beta2_distrib()]
S7::method(distrib_cdf, Beta2Distrib) <- function(distrib, q, theta,
                                                   lower.tail = TRUE,
                                                   log.p = FALSE, ...) {
  stats::pbeta(q, shape1 = theta[[1]], shape2 = theta[[2]],
               lower.tail = lower.tail, log.p = log.p)
}

#' @title Beta Quantile Function in Its Shapes
#' @name distrib_quantile.Beta2Distrib
#' @description The inverse of the incomplete beta function.
#' @param distrib A `Beta2Distrib` object.
#' @param p A numeric vector of probabilities.
#' @param theta A list with `alpha` and `beta`.
#' @param lower.tail Logical; if `TRUE`, probabilities are \eqn{P(Y \le q)}.
#' @param log.p Logical; if `TRUE`, `p` is a log-probability.
#' @param ... Unused.
#' @return A numeric vector.
#' @seealso [beta2_distrib()]
S7::method(distrib_quantile, Beta2Distrib) <- function(distrib, p, theta,
                                                        lower.tail = TRUE,
                                                        log.p = FALSE, ...) {
  stats::qbeta(p, shape1 = theta[[1]], shape2 = theta[[2]],
               lower.tail = lower.tail, log.p = log.p)
}

#' @title Beta Random Generation in Its Shapes
#' @name distrib_rng.Beta2Distrib
#' @description Delegates to [stats::rbeta()].
#' @param distrib A `Beta2Distrib` object.
#' @param n The number of draws.
#' @param theta A list with `alpha` and `beta`.
#' @return A numeric vector.
#' @seealso [beta2_distrib()]
S7::method(distrib_rng, Beta2Distrib) <- function(distrib, n, theta) {
  stats::rbeta(n, shape1 = theta[[1]], shape2 = theta[[2]])
}

#' @title Beta Analytical Gradient in Its Shapes
#' @name distrib_gradient.Beta2Distrib
#' @description
#' \deqn{\dfrac{\partial\ell}{\partial\alpha}
#'         = \log y - \psi(\alpha) + \psi(\alpha+\beta), \qquad
#'       \dfrac{\partial\ell}{\partial\beta}
#'         = \log(1-y) - \psi(\beta) + \psi(\alpha+\beta)}
#' The data enter only through \eqn{\log y} and \eqn{\log(1-y)}, which is what
#' makes every higher derivative free of them.
#' @param distrib A `Beta2Distrib` object.
#' @param y A numeric vector of observations.
#' @param theta A list with `alpha` and `beta`.
#' @param scale Either `"parameter"` or `"link"`; handled by the generic.
#' @param ... Unused.
#' @return A named list of first derivatives.
#' @seealso [beta2_distrib()]
S7::method(distrib_gradient, Beta2Distrib) <- function(distrib, y, theta,
                                                        scale = c("parameter", "link"), ...) {
  a <- theta[[1]]
  b <- theta[[2]]
  ds <- digamma(a + b)
  list(alpha = log(y) - digamma(a) + ds,
       beta = log1p(-y) - digamma(b) + ds)
}

#' Higher Derivatives of the Beta in Its Shapes
#'
#' @description
#' Orders two to four, which are free of the data: each is a difference of
#' polygamma functions at \eqn{\alpha}, \eqn{\beta} and \eqn{\alpha+\beta}.
#'
#' @details
#' Being free of the data, the observed and the expected derivatives coincide
#' at every order beyond the first, which is why this family needs no
#' expectation anywhere.
#'
#' @param theta A list with `alpha` and `beta`.
#' @param n The number of observations to recycle to.
#' @param order The derivative order, 2, 3 or 4.
#'
#' @return A named list of component vectors.
#'
#' @seealso [beta2_distrib()]
#'
#' @keywords internal
beta2_higher <- function(theta, n, order) {
  a <- theta[[1]]
  b <- theta[[2]]
  k <- order - 1L
  pa <- psigamma(a, deriv = k)
  pb <- psigamma(b, deriv = k)
  ps <- psigamma(a + b, deriv = k)
  rep_n <- function(v) rep(v, length.out = n)
  if (order == 2L) {
    return(list(alpha_alpha = rep_n(ps - pa),
                alpha_beta = rep_n(ps),
                beta_beta = rep_n(ps - pb)))
  }
  if (order == 3L) {
    return(list(alpha_alpha_alpha = rep_n(ps - pa),
                alpha_alpha_beta = rep_n(ps),
                alpha_beta_beta = rep_n(ps),
                beta_beta_beta = rep_n(ps - pb)))
  }
  list(alpha_alpha_alpha_alpha = rep_n(ps - pa),
       alpha_alpha_alpha_beta = rep_n(ps),
       alpha_alpha_beta_beta = rep_n(ps),
       alpha_beta_beta_beta = rep_n(ps),
       beta_beta_beta_beta = rep_n(ps - pb))
}

#' @title Beta Analytical Observed Hessian in Its Shapes
#' @name distrib_hessian.Beta2Distrib
#' @description
#' \deqn{\ell^{(\alpha\alpha)} = \psi'(\alpha+\beta) - \psi'(\alpha), \qquad
#'       \ell^{(\alpha\beta)} = \psi'(\alpha+\beta), \qquad
#'       \ell^{(\beta\beta)} = \psi'(\alpha+\beta) - \psi'(\beta)}
#' free of the data.
#' @param distrib A `Beta2Distrib` object.
#' @param y A numeric vector of observations.
#' @param theta A list with `alpha` and `beta`.
#' @param scale Either `"parameter"` or `"link"`; handled by the generic.
#' @param ... Unused.
#' @return A named list of second derivatives.
#' @seealso [beta2_distrib()]
S7::method(distrib_hessian, Beta2Distrib) <- function(distrib, y, theta,
                                                       scale = c("parameter", "link"), ...) {
  h <- beta2_higher(theta, length(y), 2L)
  h[hess_names(distrib@params)]
}

#' @title Beta Analytical Expected Hessian in Its Shapes
#' @name distrib_expected_hessian.Beta2Distrib
#' @description
#' Equal to the observed Hessian: it is free of the data, so there is nothing
#' to average. Fisher scoring and Newton's method therefore take the same step
#' on the parameter scale here.
#' @param distrib A `Beta2Distrib` object.
#' @param y A numeric vector of observations.
#' @param theta A list with `alpha` and `beta`.
#' @param scale Either `"parameter"` or `"link"`; handled by the generic.
#' @param approx Ignored; the expectation is exact.
#' @param nsim Ignored.
#' @param ... Unused.
#' @return A named list of expected second derivatives.
#' @seealso [beta2_distrib()]
S7::method(distrib_expected_hessian, Beta2Distrib) <- function(distrib, y, theta,
                                                                scale = c("parameter", "link"),
                                                                approx = c("bartlett", "integrate", "mc", "opg"),
                                                                nsim = 10000, ...) {
  h <- beta2_higher(theta, length(y), 2L)
  h[hess_names(distrib@params)]
}

#' @title Beta Third-Order Derivatives in Its Shapes
#' @name distrib_deriv3.Beta2Distrib
#' @description Closed form and free of the data, so `expected` changes
#'   nothing.
#' @param distrib A `Beta2Distrib` object.
#' @param y A numeric vector of observations.
#' @param theta A list with `alpha` and `beta`.
#' @param expected Logical; makes no difference here.
#' @param scale Either `"parameter"` or `"link"`; handled by the generic.
#' @param approx Ignored.
#' @param nsim Ignored.
#' @param ... Unused.
#' @return A named list of third-derivative components.
#' @seealso [beta2_distrib()]
S7::method(distrib_deriv3, Beta2Distrib) <- function(distrib, y, theta, expected = FALSE,
                                                      scale = c("parameter", "link"),
                                                      approx = c("integrate", "bartlett", "mc", "opg"),
                                                      nsim = 10000, ...) {
  beta2_higher(theta, length(y), 3L)
}

#' @title Beta Fourth-Order Derivatives in Its Shapes
#' @name distrib_deriv4.Beta2Distrib
#' @description Closed form and free of the data, so `expected` changes
#'   nothing.
#' @param distrib A `Beta2Distrib` object.
#' @param y A numeric vector of observations.
#' @param theta A list with `alpha` and `beta`.
#' @param expected Logical; makes no difference here.
#' @param scale Either `"parameter"` or `"link"`; handled by the generic.
#' @param approx Ignored.
#' @param nsim Ignored.
#' @param ... Unused.
#' @return A named list of fourth-derivative components.
#' @seealso [beta2_distrib()]
S7::method(distrib_deriv4, Beta2Distrib) <- function(distrib, y, theta, expected = FALSE,
                                                      scale = c("parameter", "link"),
                                                      approx = c("integrate", "bartlett", "mc", "opg"),
                                                      nsim = 10000, ...) {
  beta2_higher(theta, length(y), 4L)
}

#' @title Beta Response Derivatives in Its Shapes
#' @name distrib_grad_y.Beta2Distrib
#' @description
#' \eqn{\partial\ell/\partial y = (\alpha-1)/y - (\beta-1)/(1-y)}.
#' @param distrib A `Beta2Distrib` object.
#' @param y A numeric vector of observations.
#' @param theta A list with `alpha` and `beta`.
#' @param ... Unused.
#' @return A numeric vector.
#' @seealso [beta2_distrib()]
S7::method(distrib_grad_y, Beta2Distrib) <- function(distrib, y, theta, ...) {
  (theta[[1]] - 1) / y - (theta[[2]] - 1) / (1 - y)
}

#' @title Beta Second Response Derivative in Its Shapes
#' @name distrib_hess_y.Beta2Distrib
#' @description
#' \eqn{\partial^2\ell/\partial y^2 = -(\alpha-1)/y^2 - (\beta-1)/(1-y)^2}.
#' @param distrib A `Beta2Distrib` object.
#' @param y A numeric vector of observations.
#' @param theta A list with `alpha` and `beta`.
#' @param ... Unused.
#' @return A numeric vector.
#' @seealso [beta2_distrib()]
S7::method(distrib_hess_y, Beta2Distrib) <- function(distrib, y, theta, ...) {
  -(theta[[1]] - 1) / y^2 - (theta[[2]] - 1) / (1 - y)^2
}


#' Beta Distribution in Its Shapes
#'
#' @description
#' Creates a beta distribution object in its canonical parametrization, the two
#' shapes \eqn{\alpha} and \eqn{\beta}.
#'
#' @details
#' The same law as [beta1_distrib()], which carries the mean and a
#' precision: \eqn{\alpha = \mu\varphi} and \eqn{\beta = (1-\mu)\varphi}. The
#' mean parametrization is the one a regression wants; this one is the one the
#' family is usually written in and the one a conjugate analysis produces, the
#' beta being conjugate for a binomial probability.
#'
#' The data enter the log-density only through \eqn{\log y} and
#' \eqn{\log(1-y)}, both of which are linear in the parameters. Every
#' derivative beyond the first is therefore free of the data, so the observed
#' and the expected ones coincide at orders two, three and four, and Fisher
#' scoring and Newton's method take the same step on the parameter scale.
#'
#' @section The distribution:
#' \deqn{f(y) = \frac{y^{\alpha-1}(1-y)^{\beta-1}}{B(\alpha, \beta)}}
#' on \eqn{y \in (0, 1)}.
#'
#' \deqn{\mathbb{E}[Y] = \frac{\alpha}{\alpha+\beta}, \qquad \operatorname{Var}(Y) = \frac{\alpha\beta}{(\alpha+\beta)^{2}(\alpha+\beta+1)}}
#'
#' @param link_alpha Link function for \eqn{\alpha}. Defaults to the log.
#' @param link_beta Link function for \eqn{\beta}. Defaults to the log.
#'
#' @return An S7 object of class [Beta2Distrib()].
#'
#' @seealso [beta1_distrib()]
#'
#' @examples
#' d <- beta2_distrib()
#' theta <- list(alpha = 2, beta = 5)
#' distrib_pdf(d, c(0.1, 0.3, 0.7), theta)
#' c(mean = mean(d, theta), variance = variance(d, theta))
#'
#' # the same law as beta1 with mu = alpha/(alpha+beta), phi = alpha+beta
#' distrib_pdf(beta1_distrib(), 0.3, list(mu = 2 / 7, phi = 7))
#'
#' @export
beta2_distrib <- function(link_alpha = log_link(), link_beta = log_link()) {
  Beta2Distrib(
    distrib_name = "beta2",
    dimension = "univariate",
    bounds = c(0, 1),
    params = c("alpha", "beta"),
    params_interpretation = c(alpha = "shape", beta = "shape"),
    n_params = 2,
    params_bounds = list(alpha = c(0, Inf), beta = c(0, Inf)),
    link_params = list(alpha = link_alpha, beta = link_beta)
  )
}
