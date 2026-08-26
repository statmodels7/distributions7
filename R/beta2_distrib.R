#' @include distrib.R generics.R
NULL

#' @title Beta Distribution Class, the Two Shapes
#' @name Beta2Distrib
#'
#' @description
#' The S7 class of the beta family on \eqn{(0, 1)} in its canonical
#' parametrization, the two shapes \eqn{\alpha > 0} and \eqn{\beta > 0}, with
#' density \eqn{y^{\alpha-1}(1-y)^{\beta-1}/B(\alpha, \beta)}. It inherits from
#' `continuous_distrib`, so it answers every generic of the `distrib` contract;
#' the eleven methods listed below are registered on it directly and everything
#' else comes from the parent.
#'
#' Build one with [beta2_distrib()], which supplies the two link functions and
#' fills the properties in. This page documents the raw S7 constructor, which
#' takes the parent's properties and validates none of the relationships
#' between them.
#'
#' @inheritParams distrib
#'
#' @return An S7 object of class `Beta2Distrib`, inheriting from
#'   `continuous_distrib` and from `distrib`. Its properties are the parent's:
#'   `distrib_name`, `dimension`, `bounds`, `params`, `params_interpretation`,
#'   `n_params`, `params_bounds`, `link_params` and `params_smooth`. For an
#'   object built by [beta2_distrib()] they hold `"beta2"`, `"univariate"`,
#'   `c(0, 1)`, `c("alpha", "beta")`, the interpretations
#'   `c(alpha = "shape", beta = "shape")`, `2`, the domain \eqn{(0, \infty)}
#'   for both parameters, and the two links.
#'
#' @seealso [beta2_distrib()] to build one;
#'   [beta1_distrib()] for the same law in mean and precision;
#'   [distrib_pdf.Beta2Distrib()] and [distrib_gradient.Beta2Distrib()] for the
#'   closed forms this class supplies.
#'
#' @section Methods:
#' Registered on this class:
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
#'
#' @examples
#' d <- beta2_distrib()
#' S7::S7_inherits(d, continuous_distrib)
#'
#' # The properties a consumer reads to drive the family without knowing it.
#' d@params
#' d@params_interpretation
#' d@bounds
#'
#' # Both parameters are shapes, so neither is a mean: the mean is their ratio.
#' mean(d, list(alpha = 2, beta = 5))
Beta2Distrib <- S7::new_class("Beta2Distrib", parent = continuous_distrib)

# --- S7 METHODS IMPLEMENTATION ---

#' @title Beta Probability Density Function in the Shapes
#' @name distrib_pdf.Beta2Distrib
#' @description
#' Computes the beta density
#' \deqn{f(y; \alpha, \beta) = \dfrac{y^{\alpha-1}(1-y)^{\beta-1}}
#'       {B(\alpha, \beta)}, \qquad 0 < y < 1,}
#' with \eqn{B} the beta function, by calling [stats::dbeta()] at
#' `shape1 = alpha` and `shape2 = beta`. With `log = TRUE` the logarithm is
#' formed inside `dbeta()` and stays finite where the density itself
#' underflows.
#'
#' The density is unbounded at 0 when \eqn{\alpha < 1} and at 1 when
#' \eqn{\beta < 1}; at \eqn{\alpha = \beta = 1} it is the uniform.
#'
#' @param distrib A `Beta2Distrib` object, from [beta2_distrib()].
#' @param y A numeric vector of observations in \eqn{(0, 1)}. A value outside
#'   \eqn{[0, 1]} gives 0, and the endpoints give 0, a finite value or `Inf`
#'   according to the shapes.
#' @param theta A named list with components `alpha` and `beta`, each a numeric
#'   vector of length 1 or of the length of `y`. A component of length 1 is
#'   recycled. Both must be strictly positive.
#' @param log Logical of length 1. When `TRUE` the log-density is returned.
#'   Defaults to `FALSE`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of densities, of length
#'   `max(length(y), length(alpha), length(beta))`, one value per observation.
#'
#' @seealso [distrib_cdf.Beta2Distrib()] for the distribution function,
#'   [distrib_gradient.Beta2Distrib()] for the derivatives of the log-density,
#'   [distrib_pdf.Beta1Distrib()] for the same density in the mean and the
#'   precision, and [distrib_pdf()] for the generic.
#'
#' @examples
#' d <- beta2_distrib()
#' y <- c(0.1, 0.3, 0.7)
#' th <- list(alpha = 2, beta = 5)
#'
#' # The method is stats::dbeta at these two shapes.
#' all.equal(distrib_pdf(d, y, th), dbeta(y, 2, 5))
#'
#' # The same law as beta1 at mu = alpha/(alpha + beta), phi = alpha + beta.
#' all.equal(distrib_pdf(d, y, th),
#'           distrib_pdf(beta1_distrib(), y, list(mu = 2 / 7, phi = 7)))
#'
#' # Both shapes 1 is the uniform.
#' distrib_pdf(d, y, list(alpha = 1, beta = 1))
#'
#' # Near the boundary the density underflows and its logarithm does not.
#' distrib_pdf(d, 1e-40, list(alpha = 20, beta = 5))
#' distrib_pdf(d, 1e-40, list(alpha = 20, beta = 5), log = TRUE)
S7::method(distrib_pdf, Beta2Distrib) <- function(distrib, y, theta, log = FALSE, ...) {
  stats::dbeta(y, shape1 = theta[[1]], shape2 = theta[[2]], log = log)
}

#' @title Beta Cumulative Distribution Function in the Shapes
#' @name distrib_cdf.Beta2Distrib
#' @description
#' Computes the beta distribution function, the regularized incomplete beta
#' function \eqn{F(q) = I_q(\alpha, \beta)}, by calling [stats::pbeta()] at
#' `shape1 = alpha` and `shape2 = beta`. Both tails are available exactly:
#' `lower.tail = FALSE` evaluates \eqn{1 - F} without forming the difference,
#' and `log.p = TRUE` returns a logarithm that stays finite where the
#' probability itself underflows to zero.
#'
#' @param distrib A `Beta2Distrib` object, from [beta2_distrib()].
#' @param q A numeric vector of quantiles. A value at or below 0 gives a
#'   lower-tail probability of 0 and one at or above 1 gives 1.
#' @param theta A named list with components `alpha` and `beta`, each a numeric
#'   vector of length 1 or of the length of `q`. A component of length 1 is
#'   recycled. Both must be strictly positive.
#' @param lower.tail Logical of length 1. When `TRUE`, the default,
#'   probabilities are \eqn{P(Y \le q)}; when `FALSE` they are \eqn{P(Y > q)}.
#' @param log.p Logical of length 1. When `TRUE` the logarithm of the
#'   probability is returned. Defaults to `FALSE`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of probabilities in \eqn{[0, 1]}, of length
#'   `max(length(q), length(alpha), length(beta))`. With `log.p = TRUE` the
#'   values are logarithms and are non-positive.
#'
#' @seealso [distrib_quantile.Beta2Distrib()] for the inverse,
#'   [distrib_pdf.Beta2Distrib()] for the density, [distrib_grad_cdf()] for the
#'   derivatives of this function in the parameters, which the beta takes by
#'   finite difference because the derivative of an incomplete beta in its
#'   shapes is hypergeometric, and [distrib_cdf()] for the generic.
#'
#' @examples
#' d <- beta2_distrib()
#' th <- list(alpha = 2, beta = 5)
#'
#' # The method is stats::pbeta at these two shapes.
#' all.equal(distrib_cdf(d, c(0.1, 0.3, 0.7), th),
#'           pbeta(c(0.1, 0.3, 0.7), 2, 5))
#'
#' # The two tails sum to one.
#' distrib_cdf(d, 0.3, th) + distrib_cdf(d, 0.3, th, lower.tail = FALSE)
#'
#' # Both shapes 1 is the uniform, so F(q) = q.
#' distrib_cdf(d, c(0.1, 0.3, 0.7), list(alpha = 1, beta = 1))
#'
#' # Near the lower boundary the probability underflows and its log does not.
#' distrib_cdf(d, 1e-30, list(alpha = 20, beta = 5))
#' distrib_cdf(d, 1e-30, list(alpha = 20, beta = 5), log.p = TRUE)
S7::method(distrib_cdf, Beta2Distrib) <- function(distrib, q, theta,
                                                   lower.tail = TRUE,
                                                   log.p = FALSE, ...) {
  stats::pbeta(q, shape1 = theta[[1]], shape2 = theta[[2]],
               lower.tail = lower.tail, log.p = log.p)
}

#' @title Beta Quantile Function in the Shapes
#' @name distrib_quantile.Beta2Distrib
#' @description
#' Computes the beta quantile function, the inverse of the regularized
#' incomplete beta function in its argument, by calling [stats::qbeta()] at
#' `shape1 = alpha` and `shape2 = beta`. There is no elementary closed form;
#' `qbeta()` inverts the distribution function numerically. The beta
#' distribution function is strictly increasing on \eqn{(0, 1)}, so the round
#' trip through [distrib_cdf.Beta2Distrib()] returns `p`.
#'
#' @param distrib A `Beta2Distrib` object, from [beta2_distrib()].
#' @param p A numeric vector of probabilities in \eqn{[0, 1]}, or of their
#'   logarithms when `log.p = TRUE`. A value outside the range gives `NaN`
#'   with a warning; the endpoints give 0 and 1.
#' @param theta A named list with components `alpha` and `beta`, each a numeric
#'   vector of length 1 or of the length of `p`. A component of length 1 is
#'   recycled. Both must be strictly positive.
#' @param lower.tail Logical of length 1. When `TRUE`, the default, `p` is
#'   \eqn{P(Y \le q)}; when `FALSE` it is \eqn{P(Y > q)}.
#' @param log.p Logical of length 1. When `TRUE` the values in `p` are read as
#'   logarithms of probabilities. Defaults to `FALSE`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of quantiles in \eqn{[0, 1]}, of length
#'   `max(length(p), length(alpha), length(beta))`.
#'
#' @seealso [distrib_cdf.Beta2Distrib()], which this inverts;
#'   [distrib_rng.Beta2Distrib()], which does not use it;
#'   [distrib_quantile()] for the generic.
#'
#' @examples
#' d <- beta2_distrib()
#' th <- list(alpha = 2, beta = 5)
#'
#' # A central 95 percent interval, asymmetric about the mean of 2/7.
#' distrib_quantile(d, c(0.025, 0.5, 0.975), th)
#'
#' # Exact inverse: the round trip returns the probabilities it was given.
#' p <- c(0.025, 0.5, 0.975)
#' all.equal(distrib_cdf(d, distrib_quantile(d, p, th), th), p)
#'
#' # Both shapes 1 is the uniform, so Q(p) = p.
#' distrib_quantile(d, p, list(alpha = 1, beta = 1))
S7::method(distrib_quantile, Beta2Distrib) <- function(distrib, p, theta,
                                                        lower.tail = TRUE,
                                                        log.p = FALSE, ...) {
  stats::qbeta(p, shape1 = theta[[1]], shape2 = theta[[2]],
               lower.tail = lower.tail, log.p = log.p)
}

#' @title Beta Random Number Generator in the Shapes
#' @name distrib_rng.Beta2Distrib
#' @description
#' Draws `n` independent beta variates by calling [stats::rbeta()] at
#' `shape1 = alpha` and `shape2 = beta`, so the draws come from R's own beta
#' generator and depend on `.Random.seed` in the usual way. The
#' ratio-of-uniforms fallback the base class supplies is bypassed.
#'
#' @param distrib A `Beta2Distrib` object, from [beta2_distrib()].
#' @param n A single positive integer, the number of draws.
#' @param theta A named list with components `alpha` and `beta`, each a numeric
#'   vector of length 1 or of length `n`. A component of length 1 is recycled,
#'   so a vector of length `n` draws one variate per parameter setting. Both
#'   must be strictly positive.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of `n` draws in \eqn{(0, 1)}.
#'
#' @seealso [distrib_quantile.Beta2Distrib()] for the inverse-transform route,
#'   [fit_distrib()] to estimate the parameters back from a sample, and
#'   [distrib_rng()] for the generic.
#'
#' @examples
#' d <- beta2_distrib()
#'
#' # Same generator as stats::rbeta at these two shapes.
#' set.seed(2)
#' a <- distrib_rng(d, 3, list(alpha = 2, beta = 5))
#' set.seed(2)
#' identical(a, rbeta(3, 2, 5))
#'
#' # The sample moments recover the shapes through the concentration
#' # k = mean(1 - mean)/var - 1, with alpha = k mean and beta = k(1 - mean).
#' set.seed(8)
#' z <- distrib_rng(d, 2e4, list(alpha = 2, beta = 5))
#' m <- mean(z)
#' k <- m * (1 - m) / var(z) - 1
#' c(alpha = k * m, beta = k * (1 - m))
S7::method(distrib_rng, Beta2Distrib) <- function(distrib, n, theta, ...) {
  stats::rbeta(n, shape1 = theta[[1]], shape2 = theta[[2]])
}

#' @title Beta Score in the Shapes
#' @name distrib_gradient.Beta2Distrib
#' @description
#' Computes the first derivatives of the beta log-density with respect to
#' \eqn{\alpha} and \eqn{\beta}, one value per observation, in closed form:
#' \deqn{\dfrac{\partial \ell}{\partial \alpha} = \log y - \psi(\alpha)
#'         + \psi(\alpha+\beta), \qquad
#'       \dfrac{\partial \ell}{\partial \beta} = \log(1-y) - \psi(\beta)
#'         + \psi(\alpha+\beta),}
#' with \eqn{\psi} the digamma function. The beta is an exponential family in
#' the shapes with sufficient statistics \eqn{\log y} and \eqn{\log(1-y)}, so
#' each component is the corresponding statistic minus its expectation. That
#' is also why every derivative beyond this one is free of the response.
#'
#' The value is computed in plain R, this family carrying no compiled kernel,
#' and \eqn{\log(1-y)} is formed with [base::log1p()] so that it stays accurate
#' at `y` near zero.
#'
#' With `scale = "link"` the generic applies the chain rule for the links the
#' family carries before returning. This method always returns the parameter
#' scale.
#'
#' @param distrib A `Beta2Distrib` object, from [beta2_distrib()].
#' @param y A numeric vector of observations in \eqn{(0, 1)}. An endpoint makes
#'   a logarithm infinite and the score non-finite.
#' @param theta A named list with components `alpha` and `beta`, each a numeric
#'   vector of length 1 or of the length of `y`. A component of length 1 is
#'   recycled. Both must be strictly positive.
#' @param scale One of `"parameter"` (the default) or `"link"`, matched by
#'   [base::match.arg()]. Read by the generic, not by this method.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A named list of two numeric vectors, `alpha` and `beta`, each of
#'   length `max(length(y), length(alpha), length(beta))`.
#'
#' @section Notation:
#' \eqn{\ell} is the log-density of one observation and \eqn{\alpha, \beta > 0}
#' the two shapes. \eqn{\psi} is the digamma function,
#' \eqn{\psi = (\log\Gamma)'}.
#'
#' @seealso [distrib_hessian.Beta2Distrib()] for the second derivatives,
#'   [distrib_gradient.Beta1Distrib()] for the same score in the mean and the
#'   precision, [beta2_higher()] for the orders above, and
#'   [distrib_gradient()] for the generic.
#'
#' @examples
#' d <- beta2_distrib()
#' y <- c(0.1, 0.3, 0.7)
#' th <- list(alpha = 2, beta = 5)
#' g <- distrib_gradient(d, y, th)
#'
#' # The two closed forms, written out with the digamma function.
#' all.equal(g$alpha, log(y) - digamma(2) + digamma(7))
#' all.equal(g$beta, log1p(-y) - digamma(5) + digamma(7))
#'
#' # Each component is a sufficient statistic minus its expectation, so the
#' # sample mean of log(y) matches psi(alpha) - psi(alpha + beta).
#' set.seed(8)
#' z <- distrib_rng(d, 2e5, th)
#' c(sample = mean(log(z)), theory = digamma(2) - digamma(7))
#'
#' # Summed over a fitted sample the score is at the optimizer's tolerance.
#' zz <- distrib_rng(d, 2000, th)
#' fit <- fit_distrib(d, zz)
#' vapply(distrib_gradient(d, zz, as.list(coef(fit))), sum, numeric(1))
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
#' Returns the derivatives of the beta log-density in the two shapes at order
#' 2, 3 or 4. Each component is a difference of polygamma functions at
#' \eqn{\alpha}, \eqn{\beta} and \eqn{\alpha+\beta}: with \eqn{k = }`order - 1`,
#' a component naming only \eqn{\alpha} is
#' \eqn{\psi^{(k)}(\alpha+\beta) - \psi^{(k)}(\alpha)}, one naming only
#' \eqn{\beta} is \eqn{\psi^{(k)}(\alpha+\beta) - \psi^{(k)}(\beta)}, and any
#' mixed component is \eqn{\psi^{(k)}(\alpha+\beta)}.
#'
#' @details
#' The data enter the log-density only through
#' \eqn{(\alpha-1)\log y + (\beta-1)\log(1-y)}, which is linear in the two
#' parameters, so the second derivative already kills it. Every order from two
#' upwards is therefore free of the response, the observed and the expected
#' derivatives are the same numbers, and no expectation is computed anywhere in
#' this family. The three polygamma values are evaluated once each and recycled.
#'
#' @param theta A named list with components `alpha` and `beta`, each a numeric
#'   vector. Both must be strictly positive.
#' @param n The length to recycle each component to, normally `length(y)`.
#' @param order The derivative order, `2L`, `3L` or `4L`. Any other value falls
#'   through to the fourth-order branch.
#'
#' @return A named list of numeric vectors, each of length `n`: three
#'   components at order 2 (`alpha_alpha`, `alpha_beta`, `beta_beta`), four at
#'   order 3 and five at order 4, named for the distinct multi-indices.
#'
#' @seealso [beta2_distrib()] for the family, and
#'   [distrib_hessian.Beta2Distrib()], [distrib_deriv3.Beta2Distrib()] and
#'   [distrib_deriv4.Beta2Distrib()] for the methods that call it.
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

#' @title Beta Observed Hessian in the Shapes
#' @name distrib_hessian.Beta2Distrib
#' @description
#' Computes the three distinct second derivatives of the beta log-density with
#' respect to \eqn{\alpha} and \eqn{\beta}, in closed form:
#' \deqn{\ell^{(\alpha\alpha)} = \psi'(\alpha+\beta) - \psi'(\alpha), \qquad
#'       \ell^{(\alpha\beta)} = \psi'(\alpha+\beta), \qquad
#'       \ell^{(\beta\beta)} = \psi'(\alpha+\beta) - \psi'(\beta),}
#' with \eqn{\psi'} the trigamma function. **All three are free of the
#' response**, the data entering the log-density only through a term linear in
#' the parameters. The values are constant within a parameter setting and are
#' recycled to the length of `y`, and they equal
#' [distrib_expected_hessian.Beta2Distrib()] exactly.
#'
#' @param distrib A `Beta2Distrib` object, from [beta2_distrib()].
#' @param y A numeric vector of observations. Only its length is used.
#' @param theta A named list with components `alpha` and `beta`, each a numeric
#'   vector of length 1 or of the length of `y`. Both must be strictly
#'   positive.
#' @param scale One of `"parameter"` (the default) or `"link"`, matched by
#'   [base::match.arg()]. Read by the generic, not by this method.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A named list of three numeric vectors, `alpha_alpha`, `beta_beta`
#'   and `alpha_beta`, in that order, each of length `length(y)`. The three
#'   name the distinct entries of a symmetric \eqn{2 \times 2} matrix per
#'   observation.
#'
#' @section Notation:
#' \eqn{\ell^{(ij)}} is the second derivative of the log-density in parameters
#' \eqn{i} and \eqn{j}; parenthesized superscripts name derivatives.
#' \eqn{\psi'} is the trigamma function.
#'
#' @seealso [distrib_gradient.Beta2Distrib()] for the score, which is the one
#'   quantity here that does read the data;
#'   [distrib_expected_hessian.Beta2Distrib()], which returns the same numbers;
#'   [beta2_higher()], which computes this; and [distrib_hessian()] for the
#'   generic.
#'
#' @examples
#' d <- beta2_distrib()
#' y <- c(0.1, 0.3, 0.7)
#' th <- list(alpha = 2, beta = 5)
#' h <- distrib_hessian(d, y, th)
#' h
#'
#' # Written out with the trigamma function.
#' c(trigamma(7) - trigamma(2), trigamma(7) - trigamma(5), trigamma(7))
#'
#' # Free of the response, so it equals its own expectation to the bit.
#' identical(h, distrib_expected_hessian(d, y, th))
#'
#' # It is the second derivative of the log-density, so a central difference
#' # of the score reproduces it.
#' eps <- 1e-6
#' up <- distrib_gradient(d, y, list(alpha = 2 + eps, beta = 5))$alpha
#' dn <- distrib_gradient(d, y, list(alpha = 2 - eps, beta = 5))$alpha
#' all.equal((up - dn) / (2 * eps), h$alpha_alpha, tolerance = 1e-5)
S7::method(distrib_hessian, Beta2Distrib) <- function(distrib, y, theta,
                                                       scale = c("parameter", "link"), ...) {
  h <- beta2_higher(theta, length(y), 2L)
  h[hess_names(distrib@params)]
}

#' @title Beta Expected Hessian in the Shapes
#' @name distrib_expected_hessian.Beta2Distrib
#' @description
#' Returns the same three numbers as [distrib_hessian.Beta2Distrib()], the
#' observed Hessian being free of the response and so equal to its own
#' expectation:
#' \deqn{\mathbb{E}\left[\ell^{(\alpha\alpha)}\right] = \psi'(\alpha+\beta)
#'         - \psi'(\alpha), \qquad
#'       \mathbb{E}\left[\ell^{(\alpha\beta)}\right] = \psi'(\alpha+\beta),
#'       \qquad
#'       \mathbb{E}\left[\ell^{(\beta\beta)}\right] = \psi'(\alpha+\beta)
#'         - \psi'(\beta).}
#' Nothing is averaged, integrated or simulated, so `approx` and `nsim` are
#' ignored and `y` is read only for its length.
#'
#' The consequence for fitting is that Fisher scoring and Newton's method take
#' the same step on the parameter scale, both inverting the same matrix. On the
#' **link** scale they differ, the chain rule there adding a term in the score,
#' which does read the data.
#'
#' The mixed entry \eqn{\psi'(\alpha+\beta)} is positive and never zero, so the
#' two shapes are not orthogonal at any parameter setting.
#'
#' @param distrib A `Beta2Distrib` object, from [beta2_distrib()].
#' @param y A numeric vector of observations. Only its length is used.
#' @param theta A named list with components `alpha` and `beta`, each a numeric
#'   vector of length 1 or of the length of `y`. Both must be strictly
#'   positive.
#' @param scale One of `"parameter"` (the default) or `"link"`, matched by
#'   [base::match.arg()]. Read by the generic, not by this method.
#' @param approx Ignored, the expectation being exact. Accepted so that the
#'   signature matches the generic's, where it selects between the Bartlett,
#'   quadrature, Monte Carlo and outer-product routes.
#' @param nsim Ignored, for the same reason. Defaults to `10000`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A named list of three numeric vectors, `alpha_alpha`, `beta_beta`
#'   and `alpha_beta`, in that order, each of length `length(y)`.
#'
#' @section Notation:
#' The **expected information** is
#' \eqn{\mathbb{E}[-\partial^2\ell/\partial\theta\,\partial\theta^\top]}, the
#' expectation of the **observed information** under the model. The beta is a
#' regular family, so the second Bartlett identity holds and this equals the
#' variance of the score.
#'
#' @seealso [distrib_hessian.Beta2Distrib()], which returns the same numbers;
#'   [distrib_expected_hessian.Beta1Distrib()], where the same law in the mean
#'   and the precision does have an observed Hessian that reads the data;
#'   [fisher_scoring()]; and [distrib_expected_hessian()] for the generic.
#'
#' @examples
#' d <- beta2_distrib()
#' y <- c(0.1, 0.3, 0.7)
#' th <- list(alpha = 2, beta = 5)
#' lapply(distrib_expected_hessian(d, y, th), unique)
#'
#' # Identical to the observed Hessian, so the two fitting methods coincide.
#' identical(distrib_expected_hessian(d, y, th), distrib_hessian(d, y, th))
#'
#' set.seed(8)
#' z <- distrib_rng(d, 2000, th)
#' rbind(newton = coef(fit_distrib(d, z, method = "newton")),
#'       fisher = coef(fit_distrib(d, z, method = "fisher")))
#'
#' # The mixed entry is psi'(alpha + beta), positive at every setting, so the
#' # two shapes are never orthogonal.
#' distrib_expected_hessian(d, 0.5, th)$alpha_beta
S7::method(distrib_expected_hessian, Beta2Distrib) <- function(distrib, y, theta,
                                                                scale = c("parameter", "link"),
                                                                approx = c("bartlett", "integrate", "mc", "opg"),
                                                                nsim = 10000, ...) {
  h <- beta2_higher(theta, length(y), 2L)
  h[hess_names(distrib@params)]
}

#' @title Beta Third-Order Derivatives in the Shapes
#' @name distrib_deriv3.Beta2Distrib
#' @description
#' Computes the four distinct third derivatives of the beta log-density with
#' respect to \eqn{\alpha} and \eqn{\beta}, in closed form:
#' \deqn{\ell^{(\alpha\alpha\alpha)} = \psi''(\alpha+\beta) - \psi''(\alpha),
#'       \qquad
#'       \ell^{(\alpha\alpha\beta)} = \ell^{(\alpha\beta\beta)} =
#'         \psi''(\alpha+\beta), \qquad
#'       \ell^{(\beta\beta\beta)} = \psi''(\alpha+\beta) - \psi''(\beta),}
#' with \eqn{\psi''} the second derivative of the digamma function. The two
#' mixed components are equal, every mixed derivative of \eqn{\log
#' B(\alpha,\beta)} of a given order being the same polygamma of the sum.
#'
#' The values are free of the response, so `expected` selects nothing: the same
#' computation runs either way and the two results are identical to the bit.
#' `approx` and `nsim` are ignored for the same reason.
#'
#' @param distrib A `Beta2Distrib` object, from [beta2_distrib()].
#' @param y A numeric vector of observations. Only its length is used.
#' @param theta A named list with components `alpha` and `beta`, each a numeric
#'   vector of length 1 or of the length of `y`. Both must be strictly
#'   positive.
#' @param expected Logical of length 1, and without effect here, the observed
#'   and expected third derivatives being the same numbers. Defaults to
#'   `FALSE`.
#' @param scale One of `"parameter"` (the default) or `"link"`, matched by
#'   [base::match.arg()]. Read by the generic, not by this method.
#' @param approx Ignored, a closed form being available.
#' @param nsim Ignored, for the same reason. Defaults to `10000`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A named list of four numeric vectors, `alpha_alpha_alpha`,
#'   `alpha_alpha_beta`, `alpha_beta_beta` and `beta_beta_beta`, each of length
#'   `length(y)` and constant within itself when the parameters are.
#'
#' @section Notation:
#' \eqn{\ell^{(ijk)}} is the third derivative of the log-density in parameters
#' \eqn{i}, \eqn{j} and \eqn{k}. \eqn{\psi} is the digamma function and
#' \eqn{\psi^{(m)}} its \eqn{m}th derivative.
#'
#' @seealso [distrib_hessian.Beta2Distrib()] for the order below,
#'   [distrib_deriv4.Beta2Distrib()] for the order above, [beta2_higher()],
#'   which computes this, and [distrib_deriv3()] for the generic.
#'
#' @examples
#' d <- beta2_distrib()
#' y <- c(0.1, 0.3, 0.7)
#' th <- list(alpha = 2, beta = 5)
#' d3 <- distrib_deriv3(d, y, th)
#'
#' # Four constants, and the two mixed components are the same number.
#' lapply(d3, unique)
#' psigamma(7, deriv = 2)
#'
#' # Free of the response, so asking for the expectation changes nothing.
#' identical(d3, distrib_deriv3(d, y, th, expected = TRUE))
#'
#' # A central difference of the Hessian reproduces the same component.
#' eps <- 1e-6
#' up <- distrib_hessian(d, y, list(alpha = 2 + eps, beta = 5))$alpha_alpha
#' dn <- distrib_hessian(d, y, list(alpha = 2 - eps, beta = 5))$alpha_alpha
#' all.equal((up - dn) / (2 * eps), d3$alpha_alpha_alpha, tolerance = 1e-4)
S7::method(distrib_deriv3, Beta2Distrib) <- function(distrib, y, theta, expected = FALSE,
                                                      scale = c("parameter", "link"),
                                                      approx = c("integrate", "bartlett", "mc", "opg"),
                                                      nsim = 10000, ...) {
  beta2_higher(theta, length(y), 3L)
}

#' @title Beta Fourth-Order Derivatives in the Shapes
#' @name distrib_deriv4.Beta2Distrib
#' @description
#' Computes the five distinct fourth derivatives of the beta log-density with
#' respect to \eqn{\alpha} and \eqn{\beta}, in closed form:
#' \deqn{\ell^{(\alpha^4)} = \psi'''(\alpha+\beta) - \psi'''(\alpha), \qquad
#'       \ell^{(\beta^4)} = \psi'''(\alpha+\beta) - \psi'''(\beta),}
#' and every one of the three mixed components equal to
#' \eqn{\psi'''(\alpha+\beta)}, with \eqn{\psi'''} the third derivative of the
#' digamma function.
#'
#' As at third order the values are free of the response, so `expected`,
#' `approx` and `nsim` are all without effect and the result is identical to
#' the bit whichever is asked for.
#'
#' @param distrib A `Beta2Distrib` object, from [beta2_distrib()].
#' @param y A numeric vector of observations. Only its length is used.
#' @param theta A named list with components `alpha` and `beta`, each a numeric
#'   vector of length 1 or of the length of `y`. Both must be strictly
#'   positive.
#' @param expected Logical of length 1, and without effect here, the observed
#'   and expected fourth derivatives being the same numbers. Defaults to
#'   `FALSE`.
#' @param scale One of `"parameter"` (the default) or `"link"`, matched by
#'   [base::match.arg()]. Read by the generic, not by this method.
#' @param approx Ignored, a closed form being available.
#' @param nsim Ignored, for the same reason. Defaults to `10000`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A named list of five numeric vectors,
#'   `alpha_alpha_alpha_alpha`, `alpha_alpha_alpha_beta`,
#'   `alpha_alpha_beta_beta`, `alpha_beta_beta_beta` and `beta_beta_beta_beta`,
#'   each of length `length(y)` and constant within itself when the parameters
#'   are.
#'
#' @section Notation:
#' \eqn{\ell^{(ijkl)}} is the fourth derivative of the log-density in
#' parameters \eqn{i}, \eqn{j}, \eqn{k} and \eqn{l}. \eqn{\psi} is the digamma
#' function and \eqn{\psi^{(m)}} its \eqn{m}th derivative.
#'
#' @seealso [distrib_deriv3.Beta2Distrib()] for the order below,
#'   [beta2_higher()], which computes this, and [distrib_deriv4()] for the
#'   generic.
#'
#' @examples
#' d <- beta2_distrib()
#' y <- c(0.1, 0.3, 0.7)
#' th <- list(alpha = 2, beta = 5)
#' d4 <- distrib_deriv4(d, y, th)
#'
#' # Five constants, and the three mixed components are one number.
#' lapply(d4, unique)
#' psigamma(7, deriv = 3)
#'
#' # Free of the response, so asking for the expectation changes nothing.
#' identical(d4, distrib_deriv4(d, y, th, expected = TRUE))
#'
#' # A central difference of the third order reproduces it.
#' eps <- 1e-6
#' up <- distrib_deriv3(d, y, list(alpha = 2 + eps, beta = 5))$alpha_alpha_alpha
#' dn <- distrib_deriv3(d, y, list(alpha = 2 - eps, beta = 5))$alpha_alpha_alpha
#' all.equal((up - dn) / (2 * eps), d4$alpha_alpha_alpha_alpha,
#'           tolerance = 1e-3)
S7::method(distrib_deriv4, Beta2Distrib) <- function(distrib, y, theta, expected = FALSE,
                                                      scale = c("parameter", "link"),
                                                      approx = c("integrate", "bartlett", "mc", "opg"),
                                                      nsim = 10000, ...) {
  beta2_higher(theta, length(y), 4L)
}

#' @title Beta First Derivative in the Response, the Shapes
#' @name distrib_grad_y.Beta2Distrib
#' @description
#' Computes the first derivative of the beta log-density with respect to the
#' response, in closed form:
#' \deqn{\dfrac{\partial \ell}{\partial y} = \dfrac{\alpha - 1}{y}
#'       - \dfrac{\beta - 1}{1 - y}.}
#' Where both shapes exceed one it changes sign at the mode
#' \eqn{(\alpha-1)/(\alpha+\beta-2)}, so it is positive below the mode and
#' negative above it. Where a shape falls below one the density is unbounded at
#' the corresponding endpoint and the derivative does not change sign there. At
#' \eqn{\alpha = \beta = 1} it is identically zero, the density being the
#' uniform.
#'
#' Quantile residuals and the mixed derivatives of [distrib_cross_y()] read it.
#'
#' @param distrib A `Beta2Distrib` object, from [beta2_distrib()].
#' @param y A numeric vector of observations in \eqn{(0, 1)}. An endpoint makes
#'   the value infinite unless the corresponding shape is exactly 1.
#' @param theta A named list with components `alpha` and `beta`, each a numeric
#'   vector of length 1 or of the length of `y`. A component of length 1 is
#'   recycled. Both must be strictly positive.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of length
#'   `max(length(y), length(alpha), length(beta))`, one value per observation.
#'
#' @seealso [distrib_hess_y.Beta2Distrib()] for the second derivative in the
#'   response, [distrib_gradient.Beta2Distrib()] for the score in the
#'   parameters, and [distrib_grad_y()] for the generic.
#'
#' @examples
#' d <- beta2_distrib()
#' y <- c(0.1, 0.3, 0.7)
#' th <- list(alpha = 2, beta = 5)
#'
#' all.equal(distrib_grad_y(d, y, th), (2 - 1) / y - (5 - 1) / (1 - y))
#'
#' # Zero at the mode (alpha - 1)/(alpha + beta - 2).
#' mode <- (2 - 1) / (2 + 5 - 2)
#' c(mode = mode, at_mode = distrib_grad_y(d, mode, th))
#'
#' # Identically zero at alpha = beta = 1, where the density is the uniform.
#' distrib_grad_y(d, y, list(alpha = 1, beta = 1))
#'
#' # It is the derivative of the log-density, so a central difference in y
#' # reproduces it.
#' eps <- 1e-6
#' (distrib_pdf(d, y + eps, th, log = TRUE) -
#'   distrib_pdf(d, y - eps, th, log = TRUE)) / (2 * eps)
S7::method(distrib_grad_y, Beta2Distrib) <- function(distrib, y, theta, ...) {
  (theta[[1]] - 1) / y - (theta[[2]] - 1) / (1 - y)
}

#' @title Beta Second Derivative in the Response, the Shapes
#' @name distrib_hess_y.Beta2Distrib
#' @description
#' Computes the second derivative of the beta log-density with respect to the
#' response, in closed form:
#' \deqn{\dfrac{\partial^2 \ell}{\partial y^2} = -\dfrac{\alpha - 1}{y^2}
#'       - \dfrac{\beta - 1}{(1 - y)^2}.}
#' It is negative throughout, so the log-density is concave in the response,
#' whenever both shapes are at least one; a shape below one makes the
#' corresponding term positive and can turn the curvature positive near that
#' endpoint. At \eqn{\alpha = \beta = 1} it is exactly zero, the density being
#' the uniform.
#'
#' @param distrib A `Beta2Distrib` object, from [beta2_distrib()].
#' @param y A numeric vector of observations in \eqn{(0, 1)}. An endpoint makes
#'   the value infinite unless the corresponding shape is exactly 1.
#' @param theta A named list with components `alpha` and `beta`, each a numeric
#'   vector of length 1 or of the length of `y`. A component of length 1 is
#'   recycled. Both must be strictly positive.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of length
#'   `max(length(y), length(alpha), length(beta))`, one value per observation.
#'
#' @seealso [distrib_grad_y.Beta2Distrib()] for the first derivative in the
#'   response, [distrib_hessian.Beta2Distrib()] for the curvature in the
#'   parameters, and [distrib_hess_y()] for the generic.
#'
#' @examples
#' d <- beta2_distrib()
#' y <- c(0.1, 0.3, 0.7)
#' th <- list(alpha = 2, beta = 5)
#'
#' all.equal(distrib_hess_y(d, y, th),
#'           -(2 - 1) / y^2 - (5 - 1) / (1 - y)^2)
#'
#' # Concave everywhere while both shapes exceed one.
#' all(distrib_hess_y(d, y, th) < 0)
#'
#' # Exactly flat at alpha = beta = 1, where the beta is the uniform.
#' distrib_hess_y(d, y, list(alpha = 1, beta = 1))
S7::method(distrib_hess_y, Beta2Distrib) <- function(distrib, y, theta, ...) {
  -(theta[[1]] - 1) / y^2 - (theta[[2]] - 1) / (1 - y)^2
}


#' Beta Distribution, the Two Shapes
#'
#' @description
#' Builds the distribution object for the beta family on \eqn{(0, 1)} in its
#' canonical parametrization, the two shapes \eqn{\alpha > 0} and
#' \eqn{\beta > 0}. The returned object carries closed-form derivatives of the
#' log-density to fourth order, in the parameters and in the response, and
#' closed-form moments, so every generic of the toolkit answers without a
#' numerical fallback.
#'
#' The two arguments choose the links that carry each parameter to the
#' unconstrained scale an optimizer works on. Both default to the logarithm,
#' both shapes being positive.
#'
#' @param link_alpha A `link` object from `linkfunctions7` for the first shape
#'   \eqn{\alpha}. Defaults to [linkfunctions7::log_link()], which maps
#'   \eqn{(0, \infty)} onto the line and so keeps every fitted value positive.
#' @param link_beta A `link` object from `linkfunctions7` for the second shape
#'   \eqn{\beta}. Defaults to [linkfunctions7::log_link()], for the same
#'   reason.
#'
#' @details
#' # The parametrization
#'
#' The density on \eqn{y \in (0, 1)} is
#' \deqn{f(y; \alpha, \beta) = \dfrac{y^{\alpha-1}(1-y)^{\beta-1}}
#'       {B(\alpha, \beta)},}
#' with \eqn{B} the beta function, the distribution function
#' \eqn{F(q) = I_q(\alpha, \beta)} the regularized incomplete beta function,
#' and the quantile function its numerical inverse. The mean is
#' \eqn{\alpha/(\alpha+\beta)} and the variance
#' \eqn{\alpha\beta/\{(\alpha+\beta)^2(\alpha+\beta+1)\}}. At
#' \eqn{\alpha = \beta = 1} the density is the uniform, and where either shape
#' falls below one the density is unbounded at the corresponding endpoint.
#'
#' This is the same law as [beta1_distrib()], which carries the mean and a
#' precision, with \eqn{\alpha = \mu\varphi} and \eqn{\beta = (1-\mu)\varphi}.
#' The mean parametrization is the one a regression wants; this one is how the
#' family is usually written, and it is what a conjugate analysis produces, the
#' beta being conjugate for a binomial probability.
#'
#' # Where the response stops entering
#'
#' The beta is an exponential family in the shapes, with sufficient statistics
#' \eqn{\log y} and \eqn{\log(1-y)}, so the score is
#' \deqn{\dfrac{\partial \ell}{\partial \alpha} = \log y - \psi(\alpha)
#'         + \psi(\alpha+\beta), \qquad
#'       \dfrac{\partial \ell}{\partial \beta} = \log(1-y) - \psi(\beta)
#'         + \psi(\alpha+\beta),}
#' each a statistic minus its expectation. The data enter only through
#' \eqn{(\alpha-1)\log y + (\beta-1)\log(1-y)}, which is linear in the two
#' parameters, so the second derivative already kills it:
#' \deqn{\ell^{(\alpha\alpha)} = \psi'(\alpha+\beta) - \psi'(\alpha), \qquad
#'       \ell^{(\alpha\beta)} = \psi'(\alpha+\beta), \qquad
#'       \ell^{(\beta\beta)} = \psi'(\alpha+\beta) - \psi'(\beta).}
#'
#' Three consequences follow. The observed and expected Hessians are the same
#' matrix, so Fisher scoring and Newton's method take the same step on the
#' parameter scale, and asking any method for `expected = TRUE` returns the
#' same numbers. Every third and fourth derivative is likewise a difference of
#' polygamma functions and free of the data. And the mixed entry
#' \eqn{\psi'(\alpha+\beta)} is positive at every setting, so the two shapes
#' are never orthogonal.
#'
#' The derivatives of the *distribution* function in the parameters have no
#' elementary form, the derivative of an incomplete beta in its shapes being
#' hypergeometric, and are taken by finite difference on the analytic cdf.
#'
#' # Estimation
#'
#' [fit_distrib()] maximizes the log-likelihood on the link scale. Neither
#' estimate is closed form, both estimating equations involving digamma
#' functions. The method of moments supplies the starting values through the
#' concentration \eqn{k = \bar y(1-\bar y)/s^2 - 1}, with
#' \eqn{\hat\alpha = k\bar y} and \eqn{\hat\beta = k(1-\bar y)}, and the
#' example below shows them landing beside the estimates.
#'
#' @section Notation:
#' \eqn{\ell} is the log-density of one observation and \eqn{\alpha, \beta > 0}
#' the two shapes. \eqn{\psi} is the digamma function and \eqn{\psi^{(m)}} its
#' \eqn{m}th derivative. \eqn{\eta} is a parameter on the unconstrained scale
#' of its link, with \eqn{\theta = g^{-1}(\eta)}.
#'
#' @return An S7 object of class `Beta2Distrib`, inheriting from
#'   `continuous_distrib`, with `distrib_name` `"beta2"`, `dimension`
#'   `"univariate"`, `bounds` `c(0, 1)`, `params` `c("alpha", "beta")`,
#'   `n_params` `2`, `params_bounds` the domain \eqn{(0, \infty)} for both, and
#'   `link_params` the two links given here.
#'
#' @seealso
#' [beta1_distrib()] for the same law in the mean and a precision, which is the
#' parametrization a regression uses; [betabinom2_distrib()] for the beta as a
#' mixing law over a binomial probability; [gamma1_distrib()] for the family a
#' ratio of gammas turns into this one; [fit_distrib()] to estimate the
#' parameters; [check_distrib()] to validate a family of your own against the
#' same battery this one passes; [Beta2Distrib] for the class.
#'
#' @references
#' Johnson, N. L., Kotz, S. and Balakrishnan, N. (1995).
#' *Continuous Univariate Distributions*, Volume 2, 2nd edition, Chapter 25.
#' Wiley, New York.
#'
#' @examples
#' d <- beta2_distrib()
#' d
#'
#' # The density is stats::dbeta at these two shapes.
#' y <- c(0.1, 0.3, 0.7)
#' th <- list(alpha = 2, beta = 5)
#' all.equal(distrib_pdf(d, y, th), dbeta(y, 2, 5))
#'
#' # Moments: the mean is the ratio, and neither parameter is a mean.
#' c(mean = mean(d, th), var = variance(d, th),
#'   skew = skewness(d, th), kurt = kurtosis(d, th))
#' c(2 / 7, 2 * 5 / (7^2 * 8))
#'
#' # The same law as beta1 at mu = alpha/(alpha + beta), phi = alpha + beta.
#' all.equal(distrib_pdf(d, y, th),
#'           distrib_pdf(beta1_distrib(), y, list(mu = 2 / 7, phi = 7)))
#'
#' # Fitting recovers the shapes; the moment estimates start it off.
#' set.seed(8)
#' z <- distrib_rng(d, 2000, th)
#' fit <- fit_distrib(d, z)
#' m <- mean(z)
#' k <- m * (1 - m) / var(z) - 1
#' rbind(fitted  = coef(fit),
#'       moments = c(alpha = k * m, beta = k * (1 - m)))
#'
#' # Every derivative past the first is free of the response, so Newton and
#' # Fisher scoring invert the same matrix and reach the same point.
#' rbind(newton = coef(fit_distrib(d, z, method = "newton")),
#'       fisher = coef(fit_distrib(d, z, method = "fisher")))
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
