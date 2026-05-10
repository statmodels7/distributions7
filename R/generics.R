#' @title Distribution Generics
#' @description A collection of S7 generic functions for mathematical and statistical 
#' operations on probability distributions.
#' @import S7
#' @name distrib_generics
NULL

#' Probability Density Function
#'
#' @description Evaluates the probability density function (PDF) or probability mass function (PMF).
#' @param distrib A distribution object inheriting from the \code{distrib} class.
#' @param ... Additional arguments passed to the specific method (e.g., \code{y}, \code{theta}, \code{log}).
#' @export
distrib_pdf <- S7::new_generic("distrib_pdf", "distrib")

#' Cumulative Distribution Function
#'
#' @description Evaluates the cumulative distribution function (CDF) for a given distribution.
#' @param distrib A distribution object inheriting from the \code{distrib} class.
#' @param ... Additional arguments passed to the specific method (e.g., \code{q}, \code{theta}, \code{lower.tail}).
#' @export
distrib_cdf <- S7::new_generic("distrib_cdf", "distrib")

#' Quantile Function
#'
#' @description Evaluates the quantile function for a given distribution.
#' @param distrib A distribution object inheriting from the \code{distrib} class.
#' @param ... Additional arguments passed to the specific method (e.g., \code{p}, \code{theta}).
#' @export
distrib_quantile <- S7::new_generic("distrib_quantile", "distrib")

#' Random Number Generator
#'
#' @description Generates random variates from the given distribution.
#' @param distrib A distribution object inheriting from the \code{distrib} class.
#' @param ... Additional arguments passed to the specific method (e.g., \code{n}, \code{theta}).
#' @export
distrib_rng <- S7::new_generic("distrib_rng", "distrib")

#' Analytical Gradient
#'
#' @description Computes the analytical first derivatives of the log-likelihood with respect to the distribution's parameters.
#' @param distrib A distribution object inheriting from the \code{distrib} class.
#' @param ... Additional arguments passed to the specific method (e.g., \code{y}, \code{theta}, \code{par}).
#' @export
distrib_gradient <- S7::new_generic("distrib_gradient", "distrib")

#' Analytical Hessian
#'
#' @description Computes the analytical observed second derivatives (Hessian matrix) of the log-likelihood with respect to the distribution's parameters.
#' @param distrib A distribution object inheriting from the \code{distrib} class.
#' @param ... Additional arguments passed to the specific method (e.g., \code{y}, \code{theta}).
#' @export
distrib_hessian <- S7::new_generic("distrib_hessian", "distrib")

#' Analytical Expected Hessian
#'
#' @description Computes the analytical expected second derivatives of the log-likelihood with respect to the distribution's parameters.
#' @param distrib A distribution object inheriting from the \code{distrib} class.
#' @param ... Additional arguments passed to the specific method (e.g., \code{y}, \code{theta}).
#' @export
distrib_expected_hessian <- S7::new_generic("distrib_expected_hessian", "distrib")

#' Generate Random Parameters
#'
#' @description Generates sensible random parameters for a distribution based on its mathematical domain.
#' @param distrib A distribution object inheriting from the \code{distrib} class.
#' @param ... Additional arguments passed to the specific method.
#' @export
generate_random_theta <- S7::new_generic("generate_random_theta", "distrib")