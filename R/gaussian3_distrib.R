#' @include distrib.R generics.R
NULL

#' @title S7 Class for the Gaussian Distribution in Mean and Precision
#' @name Gaussian3Distrib
#'
#' @description A subclass of \code{continuous_distrib} for the Gaussian
#'   written in its mean and its \strong{precision}.
#' @inheritParams distrib
#' @return An object of class \code{Gaussian3Distrib}.
#' @seealso \code{\link{gaussian3_distrib}}, \code{\link{gaussian1_distrib}}
#'
#' @section Methods:
#' Methods implemented for this class:
#'   \code{\link[=distrib_cdf.Gaussian3Distrib]{distrib_cdf()}},
#'   \code{\link[=distrib_deriv3.Gaussian3Distrib]{distrib_deriv3()}},
#'   \code{\link[=distrib_deriv4.Gaussian3Distrib]{distrib_deriv4()}},
#'   \code{\link[=distrib_expected_hessian.Gaussian3Distrib]{distrib_expected_hessian()}},
#'   \code{\link[=distrib_grad_y.Gaussian3Distrib]{distrib_grad_y()}},
#'   \code{\link[=distrib_gradient.Gaussian3Distrib]{distrib_gradient()}},
#'   \code{\link[=distrib_hess_y.Gaussian3Distrib]{distrib_hess_y()}},
#'   \code{\link[=distrib_hessian.Gaussian3Distrib]{distrib_hessian()}},
#'   \code{\link[=distrib_pdf.Gaussian3Distrib]{distrib_pdf()}},
#'   \code{\link[=distrib_quantile.Gaussian3Distrib]{distrib_quantile()}},
#'   \code{\link[=distrib_rng.Gaussian3Distrib]{distrib_rng()}}
#'
#' Everything else is inherited from \code{\link{continuous_distrib}}.
Gaussian3Distrib <- S7::new_class("Gaussian3Distrib", parent = continuous_distrib)

# --- S7 METHODS IMPLEMENTATION ---

#' @title Gaussian Density in Mean and Precision
#' @name distrib_pdf.Gaussian3Distrib
#' @description
#' \deqn{f(y) = \sqrt{\dfrac{\tau}{2\pi}}
#'       \exp\left\{-\dfrac{\tau(y-\mu)^2}{2}\right\}}
#' @param distrib A \code{Gaussian3Distrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list with \code{mu} and \code{tau}.
#' @param log Logical; if \code{TRUE}, returns the log-density.
#' @return A numeric vector.
#' @seealso \code{\link{gaussian3_distrib}}
S7::method(distrib_pdf, Gaussian3Distrib) <- function(distrib, y, theta, log = FALSE) {
  stats::dnorm(y, mean = theta[[1]], sd = 1 / sqrt(theta[[2]]), log = log)
}

#' @title Gaussian Distribution Function in Mean and Precision
#' @name distrib_cdf.Gaussian3Distrib
#' @description The normal distribution function at standard deviation
#'   \eqn{1/\sqrt{\tau}}.
#' @param distrib A \code{Gaussian3Distrib} object.
#' @param q A numeric vector of quantiles.
#' @param theta A list with \code{mu} and \code{tau}.
#' @param lower.tail Logical; if \code{TRUE}, probabilities are \eqn{P(Y \le q)}.
#' @param log.p Logical; if \code{TRUE}, returns log-probabilities.
#' @param ... Unused.
#' @return A numeric vector.
#' @seealso \code{\link{gaussian3_distrib}}
S7::method(distrib_cdf, Gaussian3Distrib) <- function(distrib, q, theta,
                                                       lower.tail = TRUE,
                                                       log.p = FALSE, ...) {
  stats::pnorm(q, mean = theta[[1]], sd = 1 / sqrt(theta[[2]]),
               lower.tail = lower.tail, log.p = log.p)
}

#' @title Gaussian Quantile Function in Mean and Precision
#' @name distrib_quantile.Gaussian3Distrib
#' @description The normal quantile function at standard deviation
#'   \eqn{1/\sqrt{\tau}}.
#' @param distrib A \code{Gaussian3Distrib} object.
#' @param p A numeric vector of probabilities.
#' @param theta A list with \code{mu} and \code{tau}.
#' @param lower.tail Logical; if \code{TRUE}, probabilities are \eqn{P(Y \le q)}.
#' @param log.p Logical; if \code{TRUE}, \code{p} is given as a log-probability.
#' @param ... Unused.
#' @return A numeric vector.
#' @seealso \code{\link{gaussian3_distrib}}
S7::method(distrib_quantile, Gaussian3Distrib) <- function(distrib, p, theta,
                                                            lower.tail = TRUE,
                                                            log.p = FALSE, ...) {
  stats::qnorm(p, mean = theta[[1]], sd = 1 / sqrt(theta[[2]]),
               lower.tail = lower.tail, log.p = log.p)
}

#' @title Gaussian Random Generation in Mean and Precision
#' @name distrib_rng.Gaussian3Distrib
#' @description Delegates to \code{\link[stats]{rnorm}}.
#' @param distrib A \code{Gaussian3Distrib} object.
#' @param n The number of draws.
#' @param theta A list with \code{mu} and \code{tau}.
#' @return A numeric vector.
#' @seealso \code{\link{gaussian3_distrib}}
S7::method(distrib_rng, Gaussian3Distrib) <- function(distrib, n, theta) {
  stats::rnorm(n, mean = theta[[1]], sd = 1 / sqrt(theta[[2]]))
}

#' @title Gaussian Analytical Gradient in Mean and Precision
#' @name distrib_gradient.Gaussian3Distrib
#' @description
#' With \eqn{r = y - \mu},
#' \deqn{\dfrac{\partial\ell}{\partial\mu} = \tau r, \qquad
#'       \dfrac{\partial\ell}{\partial\tau} = \dfrac{1}{2\tau} - \dfrac{r^2}{2}}
#' @param distrib A \code{Gaussian3Distrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list with \code{mu} and \code{tau}.
#' @param scale Either \code{"parameter"} or \code{"link"}; handled by the generic.
#' @param ... Unused.
#' @return A named list of first derivatives.
#' @seealso \code{\link{gaussian3_distrib}}
S7::method(distrib_gradient, Gaussian3Distrib) <- function(distrib, y, theta,
                                                            scale = c("parameter", "link"), ...) {
  gaussian3_gradient_cpp(y, theta[[1]], theta[[2]])
}

#' @title Gaussian Analytical Observed Hessian in Mean and Precision
#' @name distrib_hessian.Gaussian3Distrib
#' @description
#' \deqn{\ell^{(\mu\mu)} = -\tau, \qquad \ell^{(\mu\tau)} = r, \qquad
#'       \ell^{(\tau\tau)} = -\dfrac{1}{2\tau^2}}
#' The mean block is free of the data here, as it is in every parametrization
#' of this family, and so is the precision block.
#' @param distrib A \code{Gaussian3Distrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list with \code{mu} and \code{tau}.
#' @param scale Either \code{"parameter"} or \code{"link"}; handled by the generic.
#' @param ... Unused.
#' @return A named list of second derivatives.
#' @seealso \code{\link{gaussian3_distrib}}
S7::method(distrib_hessian, Gaussian3Distrib) <- function(distrib, y, theta,
                                                           scale = c("parameter", "link"), ...) {
  gaussian3_hessian_cpp(y, theta[[1]], theta[[2]])
}

#' @title Gaussian Analytical Expected Hessian in Mean and Precision
#' @name distrib_expected_hessian.Gaussian3Distrib
#' @description
#' \deqn{\mathbb{E}[\ell^{(\mu\mu)}] = -\tau, \qquad
#'       \mathbb{E}[\ell^{(\mu\tau)}] = 0, \qquad
#'       \mathbb{E}[\ell^{(\tau\tau)}] = -\dfrac{1}{2\tau^2}}
#' Only the mixed entry differs from the observed Hessian, which is what makes
#' Fisher scoring and Newton's method take the same step on the parameter scale
#' here.
#' @param distrib A \code{Gaussian3Distrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list with \code{mu} and \code{tau}.
#' @param scale Either \code{"parameter"} or \code{"link"}; handled by the generic.
#' @param approx Ignored; the expectation is closed form.
#' @param nsim Ignored.
#' @param ... Unused.
#' @return A named list of expected second derivatives.
#' @seealso \code{\link{gaussian3_distrib}}
S7::method(distrib_expected_hessian, Gaussian3Distrib) <- function(distrib, y, theta,
                                                                    scale = c("parameter", "link"),
                                                                    approx = c("bartlett", "integrate", "mc", "opg"),
                                                                    nsim = 10000, ...) {
  gaussian3_expected_hessian_cpp(y, theta[[1]], theta[[2]])
}

#' @title Gaussian Third-Order Derivatives in Mean and Precision
#' @name distrib_deriv3.Gaussian3Distrib
#' @description
#' Closed form. Every third derivative is free of the response, so the observed
#' and the expected ones coincide and \code{expected} changes nothing.
#' @param distrib A \code{Gaussian3Distrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list with \code{mu} and \code{tau}.
#' @param expected Logical; makes no difference here.
#' @param scale Either \code{"parameter"} or \code{"link"}; handled by the generic.
#' @param approx Ignored.
#' @param nsim Ignored.
#' @param ... Unused.
#' @return A named list of third-derivative components.
#' @seealso \code{\link{gaussian3_distrib}}
S7::method(distrib_deriv3, Gaussian3Distrib) <- function(distrib, y, theta, expected = FALSE,
                                                          scale = c("parameter", "link"),
                                                          approx = c("integrate", "bartlett", "mc", "opg"),
                                                          nsim = 10000, ...) {
  gaussian3_deriv3_cpp(y, theta[[1]], theta[[2]])
}

#' @title Gaussian Fourth-Order Derivatives in Mean and Precision
#' @name distrib_deriv4.Gaussian3Distrib
#' @description
#' Closed form, and free of the response: the only non-zero component is
#' \eqn{\ell^{(\tau\tau\tau\tau)} = -3/\tau^4}.
#' @param distrib A \code{Gaussian3Distrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list with \code{mu} and \code{tau}.
#' @param expected Logical; makes no difference here.
#' @param scale Either \code{"parameter"} or \code{"link"}; handled by the generic.
#' @param approx Ignored.
#' @param nsim Ignored.
#' @param ... Unused.
#' @return A named list of fourth-derivative components.
#' @seealso \code{\link{gaussian3_distrib}}
S7::method(distrib_deriv4, Gaussian3Distrib) <- function(distrib, y, theta, expected = FALSE,
                                                          scale = c("parameter", "link"),
                                                          approx = c("integrate", "bartlett", "mc", "opg"),
                                                          nsim = 10000, ...) {
  gaussian3_deriv4_cpp(y, theta[[1]], theta[[2]])
}

#' @title Gaussian Response Derivatives in Mean and Precision
#' @name distrib_grad_y.Gaussian3Distrib
#' @description \eqn{\partial\ell/\partial y = -\tau(y-\mu)}.
#' @param distrib A \code{Gaussian3Distrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list with \code{mu} and \code{tau}.
#' @param ... Unused.
#' @return A numeric vector.
#' @seealso \code{\link{gaussian3_distrib}}
S7::method(distrib_grad_y, Gaussian3Distrib) <- function(distrib, y, theta, ...) {
  -theta[[2]] * (y - theta[[1]])
}

#' @title Gaussian Second Response Derivative in Mean and Precision
#' @name distrib_hess_y.Gaussian3Distrib
#' @description \eqn{\partial^2\ell/\partial y^2 = -\tau}, free of the response.
#' @param distrib A \code{Gaussian3Distrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list with \code{mu} and \code{tau}.
#' @param ... Unused.
#' @return A numeric vector.
#' @seealso \code{\link{gaussian3_distrib}}
S7::method(distrib_hess_y, Gaussian3Distrib) <- function(distrib, y, theta, ...) {
  rep(-theta[[2]], length.out = length(y))
}


#' Gaussian Distribution in Mean and Precision
#'
#' @description
#' Creates a Gaussian distribution object parametrized by its mean and its
#' \strong{precision} \eqn{\tau = 1/\sigma^2}.
#'
#' @details
#' This is the same law as \code{\link{gaussian1_distrib}} in different
#' coordinates, and a separate family for the same reason
#' \code{\link{gaussian2_distrib}} is: the parameter here \emph{is} the
#' precision, and that is what the estimate, the standard error and the
#' interval describe.
#'
#' It is the flattest of the three parametrizations. Every third derivative is
#' free of the response, so the observed and the expected ones coincide, and
#' the only non-zero fourth derivative is
#' \eqn{\ell^{(\tau\tau\tau\tau)} = -3/\tau^4}.
#'
#' The precision is the parametrization a Bayesian conjugate analysis uses,
#' the gamma being conjugate for \eqn{\tau} at known \eqn{\mu}.
#'
#' @section The distribution:
#' \deqn{f(y) = \sqrt{\frac{\tau}{2\pi}}\exp\!\left\{-\frac{\tau(y-\mu)^{2}}{2}\right\}}
#' on \eqn{y \in \mathbb{R}}.
#'
#' \deqn{\mathbb{E}[Y] = \mu, \qquad \operatorname{Var}(Y) = 1/\tau}
#'
#' @param link_mu Link function for \eqn{\mu}. Defaults to the identity.
#' @param link_tau Link function for \eqn{\tau}. Defaults to the log.
#'
#' @return An S7 object of class \code{\link{Gaussian3Distrib}}.
#'
#' @seealso \code{\link{gaussian1_distrib}}, \code{\link{gaussian2_distrib}}
#'
#' @examples
#' d <- gaussian3_distrib()
#' theta <- list(mu = 1, tau = 0.25)
#' distrib_pdf(d, c(0, 1, 2), theta)
#' variance(d, theta)
#'
#' # the same law as gaussian1 with sigma = 2
#' distrib_pdf(gaussian1_distrib(), 0.5, list(mu = 1, sigma = 2))
#'
#' @export
gaussian3_distrib <- function(link_mu = identity_link(),
                              link_tau = log_link()) {
  Gaussian3Distrib(
    distrib_name = "gaussian3",
    dimension = "univariate",
    bounds = c(-Inf, Inf),
    params = c("mu", "tau"),
    params_interpretation = c(mu = "mean", tau = "precision"),
    n_params = 2,
    params_bounds = list(mu = c(-Inf, Inf), tau = c(0, Inf)),
    link_params = list(mu = link_mu, tau = link_tau)
  )
}
