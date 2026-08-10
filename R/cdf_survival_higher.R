#' @include cdf_mapped_higher.R
NULL

# The families whose survival function is an exponential of something
# elementary.
#
# With S = 1 - F written as exp(L), every derivative of F in the parameters is
# d^I F = -S B_I(L), the complete Bell polynomial in L's own partial
# derivatives -- the same identity the wrappers use, read on the survival
# function instead of on the density. So a family needs only to say what L is
# and what its partials are, and all four orders follow at once. Nothing here
# is transcribed from an expansion: `bell_f_ratio()` runs the partition sum.

#' CDF Derivatives From an Exponential Survival Function
#'
#' @description
#' Returns \eqn{\partial^{I}F} for every component of the requested order,
#' given \eqn{L = \log(1-F)} and a function evaluating its partial derivatives.
#'
#' @details
#' \eqn{S = e^{L}} gives \eqn{\partial^{I}S = S\,B_{I}}, with \eqn{B_{I}} the
#' complete Bell polynomial in the partials of \eqn{L}, and \eqn{F = 1 - S}
#' turns that into \eqn{\partial^{I}F = -S\,B_{I}}. The survival function is
#' evaluated as \code{exp(L)} rather than as \code{1 - F}, which keeps the far
#' tail from cancelling.
#'
#' @param distrib An object inheriting from class \code{"distrib"}.
#' @param q A numeric vector of quantiles.
#' @param theta A named list of parameters.
#' @param order The derivative order, 1 to 4.
#' @param Lval The value of \eqn{L} at \code{q}.
#' @param Lderiv A function of a character vector of parameter names returning
#'   the corresponding partial derivative of \eqn{L}.
#'
#' @return A named list of derivative components of \eqn{F}.
#'
#' @seealso \code{\link{bell_f_ratio}}
#' @keywords internal
surv_cdf_deriv_k <- function(distrib, q, theta, order, Lval, Lderiv) {
  params <- distrib@params
  # below the support F is identically zero and so is every derivative; L is
  # still finite there and would otherwise produce a survival above one
  inside <- q > distrib@bounds[1L]
  S <- exp(Lval)
  idx <- deriv_indices(params, order)
  out <- lapply(idx, function(I) -S * bell_f_ratio(params[I], Lderiv) * inside)
  stats::setNames(out, deriv_names(params, order))
}

#' Register the Four CDF Derivative Orders of an Exponential Survival Family
#'
#' @description
#' Turns a function returning \eqn{L} and its partial-derivative evaluator into
#' the four methods.
#'
#' @param cls The S7 class.
#' @param pieces A function of \code{(distrib, q, theta)} returning a list with
#'   \code{Lval} and \code{Lderiv}.
#'
#' @return Invisibly \code{NULL}; called for the registration.
#'
#' @keywords internal
register_surv_cdf <- function(cls, pieces) {
  make <- function(o) {
    force(o)
    function(distrib, q, theta, lower.tail = TRUE, log = TRUE, ...) {
      p <- pieces(distrib, q, theta)
      params <- distrib@params
      if (!lower.tail && log) {
        # log S is L, so the components ARE its partial derivatives. Taking
        # them the general way instead would divide by S = 1 - F read off the
        # natural scale, which is exactly zero past q/mu = 37 for an
        # exponential while L is finite to the end of the range.
        inside <- q > distrib@bounds[1L]
        out <- lapply(deriv_indices(params, o),
                      function(I) rep_len(p$Lderiv(params[I]), length(q)) * inside)
        return(stats::setNames(out, deriv_names(params, o)))
      }
      tabs <- lapply(seq_len(o), function(k)
        surv_cdf_deriv_k(distrib, q, theta, k, p$Lval, p$Lderiv))
      if (!lower.tail) {
        # d^I S = -d^I F, with S itself taken as exp(L)
        return(lapply(tabs[[o]], function(v) -v))
      }
      cdf_scale_k(distrib, distrib_cdf(distrib, q, theta), tabs, o,
                  lower.tail, log)
    }
  }
  S7::method(distrib_grad_cdf, cls) <- make(1L)
  S7::method(distrib_hess_cdf, cls) <- make(2L)
  S7::method(distrib_deriv3_cdf, cls) <- make(3L)
  S7::method(distrib_deriv4_cdf, cls) <- make(4L)
  invisible(NULL)
}

#' @title Exponential Log-CDF Derivatives
#' @name distrib_grad_cdf.ExponentialDistrib
#' @description
#' Closed form at every order from the survival function
#' \eqn{S = \exp(-q/\mu)}, whose logarithm has the partial derivatives
#' \eqn{\partial^{j}L/\partial\mu^{j} = -q(-1)^{j}j!/\mu^{j+1}}.
#' @param distrib An \code{ExponentialDistrib} object.
#' @param q A numeric vector of quantiles.
#' @param theta A list containing \code{mu}.
#' @param lower.tail Logical; if \code{TRUE} (default), the lower tail.
#' @param log Logical; if \code{TRUE} (default), derivatives of the log probability.
#' @param ... Unused.
#' @return A named list, one vector per component.
#' @seealso \code{\link{exponential_distrib}}
#' @keywords internal
register_surv_cdf(ExponentialDistrib, function(distrib, q, theta) {
  mu <- theta[[1]]
  list(
    Lval = -q / mu,
    Lderiv = function(block) {
      j <- length(block)
      -q * (-1)^j * factorial(j) / mu^(j + 1)
    }
  )
})

#' @title Weibull Log-CDF Derivatives
#' @name distrib_grad_cdf.Weibull1Distrib
#' @description
#' Closed form at every order from the survival function
#' \eqn{S = \exp\{-(q/\mu)^{\sigma}\}}. Writing \eqn{h = \sigma(\log q -
#' \log\mu)} the exponent is \eqn{L = -e^{h}}, so its partial derivatives are
#' \eqn{-e^{h}} times the complete Bell polynomial in the partials of \eqn{h},
#' and those are elementary: \eqn{\partial^{j}h/\partial\mu^{j} =
#' \sigma(-1)^{j}(j-1)!/\mu^{j}}, the same without the \eqn{\sigma} when one
#' index names the shape, and zero when two do.
#' @param distrib A \code{Weibull1Distrib} object.
#' @param q A numeric vector of quantiles.
#' @param theta A list containing \code{mu} and \code{sigma}.
#' @param lower.tail Logical; if \code{TRUE} (default), the lower tail.
#' @param log Logical; if \code{TRUE} (default), derivatives of the log probability.
#' @param ... Unused.
#' @return A named list, one vector per component.
#' @seealso \code{\link{weibull1_distrib}}
#' @keywords internal
register_surv_cdf(Weibull1Distrib, function(distrib, q, theta) {
  mu <- theta[[1]]
  sigma <- theta[[2]]
  nm <- distrib@params
  h <- sigma * (log(q) - log(mu))
  # h is linear in the shape and a logarithm in the scale, so a block naming
  # the shape twice is zero and one naming it once drops the factor sigma
  hderiv <- function(block) {
    ns <- sum(block == nm[2L])
    j <- sum(block == nm[1L])
    if (ns >= 2L) return(0 * q)
    base <- if (j == 0L) rep_len(1, length(q)) else
      (-1)^j * factorial(j - 1L) / mu^j
    if (ns == 1L) {
      if (j == 0L) log(q) - log(mu) else base
    } else {
      sigma * base
    }
  }
  eh <- exp(h)
  list(
    Lval = -eh,
    Lderiv = function(block) -eh * bell_f_ratio(block, hderiv)
  )
})
