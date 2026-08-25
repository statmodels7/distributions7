#' @include betabinom1_distrib.R betabinom2_distrib.R reparametrize.R
NULL

# The beta-binomial in its mean and dispersion is the shape parametrization at
# a = mu/sigma and b = (1 - mu)/sigma, and the shape parametrization already
# carries closed derivatives at every order. So orders three and four here are
# the partition sum of chain_derivatives() over that map, with the map's own
# partials written out below -- the same construction reparametrize() runs,
# used by a family that is written by hand for its first two orders.
#
# Every partial of the map with two or more mu is exactly zero, both shapes
# being linear in mu at fixed sigma, so the table is sparse.

#' @rdname reparam_map_derivs
#' @keywords internal
md_betabinom1 <- function(psi) {
  m <- psi[[1]]
  s <- psi[[2]]
  one <- rep_len(1, max(length(m), length(s)))
  q <- 1 - m
  list(
    list("1" = one / s, "2" = -m / s^2,
         "1,2" = -one / s^2, "2,2" = 2 * m / s^3,
         "1,2,2" = 2 * one / s^3, "2,2,2" = -6 * m / s^4,
         "1,2,2,2" = -6 * one / s^4, "2,2,2,2" = 24 * m / s^5),
    list("1" = -one / s, "2" = -q / s^2,
         "1,2" = one / s^2, "2,2" = 2 * q / s^3,
         "1,2,2" = -2 * one / s^3, "2,2,2" = -6 * q / s^4,
         "1,2,2,2" = 6 * one / s^4, "2,2,2,2" = 24 * q / s^5)
  )
}

#' Derivative Components of the Beta-Binomial in Mean and Dispersion
#'
#' @description
#' Returns the components of a derivative of order one to four of the
#' beta-binomial log-mass with respect to \eqn{\mu} and \eqn{\sigma}, by
#' carrying the shape parametrization's own derivatives through the map
#' \eqn{(\alpha, \beta) = (\mu/\sigma, (1-\mu)/\sigma)}.
#'
#' @details
#' The shape parametrization [betabinom2_distrib()] carries closed derivatives
#' at every order, each a difference of polygammas, and this parametrization is
#' that one composed with a map. [chain_derivatives()] runs the Faa di Bruno
#' partition sum over the map, with the map's own partials supplied by
#' [md_betabinom1()] as a keyed table; a key absent from the table is an exact
#' zero. Both shapes are linear in \eqn{\mu} at fixed \eqn{\sigma}, so every
#' partial carrying two or more \eqn{\mu} vanishes and the sum is short.
#'
#' The construction is the one [reparametrize()] runs, used here by a family
#' whose first two orders are written out by hand in a compiled kernel. A
#' fresh `BetaBinom2Distrib` object is built on each call, which costs one S7
#' construction per evaluation.
#'
#' @param distrib A [BetaBinom1Distrib()] object, read for its `size` and
#'   `params`.
#' @param y A numeric vector of counts in \eqn{\{0, \dots, n\}}.
#' @param theta A named list with components `mu` and `sigma`, each a numeric
#'   vector of length 1 or of the length of `y`. `mu` must lie in
#'   \eqn{(0, 1)} and `sigma` be strictly positive.
#' @param order The derivative order, an integer from 1 to 4.
#'
#' @return A named list of component vectors, one per distinct multi-index of
#'   the given order and keyed as [deriv_names()] keys them, so four
#'   components at order 3 and five at order 4. Each has the recycled length of
#'   the inputs.
#'
#' @seealso [distrib_deriv3.BetaBinom1Distrib()] and
#'   [distrib_deriv4.BetaBinom1Distrib()], which call this;
#'   [chain_derivatives()] for the partition sum;
#'   [betabinom1_distrib()] for the family.
#' @keywords internal
betabinom1_components <- function(distrib, y, theta, order) {
  m <- theta[[1]]
  s <- theta[[2]]
  chain_derivatives(
    parent = betabinom2_distrib(size = distrib@size),
    y = y,
    th_par = list(alpha = m / s, beta = (1 - m) / s),
    maps = md_betabinom1(theta[1:2]),
    new_params = distrib@params,
    order = order
  )
}


#' @title Beta-Binomial Third-Order Derivatives in Mean and Dispersion
#' @name distrib_deriv3.BetaBinom1Distrib
#' @description
#' Computes the four distinct third derivatives of the beta-binomial log-mass
#' in \eqn{\mu} and \eqn{\sigma}, **in closed form**. The shape parametrization
#' carries closed derivatives at every order, each a difference of polygammas,
#' and this parametrization is that one at \eqn{\alpha = \mu/\sigma} and
#' \eqn{\beta = (1-\mu)/\sigma}, so the Faa di Bruno partition sum of
#' [chain_derivatives()] over the map delivers them. Both shapes are linear in
#' \eqn{\mu} at fixed \eqn{\sigma}, so every partial of the map carrying two or
#' more \eqn{\mu} vanishes and the sum is short.
#'
#' With `expected = TRUE` the method calls [expected_derivative()] instead.
#' That is the one place on this page where `approx` and `nsim` are read; on
#' the observed branch both are ignored.
#'
#' @param distrib A `BetaBinom1Distrib` object, from [betabinom1_distrib()].
#' @param y A numeric vector of counts in \eqn{\{0, \dots, n\}}. With
#'   `expected = TRUE` only its length is read.
#' @param theta A named list with components `mu` and `sigma`, each a numeric
#'   vector of length 1 or of the length of `y`. A component of length 1 is
#'   recycled. `mu` must lie in \eqn{(0, 1)} and `sigma` be strictly positive.
#' @param expected Logical of length 1. When `TRUE` the expectation under the
#'   model is returned in place of the value at the data. Defaults to `FALSE`.
#' @param scale One of `"parameter"` (the default) or `"link"`, matched by
#'   [base::match.arg()]. Read by the generic, not by this method.
#' @param approx One of `"integrate"` (the default here), `"bartlett"`, `"mc"`
#'   or `"opg"`, the strategy [expected_derivative()] uses. Read only when
#'   `expected = TRUE`.
#' @param nsim A single positive integer, the sample size when
#'   `approx = "mc"`. Read only when `expected = TRUE`. Defaults to `10000`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A named list of four numeric vectors, `mu_mu_mu`, `mu_mu_sigma`,
#'   `mu_sigma_sigma` and `sigma_sigma_sigma`, each of length
#'   `max(length(y), length(mu), length(sigma))`. The four name the distinct
#'   entries of a symmetric third-order array over two parameters.
#'
#' @section Notation:
#' \eqn{\ell} is the log-mass of one observation, \eqn{\mu \in (0,1)} the mean
#' proportion, \eqn{\sigma > 0} the dispersion and \eqn{n} the trial count.
#' \eqn{\alpha} and \eqn{\beta} are the two beta shapes the family is written
#' in internally.
#'
#' @seealso [distrib_hessian.BetaBinom1Distrib()] for the order below,
#'   [distrib_deriv4.BetaBinom1Distrib()] for the order above,
#'   [betabinom1_components()] for the assembly, and [distrib_deriv3()] for the
#'   generic.
#'
#' @examples
#' d <- betabinom1_distrib(size = 10)
#' th <- list(mu = 0.3, sigma = 0.5)
#' d3 <- distrib_deriv3(d, 0:10, th)
#' names(d3)
#'
#' # A central difference of the Hessian reproduces the pure-mu component,
#' # which is what says the partition sum over the map is right.
#' eps <- 1e-5
#' up <- distrib_hessian(d, 0:10, list(mu = 0.3 + eps, sigma = 0.5))$mu_mu
#' dn <- distrib_hessian(d, 0:10, list(mu = 0.3 - eps, sigma = 0.5))$mu_mu
#' all.equal((up - dn) / (2 * eps), d3$mu_mu_mu, tolerance = 1e-6)
#'
#' # The expected branch is a numerical expectation of the observed one, and
#' # the mass-weighted sum over the support reaches it.
#' w <- distrib_pdf(d, 0:10, th)
#' rbind(expected = vapply(distrib_deriv3(d, 0, th, expected = TRUE),
#'                         function(v) v[1], numeric(1)),
#'       summed = vapply(d3, function(v) sum(w * v), numeric(1)))
S7::method(distrib_deriv3, BetaBinom1Distrib) <- function(distrib, y, theta,
                                                           expected = FALSE,
                                                           scale = c("parameter", "link"),
                                                           approx = c("integrate", "bartlett", "mc", "opg"),
                                                           nsim = 10000, ...) {
  if (expected) {
    return(expected_derivative(distrib, y, theta, order = 3L,
                               approx = match.arg(approx), nsim = nsim))
  }
  betabinom1_components(distrib, y, theta, 3L)
}

#' @title Beta-Binomial Fourth-Order Derivatives in Mean and Dispersion
#' @name distrib_deriv4.BetaBinom1Distrib
#' @description
#' Computes the five distinct fourth derivatives of the beta-binomial log-mass
#' in \eqn{\mu} and \eqn{\sigma}, **in closed form**, by the construction
#' [distrib_deriv3.BetaBinom1Distrib()] describes carried one order further:
#' the shape parametrization's fourth derivatives, each a difference of
#' \eqn{\psi^{(3)}}, carried through the map
#' \eqn{(\alpha, \beta) = (\mu/\sigma, (1-\mu)/\sigma)} by the partition sum of
#' [chain_derivatives()].
#'
#' With `expected = TRUE` the method calls [expected_derivative()] instead.
#' That is the one place on this page where `approx` and `nsim` are read.
#'
#' @param distrib A `BetaBinom1Distrib` object, from [betabinom1_distrib()].
#' @param y A numeric vector of counts in \eqn{\{0, \dots, n\}}. With
#'   `expected = TRUE` only its length is read.
#' @param theta A named list with components `mu` and `sigma`, each a numeric
#'   vector of length 1 or of the length of `y`. A component of length 1 is
#'   recycled. `mu` must lie in \eqn{(0, 1)} and `sigma` be strictly positive.
#' @param expected Logical of length 1. When `TRUE` the expectation under the
#'   model is returned in place of the value at the data. Defaults to `FALSE`.
#' @param scale One of `"parameter"` (the default) or `"link"`, matched by
#'   [base::match.arg()]. Read by the generic, not by this method.
#' @param approx One of `"integrate"` (the default here), `"bartlett"`, `"mc"`
#'   or `"opg"`, the strategy [expected_derivative()] uses. Read only when
#'   `expected = TRUE`.
#' @param nsim A single positive integer, the sample size when
#'   `approx = "mc"`. Read only when `expected = TRUE`. Defaults to `10000`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A named list of five numeric vectors, `mu_mu_mu_mu`,
#'   `mu_mu_mu_sigma`, `mu_mu_sigma_sigma`, `mu_sigma_sigma_sigma` and
#'   `sigma_sigma_sigma_sigma`, each of length
#'   `max(length(y), length(mu), length(sigma))`.
#'
#' @section Notation:
#' \eqn{\ell} is the log-mass of one observation, \eqn{\mu \in (0,1)} the mean
#' proportion, \eqn{\sigma > 0} the dispersion and \eqn{n} the trial count.
#' \eqn{\alpha} and \eqn{\beta} are the two beta shapes the family is written
#' in internally.
#'
#' @seealso [distrib_deriv3.BetaBinom1Distrib()] for the order below,
#'   [distrib_hessian.BetaBinom1Distrib()] for the second order,
#'   [betabinom1_components()] for the assembly, and [distrib_deriv4()] for the
#'   generic.
#'
#' @examples
#' d <- betabinom1_distrib(size = 10)
#' th <- list(mu = 0.3, sigma = 0.5)
#' d4 <- distrib_deriv4(d, 0:10, th)
#' names(d4)
#'
#' # A central difference of the third order reproduces a mixed component.
#' eps <- 1e-5
#' up <- distrib_deriv3(d, 0:10, list(mu = 0.3, sigma = 0.5 + eps))$mu_mu_sigma
#' dn <- distrib_deriv3(d, 0:10, list(mu = 0.3, sigma = 0.5 - eps))$mu_mu_sigma
#' all.equal((up - dn) / (2 * eps), d4$mu_mu_sigma_sigma, tolerance = 1e-5)
#'
#' # The expected branch is a numerical expectation, and the mass-weighted
#' # sum over the support reaches it.
#' w <- distrib_pdf(d, 0:10, th)
#' rbind(expected = vapply(distrib_deriv4(d, 0, th, expected = TRUE),
#'                         function(v) v[1], numeric(1)),
#'       summed = vapply(d4, function(v) sum(w * v), numeric(1)))
S7::method(distrib_deriv4, BetaBinom1Distrib) <- function(distrib, y, theta,
                                                           expected = FALSE,
                                                           scale = c("parameter", "link"),
                                                           approx = c("integrate", "bartlett", "mc", "opg"),
                                                           nsim = 10000, ...) {
  if (expected) {
    return(expected_derivative(distrib, y, theta, order = 4L,
                               approx = match.arg(approx), nsim = nsim))
  }
  betabinom1_components(distrib, y, theta, 4L)
}
