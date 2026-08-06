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
#' The components of any order from one to four, obtained by carrying the
#' shape parametrization's derivatives through the map
#' \eqn{(a, b) = (\mu/\sigma, (1-\mu)/\sigma)}.
#'
#' @param distrib A \code{\link{BetaBinom1Distrib}} object.
#' @param y A numeric vector of counts.
#' @param theta A list containing \code{mu} and \code{sigma}.
#' @param order The derivative order, 1 to 4.
#'
#' @return A named list of component vectors, keyed as
#'   \code{\link{deriv_names}}.
#'
#' @seealso \code{\link{betabinom1_distrib}}, \code{\link{chain_derivatives}}
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


#' @title Beta-Binomial Third and Fourth Derivatives in Mean and Dispersion
#' @name distrib_deriv3.BetaBinom1Distrib
#' @description
#' Closed form at both orders. The shape parametrization carries closed
#' derivatives at every order, and this one is that one at
#' \eqn{a = \mu/\sigma} and \eqn{b = (1-\mu)/\sigma}, so the partition sum of
#' \code{\link{chain_derivatives}} over the map delivers them. Every partial
#' of the map with two or more \eqn{\mu} vanishes, both shapes being linear in
#' \eqn{\mu} at fixed \eqn{\sigma}.
#' @param distrib A \code{BetaBinom1Distrib} object.
#' @param y A numeric vector of counts.
#' @param theta A list containing \code{mu} and \code{sigma}.
#' @param expected Logical; if \code{TRUE}, the expected derivatives.
#' @param scale Either \code{"parameter"} or \code{"link"}; handled by the generic.
#' @param approx The approximation used when \code{expected} is \code{TRUE}.
#' @param nsim Monte Carlo draws when \code{approx = "mc"}.
#' @param ... Unused.
#' @return A named list of third-derivative components.
#' @seealso \code{\link{betabinom1_distrib}}
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

#' @rdname distrib_deriv3.BetaBinom1Distrib
#' @name distrib_deriv4.BetaBinom1Distrib
#' @return A named list of fourth-derivative components.
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
