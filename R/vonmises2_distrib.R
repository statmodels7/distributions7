#' @include distrib.R generics.R vonmises1_distrib.R moments.R
NULL

# The von Mises in its mean resultant length. The concentration and the
# resultant length are related by rho = A(kappa) = I_1(kappa)/I_0(kappa), a
# strictly increasing bijection from (0, Inf) onto (0, 1) whose inverse has no
# closed form. numericals7::bessel_i_ratio_inverse() obtains it by root
# finding and differentiates it by the inverse function rule, with A' to
# A'''' from the Bessel recurrences rather than four more evaluations.

#' @title von Mises Distribution Class, Mean Resultant Length
#' @name VonMises2Distrib
#'
#' @description
#' The S7 class of the von Mises family written in its mean direction
#' \eqn{\mu} and its **mean resultant length** \eqn{\rho \in (0, 1)}, the
#' quantity circular statistics reports and one minus the circular variance.
#' It inherits from `continuous_distrib`; the nine methods listed below are
#' registered on it directly.
#'
#' The concentration of [vonmises1_distrib()] is recovered as
#' \eqn{\kappa = A^{-1}(\rho)} with \eqn{A(\kappa) = I_1(\kappa)/I_0(\kappa)},
#' a strictly increasing bijection from \eqn{(0, \infty)} onto \eqn{(0, 1)}
#' whose inverse has no closed form. That is why the family is written out
#' here instead of through [reparametrize()].
#'
#' Build one with [vonmises2_distrib()], which supplies the two link functions
#' and fills the properties in. This page documents the raw S7 constructor,
#' which takes the parent's properties and validates none of the relationships
#' between them.
#'
#' @inheritParams distrib
#'
#' @return An S7 object of class `VonMises2Distrib`, inheriting from
#'   `continuous_distrib` and from `distrib`. Its properties are the parent's:
#'   `distrib_name`, `dimension`, `bounds`, `params`, `params_interpretation`,
#'   `n_params`, `params_bounds`, `link_params` and `params_smooth`. For an
#'   object built by [vonmises2_distrib()] they hold `"von mises2"`,
#'   `"univariate"`, `c(-pi, pi)`, `c("mu", "rho")`, the interpretations
#'   `c(mu = "mean direction", rho = "mean resultant length")`, `2`, and the
#'   domains \eqn{(-\pi, \pi)} and \eqn{(0, 1)}.
#'
#' @seealso [vonmises2_distrib()] to build one;
#'   [vonmises1_distrib()] for the same law in the concentration;
#'   [numericals7::bessel_i_ratio_inverse()] for the map;
#'   [distrib_expected_hessian.VonMises2Distrib()] for the information.
#'
#' @section Methods:
#' Registered on this class:
#'   [`distrib_cdf()`][distrib_cdf.VonMises2Distrib],
#'   [`distrib_deriv3()`][distrib_deriv3.VonMises2Distrib],
#'   [`distrib_deriv4()`][distrib_deriv4.VonMises2Distrib],
#'   [`distrib_expected_hessian()`][distrib_expected_hessian.VonMises2Distrib],
#'   [`distrib_gradient()`][distrib_gradient.VonMises2Distrib],
#'   [`distrib_hessian()`][distrib_hessian.VonMises2Distrib],
#'   [`distrib_pdf()`][distrib_pdf.VonMises2Distrib],
#'   [`distrib_rng()`][distrib_rng.VonMises2Distrib],
#'   [`mean()`][mean.VonMises2Distrib].
#'
#' The **quantile** and the response derivatives come from
#' [continuous_distrib()], the quantile by root finding on this class's own
#' distribution function. The remaining moments come from [variance()] and its
#' siblings, numerically, and are the ordinary moments of \eqn{Y} as a number.
#'
#' @examples
#' d <- vonmises2_distrib()
#' S7::S7_inherits(d, continuous_distrib)
#'
#' # The resultant length is bounded, which is what makes it readable.
#' d@params
#' d@params_bounds$rho
#' d@params_interpretation
#'
#' # The direction rides a bounded link and the resultant length a logit.
#' vapply(d@link_params, function(l) l@link_name, character(1))
VonMises2Distrib <- S7::new_class("VonMises2Distrib", parent = continuous_distrib)

#' The Pieces a von Mises Derivative in the Resultant Length Needs
#'
#' @description
#' Evaluates the concentration \eqn{\kappa = A^{-1}(\rho)}, its four
#' derivatives in \eqn{\rho}, and the derivatives of
#' \eqn{A(\kappa) = I_1(\kappa)/I_0(\kappa)} at that concentration, once per
#' call, so that a density or a derivative method shares them.
#'
#' @details
#' The inverse of \eqn{A} has no closed form.
#' [numericals7::bessel_i_ratio_inverse()] obtains \eqn{\kappa} by root finding
#' on \eqn{\log\kappa} and differentiates it by the inverse function rule, so
#' \eqn{\mathrm{d}\kappa/\mathrm{d}\rho = 1/A'(\kappa)} and the higher
#' derivatives follow from \eqn{A'} to \eqn{A^{(4)}}. Those come from the
#' Bessel recurrences and from the same two evaluations \eqn{A} already needs,
#' so no further Bessel call is made at any order.
#'
#' The map is very steep near \eqn{\rho = 1}: measured,
#' \eqn{\mathrm{d}\kappa/\mathrm{d}\rho} is 2.00 at \eqn{\rho = 0.01}, 3.14 at
#' 0.5, 199 at 0.95 and \eqn{5.0\times10^{5}} at 0.999.
#'
#' @param theta A named list with components `mu` and `rho`, each a numeric
#'   vector of length 1 or of a common length. `rho` must lie in \eqn{(0, 1)};
#'   only it is read here, `mu` not entering the map.
#'
#' @return A named list with `kappa`, the concentration; `kd`, the result of
#'   [numericals7::bessel_i_ratio_inverse()], carrying `kappa` and its
#'   derivatives `d1` to `d4` in \eqn{\rho}; and `ad`, the result of
#'   [numericals7::bessel_i_ratio_derivs()] at that concentration, carrying `A`
#'   and its derivatives `d1` to `d4` in \eqn{\kappa}.
#'
#' @seealso [distrib_gradient.VonMises2Distrib()] for the first consumer,
#'   [numericals7::bessel_i_ratio_inverse()] for the root finding, and
#'   [vonmises2_distrib()] for the family.
#'
#' @keywords internal
vm2_parts <- function(theta) {
  kd <- numericals7::bessel_i_ratio_inverse(theta[[2]])
  list(kappa = kd$kappa, kd = kd, ad = numericals7::bessel_i_ratio_derivs(kd$kappa))
}

# --- S7 METHODS IMPLEMENTATION ---

#' @title von Mises Density in the Resultant Length
#' @name distrib_pdf.VonMises2Distrib
#' @description
#' Computes the von Mises density
#' \deqn{f(y; \mu, \rho) = \dfrac{e^{\kappa\cos(y-\mu)}}{2\pi I_0(\kappa)},
#'       \qquad \kappa = A^{-1}(\rho),}
#' by inverting the map once and calling
#' [distrib_pdf.VonMises1Distrib()] at the implied concentration. The two
#' parametrizations are the same law, so the value is identical to that
#' family's at \eqn{\kappa}.
#'
#' @param distrib A `VonMises2Distrib` object, from [vonmises2_distrib()].
#' @param y A numeric vector of angles. A value outside \eqn{[-\pi, \pi)} is
#'   off the support and gives a density of 0.
#' @param theta A named list with components `mu` and `rho`, each a numeric
#'   vector of length 1 or of the length of `y`. `mu` must lie in
#'   \eqn{(-\pi, \pi)} and `rho` in \eqn{(0, 1)}.
#' @param log Logical of length 1. When `TRUE` the log-density is returned.
#'   Defaults to `FALSE`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#' @param threads A single positive integer, passed on to
#'   [numericals7::log_bessel_i()] through the concentration parametrization.
#'   Defaults to `1L`.
#'
#' @return A numeric vector of densities, of length
#'   `max(length(y), length(mu), length(rho))`, one value per observation.
#'
#' @section Notation:
#' \eqn{\mu} is the mean direction, \eqn{\rho \in (0,1)} the mean resultant
#' length, \eqn{\kappa} the concentration, \eqn{I_m} the modified Bessel
#' function of the first kind of order \eqn{m}, and
#' \eqn{A(\kappa) = I_1(\kappa)/I_0(\kappa)}.
#'
#' @seealso [distrib_pdf.VonMises1Distrib()], which this calls;
#'   [distrib_cdf.VonMises2Distrib()] for the distribution function; and
#'   [distrib_pdf()] for the generic.
#'
#' @examples
#' d2 <- vonmises2_distrib()
#' th <- list(mu = 0.5, rho = 0.7)
#' y <- c(-1, 0, 0.5, 2)
#' distrib_pdf(d2, y, th)
#'
#' # It is the same law as the concentration parametrization, at the
#' # concentration this resultant length implies.
#' k <- numericals7::bessel_i_ratio_inverse(0.7)$kappa
#' k
#' all.equal(distrib_pdf(d2, y, th),
#'           distrib_pdf(vonmises1_distrib(), y, list(mu = 0.5, kappa = k)))
#'
#' # And that concentration maps back to the resultant length given.
#' numericals7::bessel_i_ratio(k)
S7::method(distrib_pdf, VonMises2Distrib) <- function(distrib, y, theta,
                                                     log = FALSE, ...,
                                                     threads = 1L) {
  distrib_pdf(vonmises1_distrib(), y,
              list(mu = theta[[1]], kappa = vm2_parts(theta)$kappa), log = log,
              threads = threads)
}

#' @title von Mises Random Generation in the Resultant Length
#' @name distrib_rng.VonMises2Distrib
#' @description
#' Draws `n` independent angles by inverting the map once and calling
#' [distrib_rng.VonMises1Distrib()] at the implied concentration
#' \eqn{\kappa = A^{-1}(\rho)}. The generator is the rejection algorithm of
#' Best and Fisher (1979), which involves no Bessel function; the only Bessel
#' work here is the single inversion of \eqn{A}.
#'
#' @param distrib A `VonMises2Distrib` object, from [vonmises2_distrib()].
#' @param n A single positive integer, the number of draws.
#' @param theta A named list with components `mu` and `rho`, each a numeric
#'   vector of length 1. `mu` must lie in \eqn{(-\pi, \pi)} and `rho` in
#'   \eqn{(0, 1)}.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of `n` angles in \eqn{[-\pi, \pi)}.
#'
#' @seealso [distrib_rng.VonMises1Distrib()], which this calls;
#'   [fit_distrib()] to estimate the parameters back from a sample; and
#'   [distrib_rng()] for the generic.
#'
#' @examples
#' d2 <- vonmises2_distrib()
#' th <- list(mu = 0.5, rho = 0.7)
#' set.seed(1)
#' z <- distrib_rng(d2, 3e5, th)
#'
#' # This parametrization is the one a sample reads back directly: the mean
#' # resultant length of the draws is rho, and the circular mean is mu.
#' c(resultant = sqrt(mean(cos(z))^2 + mean(sin(z))^2), rho = 0.7)
#' c(circular_mean = atan2(mean(sin(z)), mean(cos(z))), mu = 0.5)
S7::method(distrib_rng, VonMises2Distrib) <- function(distrib, n, theta, ...) {
  distrib_rng(vonmises1_distrib(), n,
              list(mu = theta[[1]], kappa = vm2_parts(theta)$kappa))
}

#' @title von Mises Score in the Resultant Length
#' @name distrib_gradient.VonMises2Distrib
#' @description
#' Computes the first derivatives of the log-density with respect to the mean
#' direction \eqn{\mu} and the mean resultant length \eqn{\rho}, one value per
#' observation, in closed form:
#' \deqn{\dfrac{\partial\ell}{\partial\mu} = \kappa\sin(y-\mu), \qquad
#'       \dfrac{\partial\ell}{\partial\rho}
#'         = \left\{\cos(y-\mu) - A(\kappa)\right\}\kappa'(\rho),}
#' with \eqn{\kappa = A^{-1}(\rho)} and
#' \eqn{\kappa'(\rho) = 1/A'(\kappa)} from the inverse function rule.
#'
#' The map touches the **second parameter only**, so the chain rule is the
#' one-variable one: the direction's component is unchanged from the
#' concentration parametrization, and the second is that family's multiplied by
#' a single factor. No multivariate expansion and no cancellation are involved
#' at any order.
#'
#' With `scale = "link"` the generic applies the chain rule for the links the
#' family carries. This method always returns the parameter scale.
#'
#' @param distrib A `VonMises2Distrib` object, from [vonmises2_distrib()].
#' @param y A numeric vector of angles in \eqn{[-\pi, \pi)}.
#' @param theta A named list with components `mu` and `rho`, each a numeric
#'   vector of length 1 or of the length of `y`. `mu` must lie in
#'   \eqn{(-\pi, \pi)} and `rho` in \eqn{(0, 1)}.
#' @param scale One of `"parameter"` (the default) or `"link"`, matched by
#'   [base::match.arg()]. Read by the generic, not by this method.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A named list of two numeric vectors, `mu` and `rho`, each of length
#'   `max(length(y), length(mu), length(rho))`.
#'
#' @section Notation:
#' \eqn{\ell} is the log-density of one observation, \eqn{\mu} the mean
#' direction, \eqn{\rho \in (0,1)} the mean resultant length, \eqn{\kappa} the
#' concentration and \eqn{A(\kappa) = I_1(\kappa)/I_0(\kappa)}, so that
#' \eqn{\rho = A(\kappa)}.
#'
#' @seealso [distrib_hessian.VonMises2Distrib()] for the second derivatives,
#'   [distrib_expected_hessian.VonMises2Distrib()] for their expectation,
#'   [distrib_gradient.VonMises1Distrib()] for the same quantity in the
#'   concentration, and [distrib_gradient()] for the generic.
#'
#' @examples
#' d2 <- vonmises2_distrib()
#' y <- c(-1, 0, 0.5, 2)
#' th <- list(mu = 0.5, rho = 0.7)
#' g <- distrib_gradient(d2, y, th)
#'
#' # numDeriv on the summed log-density reproduces the summed score.
#' fn <- function(v)
#'   sum(distrib_pdf(d2, y, list(mu = v[1], rho = v[2]), log = TRUE))
#' rbind(numeric = numDeriv::grad(fn, c(0.5, 0.7)),
#'       closed = vapply(g, sum, numeric(1)))
#'
#' # The direction component is unchanged from the concentration
#' # parametrization, the map touching the second parameter only.
#' k <- numericals7::bessel_i_ratio_inverse(0.7)$kappa
#' all.equal(g$mu,
#'           distrib_gradient(vonmises1_distrib(), y,
#'                            list(mu = 0.5, kappa = k))$mu)
S7::method(distrib_gradient, VonMises2Distrib) <- function(distrib, y, theta,
                                                           scale = c("parameter", "link"), ...) {
  p <- vm2_parts(theta)
  d <- y - theta[[1]]
  list(mu = p$kappa * sin(d),
       rho = (cos(d) - p$ad$A) * p$kd$d1)
}

#' @title von Mises Observed Hessian in the Resultant Length
#' @name distrib_hessian.VonMises2Distrib
#' @description
#' Computes the three distinct second derivatives of the log-density in
#' \eqn{\mu} and \eqn{\rho}, one value per observation, in closed form. The
#' concentration parametrization's second derivatives are carried through the
#' one-variable chain rule,
#' \deqn{\ell^{(\rho\rho)} = \ell^{(\kappa\kappa)}\{\kappa'(\rho)\}^2
#'                          + \ell^{(\kappa)}\kappa''(\rho), \qquad
#'       \ell^{(\mu\rho)} = \ell^{(\mu\kappa)}\kappa'(\rho),}
#' with \eqn{\ell^{(\kappa\kappa)} = -A'(\kappa)},
#' \eqn{\ell^{(\mu\kappa)} = \sin(y-\mu)} and
#' \eqn{\ell^{(\mu\mu)} = -\kappa\cos(y-\mu)} unchanged.
#'
#' Unlike in the concentration parametrization, the pure second derivative is
#' **not** free of the data: the term in \eqn{\kappa''} carries
#' \eqn{\cos(y-\mu)}, which the map's curvature brings in.
#'
#' @param distrib A `VonMises2Distrib` object, from [vonmises2_distrib()].
#' @param y A numeric vector of angles in \eqn{[-\pi, \pi)}.
#' @param theta A named list with components `mu` and `rho`, each a numeric
#'   vector of length 1 or of the length of `y`. `mu` must lie in
#'   \eqn{(-\pi, \pi)} and `rho` in \eqn{(0, 1)}.
#' @param scale One of `"parameter"` (the default) or `"link"`, matched by
#'   [base::match.arg()]. Read by the generic, not by this method.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A named list of three numeric vectors, `mu_mu`, `rho_rho` and
#'   `mu_rho`, each of length `max(length(y), length(mu), length(rho))`.
#'
#' @section Notation:
#' \eqn{\ell} is the log-density of one observation, \eqn{\mu} the mean
#' direction, \eqn{\rho \in (0,1)} the mean resultant length, \eqn{\kappa} the
#' concentration and \eqn{A(\kappa) = I_1(\kappa)/I_0(\kappa)}.
#'
#' @seealso [distrib_gradient.VonMises2Distrib()] for the score,
#'   [distrib_expected_hessian.VonMises2Distrib()] for the expectation of this
#'   quantity, [distrib_hessian.VonMises1Distrib()] for the same quantity in
#'   the concentration, and [distrib_hessian()] for the generic.
#'
#' @examples
#' d2 <- vonmises2_distrib()
#' y <- c(-1, 0, 0.5, 2)
#' th <- list(mu = 0.5, rho = 0.7)
#' h <- distrib_hessian(d2, y, th)
#' names(h)
#'
#' # numDeriv on the summed log-density reproduces the summed matrix.
#' fn <- function(v)
#'   sum(distrib_pdf(d2, y, list(mu = v[1], rho = v[2]), log = TRUE))
#' H <- numDeriv::hessian(fn, c(0.5, 0.7))
#' rbind(numeric = c(H[1, 1], H[2, 2], H[1, 2]),
#'       closed = c(sum(h$mu_mu), sum(h$rho_rho), sum(h$mu_rho)))
#'
#' # The pure second derivative varies with the data here, where in the
#' # concentration parametrization it does not.
#' h$rho_rho
S7::method(distrib_hessian, VonMises2Distrib) <- function(distrib, y, theta,
                                                          scale = c("parameter", "link"), ...) {
  p <- vm2_parts(theta)
  d <- y - theta[[1]]
  list(
    mu_mu = -p$kappa * cos(d),
    rho_rho = -p$ad$d1 * p$kd$d1^2 + (cos(d) - p$ad$A) * p$kd$d2,
    mu_rho = sin(d) * p$kd$d1
  )
}

#' @title von Mises Expected Hessian in the Resultant Length
#' @name distrib_expected_hessian.VonMises2Distrib
#' @description
#' Returns the expectation of the observed Hessian under the model, in closed
#' form and with no quadrature or simulation. Since
#' \eqn{\mathbb{E}[\cos(Y-\mu)] = A(\kappa)} and
#' \eqn{\mathbb{E}[\sin(Y-\mu)] = 0}, the term carrying the map's second
#' derivative drops out and
#' \deqn{\mathbb{E}[\ell^{(\mu\mu)}] = -\kappa A(\kappa), \qquad
#'       \mathbb{E}[\ell^{(\mu\rho)}] = 0, \qquad
#'       \mathbb{E}[\ell^{(\rho\rho)}] = -\dfrac{1}{A'(\kappa)}.}
#'
#' The last equality is the reparametrization identity: the information in
#' \eqn{\kappa} is \eqn{A'(\kappa)} and \eqn{\kappa'(\rho) = 1/A'(\kappa)}, so
#' \eqn{A'(\kappa)\{\kappa'(\rho)\}^2 = 1/A'(\kappa)}. The information in
#' \eqn{\rho} is therefore the **reciprocal** of the information in
#' \eqn{\kappa}, which a one-to-one reparametrization of a single parameter
#' whose Jacobian is that reciprocal must give. The two parameters stay
#' orthogonal.
#'
#' `approx` and `nsim` are ignored, the expectation being exact.
#'
#' @param distrib A `VonMises2Distrib` object, from [vonmises2_distrib()].
#' @param y A numeric vector of angles. Only its length is read, the
#'   expectation not depending on the data; the values are ignored.
#' @param theta A named list with components `mu` and `rho`, each a numeric
#'   vector of length 1. `rho` must lie in \eqn{(0, 1)}.
#' @param scale One of `"parameter"` (the default) or `"link"`, matched by
#'   [base::match.arg()]. Read by the generic, not by this method.
#' @param approx Ignored here, the expectation being exact. Accepted so that
#'   the signature matches the generic's, where it selects between
#'   `"bartlett"`, `"integrate"`, `"mc"` and `"opg"`.
#' @param nsim Ignored here, for the same reason. Defaults to `10000`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A named list of three numeric vectors, `mu_mu`, `rho_rho` and
#'   `mu_rho`, each of length `length(y)` and constant along it. `mu_rho` is
#'   exactly zero.
#'
#' @section Notation:
#' \eqn{\ell} is the log-density of one observation, \eqn{\mu} the mean
#' direction, \eqn{\rho \in (0,1)} the mean resultant length, \eqn{\kappa} the
#' concentration and \eqn{A(\kappa) = I_1(\kappa)/I_0(\kappa)}.
#'
#' @seealso [distrib_hessian.VonMises2Distrib()] for the quantity this is the
#'   expectation of, [distrib_expected_hessian.VonMises1Distrib()] for the same
#'   quantity in the concentration, and [distrib_expected_hessian()] for the
#'   generic.
#'
#' @examples
#' d2 <- vonmises2_distrib()
#' th <- list(mu = 0.5, rho = 0.7)
#' eh <- distrib_expected_hessian(d2, c(-1, 0, 0.5, 2), th)
#' vapply(eh, function(v) v[1], numeric(1))
#'
#' # The reparametrization identity: the information in rho is the reciprocal
#' # of the information in kappa.
#' k <- numericals7::bessel_i_ratio_inverse(0.7)$kappa
#' Ap <- numericals7::bessel_i_ratio_derivs(k)$d1
#' c(supplied = eh$rho_rho[1], reciprocal = -1 / Ap)
#'
#' # Averaging the observed Hessian over draws reaches the same three numbers.
#' set.seed(1)
#' z <- distrib_rng(d2, 3e5, th)
#' vapply(distrib_hessian(d2, z, th), mean, numeric(1))
S7::method(distrib_expected_hessian, VonMises2Distrib) <- function(distrib, y, theta,
                                                                    scale = c("parameter", "link"),
                                                                    approx = c("bartlett", "integrate", "mc", "opg"),
                                                                    nsim = 10000, ...) {
  p <- vm2_parts(theta)
  n <- length(y)
  list(
    mu_mu = rep(-p$kappa * p$ad$A, length.out = n),
    rho_rho = rep(-p$ad$d1 * p$kd$d1^2, length.out = n),
    mu_rho = rep(0, length.out = n)
  )
}

#' @title Mean of a von Mises in the Resultant Length
#' @name mean.VonMises2Distrib
#' @description
#' Returns the ordinary expectation of \eqn{Y} as a number on
#' \eqn{[-\pi, \pi)}, obtained numerically by delegating to
#' [mean.distrib()] at the implied concentration.
#'
#' It is **not** a circular quantity, and neither parameter describes it.
#' \eqn{\mu} is the mean *direction* and \eqn{\rho} the mean resultant length;
#' \eqn{\mathbb{E}[Y]} differs from \eqn{\mu} whenever \eqn{\mu \ne 0}, because
#' the interval is cut at \eqn{\pm\pi} and the density is not symmetric about
#' \eqn{\mu} on it. Compute the circular mean of a sample as
#' `atan2(mean(sin(z)), mean(cos(z)))`, which recovers \eqn{\mu}.
#'
#' @param x A `VonMises2Distrib` object, from [vonmises2_distrib()]. The
#'   argument is named `x` because the generic is [base::mean()].
#' @param theta A named list with components `mu` and `rho`, each a numeric
#'   vector of length 1. `rho` must lie in \eqn{(0, 1)}.
#' @param ... Passed on to [mean.distrib()], which reads the quadrature's
#'   settings.
#'
#' @return A numeric vector of length 1, the ordinary mean of \eqn{Y}.
#'
#' @seealso [mean.distrib()], which this calls;
#'   [vonmises1_distrib()], whose page discusses the same distinction; and
#'   [variance()] for the other ordinary moments.
#' @keywords internal
S7::method(mean, VonMises2Distrib) <- function(x, theta, ...) {
  mean(vonmises1_distrib(), list(mu = theta[[1]], kappa = vm2_parts(theta)$kappa), ...)
}


#' von Mises Distribution, Mean Direction and Mean Resultant Length
#'
#' @description
#' Builds the distribution object for the von Mises family parametrized by its
#' mean direction \eqn{\mu} and its **mean resultant length**
#' \eqn{\rho \in (0, 1)}. The returned object carries closed-form derivatives
#' of the log-density to fourth order, a closed-form expected information, and
#' the same Bessel-series distribution function the concentration
#' parametrization uses.
#'
#' The concentration \eqn{\kappa} of [vonmises1_distrib()] is unbounded and
#' hard to read. The resultant length is bounded, is the quantity circular
#' statistics reports, and is one minus the circular variance; a sample reads
#' it back directly as
#' `sqrt(mean(cos(z))^2 + mean(sin(z))^2)`.
#'
#' @param link_mu A `link` object from `linkfunctions7` for the mean direction
#'   \eqn{\mu}. Defaults to `linkfunctions7::bounded_link(lwr = -pi, upr = pi)`.
#'   The chart keeps \eqn{\mu} identified at the cost that a fit cannot walk
#'   across \eqn{\pm\pi}; see [vonmises1_distrib()].
#' @param link_rho A `link` object from `linkfunctions7` for the resultant
#'   length \eqn{\rho}. Defaults to [linkfunctions7::logit_link()], the natural
#'   link onto \eqn{(0, 1)}.
#'
#' @details
#' # The parametrization, and why it is a family of its own
#'
#' The density on \eqn{y \in [-\pi, \pi)} is
#' \deqn{f(y; \mu, \rho) = \frac{e^{\kappa\cos(y-\mu)}}{2\pi I_{0}(\kappa)},
#'       \qquad \kappa = A^{-1}(\rho),}
#' with \eqn{A(\kappa) = I_1(\kappa)/I_0(\kappa)}, and
#' \eqn{\mathbb{E}[\cos(Y-\mu)] = \rho}. The map \eqn{\rho = A(\kappa)} is a
#' strictly increasing bijection from \eqn{(0, \infty)} onto \eqn{(0, 1)}.
#'
#' **Its inverse has no closed form**, which is why this is a family of its own
#' and not a [reparametrize()] of the other. \eqn{\kappa} is obtained by root
#' finding on \eqn{\log\kappa}, and its four derivatives come from the inverse
#' function rule applied to \eqn{A'} through \eqn{A^{(4)}}, which the Bessel
#' recurrences give from the same two evaluations \eqn{A} already needs.
#'
#' The map is steep near \eqn{\rho = 1}: measured,
#' \eqn{\mathrm{d}\kappa/\mathrm{d}\rho} is 2.00 at \eqn{\rho = 0.01}, 3.14 at
#' 0.5, 199 at 0.95 and \eqn{5.0\times10^{5}} at 0.999. A fit whose data are
#' nearly all in one direction is therefore far better conditioned in
#' \eqn{\kappa} than in \eqn{\rho}.
#'
#' # Derivatives and information
#'
#' The map touches the **second parameter only**, so every chain rule is the
#' one-variable one and the derivatives are exact at every order. Every
#' component carrying at least one \eqn{\mu} collapses to a single term, the
#' concentration parametrization's \eqn{\mu}-derivatives being linear in
#' \eqn{\kappa}.
#'
#' The expected information is closed form and the two parameters are
#' orthogonal, as in the concentration parametrization, with
#' \deqn{\mathbb{E}[\ell^{(\mu\mu)}] = -\kappa\rho, \qquad
#'       \mathbb{E}[\ell^{(\mu\rho)}] = 0, \qquad
#'       \mathbb{E}[\ell^{(\rho\rho)}] = -\dfrac{1}{A'(\kappa)},}
#' the last being the reciprocal of the information in \eqn{\kappa}, since the
#' Jacobian of the map is that reciprocal.
#'
#' # The moments are not the parameters
#'
#' [mean.VonMises2Distrib()] returns the ordinary expectation of \eqn{Y} on
#' \eqn{[-\pi, \pi)}, which differs from \eqn{\mu} whenever \eqn{\mu \ne 0};
#' [vonmises1_distrib()] discusses the same distinction. What this
#' parametrization buys is that the **second** parameter is a quantity a sample
#' reads back directly, which the concentration is not.
#'
#' @section Notation:
#' \eqn{\ell} is the log-density of one observation, \eqn{\mu} the mean
#' direction, \eqn{\rho \in (0,1)} the mean resultant length, \eqn{\kappa} the
#' concentration, \eqn{I_m} the modified Bessel function of the first kind of
#' order \eqn{m}, and \eqn{A(\kappa) = I_1(\kappa)/I_0(\kappa)}. \eqn{\eta} is
#' a parameter on the unconstrained scale of its link, with
#' \eqn{\theta = g^{-1}(\eta)}.
#'
#' @return An S7 object of class `VonMises2Distrib`, inheriting from
#'   `continuous_distrib`, with `distrib_name` `"von mises2"`, `dimension`
#'   `"univariate"`, `bounds` `c(-pi, pi)`, `params` `c("mu", "rho")`,
#'   `n_params` `2`, `params_bounds` the domains \eqn{(-\pi, \pi)} and
#'   \eqn{(0, 1)}, and `link_params` the two links given here.
#'
#' @references
#' Mardia, K. V. and Jupp, P. E. (2000). *Directional Statistics*, Chapter 3.
#' Wiley, Chichester.
#'
#' @importFrom linkfunctions7 bounded_link logit_link
#'
#' @examples
#' d <- vonmises2_distrib()
#' d
#'
#' # The resultant length is bounded, which is what makes it readable.
#' d@params_bounds$rho
#'
#' # The same law as the concentration parametrization at the implied kappa.
#' th <- list(mu = 0.5, rho = 0.7)
#' k <- numericals7::bessel_i_ratio_inverse(0.7)$kappa
#' all.equal(distrib_pdf(d, c(-1, 0, 1), th),
#'           distrib_pdf(vonmises1_distrib(), c(-1, 0, 1),
#'                       list(mu = 0.5, kappa = k)))
#'
#' # A sample reads the second parameter back directly, which is the point.
#' set.seed(1)
#' z <- distrib_rng(d, 3e5, th)
#' c(resultant = sqrt(mean(cos(z))^2 + mean(sin(z))^2), rho = 0.7)
#'
#' # The map is steep near one, so a nearly deterministic direction is better
#' # conditioned in kappa than in rho.
#' vapply(c(0.01, 0.5, 0.95, 0.999),
#'        function(r) numericals7::bessel_i_ratio_inverse(r)$d1, numeric(1))
#'
#' # Fitting recovers both parameters.
#' set.seed(3)
#' coef(fit_distrib(d, distrib_rng(d, 2000, list(mu = 0.8, rho = 0.6))))
#'
#' @seealso
#' [vonmises1_distrib()] for the same law in the concentration, which is
#' better conditioned at a nearly deterministic direction;
#' [numericals7::bessel_i_ratio_inverse()] for the map and its derivatives;
#' [beta1_distrib()] for another family on \eqn{(0, 1)};
#' [fit_distrib()] to estimate the parameters; [check_distrib()] to validate a
#' family of your own against the same battery this one passes;
#' [VonMises2Distrib] for the class.
#' @export
vonmises2_distrib <- function(link_mu = bounded_link(lwr = -pi, upr = pi),
                              link_rho = logit_link()) {
  VonMises2Distrib(
    distrib_name = "von mises2",
    dimension = "univariate",
    bounds = c(-pi, pi),
    params = c("mu", "rho"),
    params_interpretation = c(mu = "mean direction",
                              rho = "mean resultant length"),
    n_params = 2,
    params_bounds = list(mu = c(-pi, pi), rho = c(0, 1)),
    link_params = list(mu = link_mu, rho = link_rho)
  )
}


#' @title von Mises Third-Order Derivatives in the Resultant Length
#' @name distrib_deriv3.VonMises2Distrib
#' @description
#' Computes the four distinct third derivatives of the log-density in
#' \eqn{\mu} and \eqn{\rho}, in closed form. Every component carrying at least
#' one \eqn{\mu} collapses to a **single term**
#' \eqn{D_a \kappa^{(b)}(\rho)}, with \eqn{D_a} the \eqn{a}-th
#' \eqn{\mu}-derivative of \eqn{\cos(y-\mu)} and \eqn{\kappa^{(b)}} the
#' \eqn{b}-th derivative of \eqn{A^{-1}}: the concentration parametrization's
#' \eqn{\mu}-derivatives are linear in \eqn{\kappa}, so the composition has
#' nothing to expand. The pure-\eqn{\rho} component carries the full
#' one-variable Faa di Bruno on \eqn{\log I_0}, written out.
#'
#' With `expected = TRUE` the method calls [expected_derivative()], which is
#' the one place on this page where `approx` and `nsim` are read.
#'
#' @param distrib A `VonMises2Distrib` object, from [vonmises2_distrib()].
#' @param y A numeric vector of angles in \eqn{[-\pi, \pi)}. With
#'   `expected = TRUE` only its length is read.
#' @param theta A named list with components `mu` and `rho`, each a numeric
#'   vector of length 1 or of the length of `y`. `rho` must lie in
#'   \eqn{(0, 1)}.
#' @param scale One of `"parameter"` (the default) or `"link"`, matched by
#'   [base::match.arg()]. Read by the generic, not by this method. **Note the
#'   argument order**: `scale` precedes `expected` here, so both are best given
#'   by name.
#' @param expected Logical of length 1. When `TRUE` the expectation under the
#'   model is returned in place of the value at the data, computed numerically.
#'   Defaults to `FALSE`.
#' @param approx One of `"integrate"` (the default here), `"bartlett"`, `"mc"`
#'   or `"opg"`. Read only when `expected = TRUE`.
#' @param nsim A single positive integer, the sample size when
#'   `approx = "mc"`. Read only when `expected = TRUE`. Defaults to `10000`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A named list of four numeric vectors, `mu_mu_mu`, `mu_mu_rho`,
#'   `mu_rho_rho` and `rho_rho_rho`, each of length `length(y)`.
#'
#' @section Notation:
#' \eqn{\ell} is the log-density of one observation, \eqn{\mu} the mean
#' direction, \eqn{\rho \in (0,1)} the mean resultant length, \eqn{\kappa} the
#' concentration and \eqn{A(\kappa) = I_1(\kappa)/I_0(\kappa)}.
#'
#' @seealso [distrib_hessian.VonMises2Distrib()] for the order below,
#'   [distrib_deriv4.VonMises2Distrib()] for the order above,
#'   [vm2_parts()] for the map's derivatives, and [distrib_deriv3()] for the
#'   generic.
#'
#' @examples
#' d2 <- vonmises2_distrib()
#' y <- c(-1, 0, 0.5, 2)
#' th <- list(mu = 0.5, rho = 0.7)
#' d3 <- distrib_deriv3(d2, y, th)
#' names(d3)
#'
#' # A central difference of the Hessian reproduces the pure-rho component,
#' # which is the one carrying the whole chain rule.
#' eps <- 1e-5
#' up <- distrib_hessian(d2, y, list(mu = 0.5, rho = 0.7 + eps))$rho_rho
#' dn <- distrib_hessian(d2, y, list(mu = 0.5, rho = 0.7 - eps))$rho_rho
#' all.equal((up - dn) / (2 * eps), d3$rho_rho_rho, tolerance = 1e-5)
#'
#' # Unlike the concentration parametrization, no component is exactly zero:
#' # the map's higher derivatives bring the data into all four.
#' vapply(d3, function(v) v[1], numeric(1))
S7::method(distrib_deriv3, VonMises2Distrib) <- function(distrib, y, theta,
                                                         scale = c("parameter", "link"),
                                                         expected = FALSE,
                                                         approx = c("integrate", "bartlett", "mc", "opg"),
                                                         nsim = 10000, ...) {
  if (expected) {
    return(expected_derivative(distrib, y, theta, order = 3L,
                               approx = match.arg(approx), nsim = nsim))
  }
  p <- vm2_parts(theta)
  d <- y - theta[[1]]
  A <- p$ad
  k <- p$kd
  # d^3/drho^3 of log I_0(kappa(rho)), Faa di Bruno written out
  phi3 <- A$d2 * k$d1^3 + 3 * A$d1 * k$d1 * k$d2 + A$A * k$d3
  list(mu_mu_mu = -sin(d) * p$kappa,
       mu_mu_rho = -cos(d) * k$d1,
       mu_rho_rho = sin(d) * k$d2,
       rho_rho_rho = cos(d) * k$d3 - phi3)
}

#' @title von Mises Fourth-Order Derivatives in the Resultant Length
#' @name distrib_deriv4.VonMises2Distrib
#' @description
#' Computes the five distinct fourth derivatives of the log-density in
#' \eqn{\mu} and \eqn{\rho}, in closed form, by the construction
#' [distrib_deriv3.VonMises2Distrib()] describes carried one order further: a
#' single term \eqn{D_a \kappa^{(b)}(\rho)} for every component carrying a
#' \eqn{\mu}, and the fourth-order one-variable Faa di Bruno on
#' \eqn{\log I_0} for the pure-\eqn{\rho} one.
#'
#' With `expected = TRUE` the method calls [expected_derivative()], which is
#' the one place on this page where `approx` and `nsim` are read.
#'
#' @param distrib A `VonMises2Distrib` object, from [vonmises2_distrib()].
#' @param y A numeric vector of angles in \eqn{[-\pi, \pi)}. With
#'   `expected = TRUE` only its length is read.
#' @param theta A named list with components `mu` and `rho`, each a numeric
#'   vector of length 1 or of the length of `y`. `rho` must lie in
#'   \eqn{(0, 1)}.
#' @param scale One of `"parameter"` (the default) or `"link"`, matched by
#'   [base::match.arg()]. Read by the generic, not by this method. **Note the
#'   argument order**: `scale` precedes `expected` here, so both are best given
#'   by name.
#' @param expected Logical of length 1. When `TRUE` the expectation under the
#'   model is returned in place of the value at the data, computed numerically.
#'   Defaults to `FALSE`.
#' @param approx One of `"integrate"` (the default here), `"bartlett"`, `"mc"`
#'   or `"opg"`. Read only when `expected = TRUE`.
#' @param nsim A single positive integer, the sample size when
#'   `approx = "mc"`. Read only when `expected = TRUE`. Defaults to `10000`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A named list of five numeric vectors, `mu_mu_mu_mu`,
#'   `mu_mu_mu_rho`, `mu_mu_rho_rho`, `mu_rho_rho_rho` and `rho_rho_rho_rho`,
#'   each of length `length(y)`.
#'
#' @section Notation:
#' \eqn{\ell} is the log-density of one observation, \eqn{\mu} the mean
#' direction, \eqn{\rho \in (0,1)} the mean resultant length, \eqn{\kappa} the
#' concentration and \eqn{A(\kappa) = I_1(\kappa)/I_0(\kappa)}.
#'
#' @seealso [distrib_deriv3.VonMises2Distrib()] for the order below and the
#'   construction, [vm2_parts()] for the map's derivatives, and
#'   [distrib_deriv4()] for the generic.
#'
#' @examples
#' d2 <- vonmises2_distrib()
#' y <- c(-1, 0, 0.5, 2)
#' th <- list(mu = 0.5, rho = 0.7)
#' d4 <- distrib_deriv4(d2, y, th)
#' names(d4)
#'
#' # A central difference of the third order reproduces the pure-direction
#' # component.
#' eps <- 1e-5
#' up <- distrib_deriv3(d2, y, list(mu = 0.5 + eps, rho = 0.7))$mu_mu_mu
#' dn <- distrib_deriv3(d2, y, list(mu = 0.5 - eps, rho = 0.7))$mu_mu_mu
#' all.equal((up - dn) / (2 * eps), d4$mu_mu_mu_mu, tolerance = 1e-5)
S7::method(distrib_deriv4, VonMises2Distrib) <- function(distrib, y, theta,
                                                         scale = c("parameter", "link"),
                                                         expected = FALSE,
                                                         approx = c("integrate", "bartlett", "mc", "opg"),
                                                         nsim = 10000, ...) {
  if (expected) {
    return(expected_derivative(distrib, y, theta, order = 4L,
                               approx = match.arg(approx), nsim = nsim))
  }
  p <- vm2_parts(theta)
  d <- y - theta[[1]]
  A <- p$ad
  k <- p$kd
  phi4 <- A$d3 * k$d1^4 + 6 * A$d2 * k$d1^2 * k$d2 +
    A$d1 * (3 * k$d2^2 + 4 * k$d1 * k$d3) + A$A * k$d4
  list(mu_mu_mu_mu = cos(d) * p$kappa,
       mu_mu_mu_rho = -sin(d) * k$d1,
       mu_mu_rho_rho = -cos(d) * k$d2,
       mu_rho_rho_rho = sin(d) * k$d3,
       rho_rho_rho_rho = cos(d) * k$d4 - phi4)
}


#' @title von Mises Distribution Function in the Resultant Length
#' @name distrib_cdf.VonMises2Distrib
#' @description
#' Computes \eqn{F(q) = P(Y \le q)} on \eqn{[-\pi, \pi)} from the Fourier
#' series of the concentration parametrization, read at
#' \eqn{\kappa = A^{-1}(\rho)}. The map touches the second parameter only and
#' the response not at all, so the distribution function is the other family's
#' at that concentration.
#'
#' What it replaces is the base class's quadrature, one integration per
#' observation. See [vm_cdf()] for the series, its measured term count and the
#' blocking that keeps the intermediate small.
#'
#' @param distrib A `VonMises2Distrib` object, from [vonmises2_distrib()].
#' @param q A numeric vector of angles. Below \eqn{-\pi} the value is 0 and at
#'   or above \eqn{\pi} it is 1.
#' @param theta A named list with components `mu` and `rho`, each a numeric
#'   vector of length 1 or of the length of `q`. `mu` must lie in
#'   \eqn{(-\pi, \pi)} and `rho` in \eqn{(0, 1)}.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'   This method takes **no** `lower.tail` or `log.p`: the upper tail is
#'   `1 - F(q)` and the logarithm is `log(F(q))`.
#'
#' @return A numeric vector of probabilities in \eqn{[0, 1]}, of length
#'   `max(length(q), length(mu), length(rho))`.
#'
#' @seealso [vm_cdf()] for the series,
#'   [distrib_cdf.VonMises1Distrib()] for the same quantity in the
#'   concentration, and [distrib_cdf()] for the generic.
#'
#' @examples
#' d2 <- vonmises2_distrib()
#' th <- list(mu = 0.5, rho = 0.7)
#' y <- c(-1, 0, 0.5, 2)
#'
#' # The series agrees with a direct quadrature of the density.
#' rbind(series = distrib_cdf(d2, y, th),
#'       quadrature = vapply(y, function(v)
#'         integrate(function(u) distrib_pdf(d2, u, th), -pi, v)$value,
#'         numeric(1)))
#'
#' # And with the concentration parametrization at the implied concentration.
#' k <- numericals7::bessel_i_ratio_inverse(0.7)$kappa
#' all.equal(distrib_cdf(d2, y, th),
#'           distrib_cdf(vonmises1_distrib(), y, list(mu = 0.5, kappa = k)))
#' @keywords internal
S7::method(distrib_cdf, VonMises2Distrib) <- function(distrib, q, theta,
                                                      ...) {
  vm_cdf(q, theta[[1]], vm2_parts(theta)$kappa)
}
