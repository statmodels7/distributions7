#' @include cdf_higher.R
NULL

# Third and fourth derivatives of the distribution function for the families
# that reach them through another family.
#
# A family written as a map of another carries the parent's cdf derivatives
# through one Faa di Bruno pass. Orders one and two do this already, in
# `chain_cdf_deriv`, which is written out; at three and four the same sum comes
# from `chain_assemble`, the enumeration the reparametrized parameter
# derivatives already run on, so there is no second copy of the partition sum.

#' The Chain Rule of Any Order on a Parent's CDF Derivatives
#'
#' @description
#' Carries the parent's derivatives of \eqn{F} in its own parameters onto the
#' new ones, for a map given as keyed partial tables.
#'
#' @details
#' The parent's tables are taken on the natural scale, \code{lower.tail = TRUE}
#' and \code{log = FALSE}, because the chain rule applies to \eqn{F} itself; the
#' tail and the logarithm are put on afterwards by \code{\link{cdf_scale_k}}.
#'
#' @param parent The distribution being mapped.
#' @param q A numeric vector of quantiles.
#' @param th_par The parent's parameters at the new ones.
#' @param maps The map's keyed partial tables.
#' @param new_params The new parameter names.
#' @param order The derivative order, 1 to 4.
#'
#' @return A named list of derivative components of \eqn{F}.
#'
#' @seealso \code{\link{chain_assemble}}, \code{\link{chain_cdf_deriv}}
#' @keywords internal
chain_cdf_deriv_k <- function(parent, q, th_par, maps, new_params, order) {
  D <- lapply(seq_len(order), function(k) {
    switch(k,
      distrib_grad_cdf(parent, q, th_par, lower.tail = TRUE, log = FALSE),
      distrib_hess_cdf(parent, q, th_par, lower.tail = TRUE, log = FALSE),
      distrib_deriv3_cdf(parent, q, th_par, lower.tail = TRUE, log = FALSE),
      distrib_deriv4_cdf(parent, q, th_par, lower.tail = TRUE, log = FALSE)
    )
  })
  chain_assemble(D, parent@params, maps, new_params, order, length(q))
}

#' Third and Fourth Log-CDF Derivatives of a Mapped Family
#'
#' @description
#' The chain rule on the parent's when the parent's are exact at every order up
#' to the one asked for, and the stencil otherwise.
#'
#' @details
#' The gate is the one orders one and two use. A chain rule carrying a
#' differenced quantity would report a closed form and deliver the parent's
#' noise, and the truncation wrapper reads that distinction to choose its own
#' route.
#'
#' @param distrib The mapped distribution.
#' @param parent The distribution being mapped.
#' @param th_par The parent's parameters at the new ones.
#' @param maps The map's keyed partial tables.
#' @param q A numeric vector of quantiles.
#' @param theta A named list of the new parameters.
#' @param order The derivative order, 3 or 4.
#' @param lower.tail Logical; whether the lower tail is wanted.
#' @param log Logical; whether derivatives of the log probability are wanted.
#' @param q_par The point at which to evaluate the parent, when the two
#'   families are related by a monotone transformation of the response as well
#'   as by a map of the parameters. A lognormal is a gaussian at \eqn{\log q},
#'   and since the transformation carries no parameter the derivatives in
#'   \eqn{\theta} are the parent's with the point substituted.
#'
#' @return A named list of derivative component vectors.
#'
#' @seealso \code{\link{mapped_cdf_deriv}}
#' @keywords internal
mapped_cdf_deriv_k <- function(distrib, parent, th_par, maps, q, theta, order,
                               lower.tail, log, q_par = q) {
  exact <- all(vapply(seq_len(order),
                      function(k) has_exact_cdf_deriv(parent, k), logical(1)))
  tabs <- if (exact) {
    lapply(seq_len(order),
           function(k) chain_cdf_deriv_k(parent, q_par, th_par, maps,
                                         distrib@params, k))
  } else {
    cdf_tables(distrib, q, theta, order)
  }
  cdf_scale_k(distrib, distrib_cdf(distrib, q, theta), tabs, order,
              lower.tail, log)
}

#' Register the Two New Orders on a Mapped Family
#'
#' @description
#' Turns the parent and the map into the two methods, so that a family states
#' its map once instead of twice.
#'
#' @param cls The S7 class.
#' @param parent_fn A function of no arguments returning the parent.
#' @param th_fn A function of \code{theta} returning the parent's parameters.
#' @param md_fn The map's table function.
#' @param q_fn The transformation of the response, when the parent is the same
#'   law on a transformed scale. The identity by default.
#' @param orders The orders to register, 3 and 4 by default. A family whose
#'   written-out route stops below the fourth order takes the rest here.
#'
#' @return Invisibly \code{NULL}; called for the registration.
#'
#' @keywords internal
register_mapped_cdf_k <- function(cls, parent_fn, th_fn, md_fn,
                                  q_fn = identity, orders = 3:4) {
  make <- function(o) {
    force(o)
    function(distrib, q, theta, lower.tail = TRUE, log = TRUE, ...) {
      mapped_cdf_deriv_k(distrib, parent_fn(), th_fn(theta), md_fn(theta),
                         q, theta, o, lower.tail, log, q_par = q_fn(q))
    }
  }
  gens <- list(distrib_grad_cdf, distrib_hess_cdf,
               distrib_deriv3_cdf, distrib_deriv4_cdf)
  for (o in orders) S7::method(gens[[o]], cls) <- make(o)
  invisible(NULL)
}

#' @title Third and Fourth Log-CDF Derivatives of a Reparametrized Distribution
#' @name distrib_deriv3_cdf.ReparamContinuousDistrib
#' @description
#' The chain rule on the parent's, exact whenever the parent's are and the
#' stencil otherwise, as at the two orders below.
#' @param distrib A reparametrized distribution.
#' @param q A numeric vector of quantiles.
#' @param theta A named list of the new parameters.
#' @param lower.tail Logical; if \code{TRUE} (default), the lower tail.
#' @param log Logical; if \code{TRUE} (default), derivatives of the log probability.
#' @param ... Unused.
#' @return A named list, one vector per component.
#' @seealso \code{\link{reparametrize}}
#' @keywords internal
S7::method(distrib_deriv3_cdf, ReparamContinuousDistrib) <-
  function(distrib, q, theta, lower.tail = TRUE, log = TRUE, ...) {
    mapped_cdf_deriv_k(distrib, distrib@parent_distrib,
                       reparam_theta(distrib, theta),
                       reparam_tables(distrib, theta),
                       q, theta, 3L, lower.tail, log)
  }

#' @rdname distrib_deriv3_cdf.ReparamContinuousDistrib
#' @name distrib_deriv4_cdf.ReparamContinuousDistrib
#' @keywords internal
S7::method(distrib_deriv4_cdf, ReparamContinuousDistrib) <-
  function(distrib, q, theta, lower.tail = TRUE, log = TRUE, ...) {
    mapped_cdf_deriv_k(distrib, distrib@parent_distrib,
                       reparam_theta(distrib, theta),
                       reparam_tables(distrib, theta),
                       q, theta, 4L, lower.tail, log)
  }

# the three families the package writes as a map of another
register_mapped_cdf_k(Gaussian2Distrib, gaussian1_distrib,
                      function(theta) list(mu = theta[[1]],
                                           sigma = sqrt(theta[[2]])),
                      md_gaussian2)
register_mapped_cdf_k(Gaussian3Distrib, gaussian1_distrib,
                      function(theta) list(mu = theta[[1]],
                                           sigma = 1 / sqrt(theta[[2]])),
                      md_gaussian3)
# the second inverse-gaussian parametrization also takes its Hessian here: the
# written-out route registered only the gradient, which was right while the
# parent differenced its own second order and is not now that it does not
register_mapped_cdf_k(InvGauss2Distrib, invgauss1_distrib,
                      function(theta) list(mu = theta[[1]],
                                           phi = 1 / theta[[2]]),
                      md_invgauss2, orders = 2:4)

# the second Laplace parametrization carries the rate, so it is the first at
# sigma = 1/lambda
register_mapped_cdf_k(Laplace2Distrib, laplace_distrib,
                      function(theta) list(mu = theta[[1]],
                                           sigma = 1 / theta[[2]]),
                      md_laplace2)

# a lognormal is a gaussian at log q, and the transformation carries no
# parameter, so the derivatives in the parameters are the gaussian's with the
# point substituted and the scale read off the variance
register_mapped_cdf_k(Lognormal1Distrib, gaussian1_distrib,
                      function(theta) list(mu = theta[[1]],
                                           sigma = sqrt(theta[[2]])),
                      md_gaussian2, q_fn = log)

# the Gumbel is location-scale, like the four already carrying these orders
S7::method(distrib_deriv3_cdf, GumbelDistrib) <- loc_scale_deriv_cdf_k(3L)
S7::method(distrib_deriv4_cdf, GumbelDistrib) <- loc_scale_deriv_cdf_k(4L)


# --- families location-scale in their first two parameters -------------------

#' Higher Log-CDF Derivatives When Only Some Parameters Are Location-Scale
#'
#' @description
#' The higher-order companion of \code{\link{partial_loc_scale_hess_cdf}}: the
#' components over the location and the scale from the location-scale
#' construction, and the components naming a shape parameter from the stencil.
#'
#' @param order The derivative order, 3 or 4.
#'
#' @return A function suitable for registering as an S7 method.
#'
#' @seealso \code{\link{loc_scale_cdf_deriv_k}}
#' @keywords internal
partial_loc_scale_deriv_cdf_k <- function(order) {
  function(distrib, q, theta, lower.tail = TRUE, log = TRUE, ...) {
    params <- distrib@params
    tabs <- lapply(seq_len(order), function(k) {
      # the stencil is taken over every component and the closed ones then
      # replace it. Computing the location and scale twice is waste, and it is
      # immaterial here: the whole cdf surface costs milliseconds at a thousand
      # quantiles, so the alternative -- widening the stencil's signature to
      # take a subset -- would buy nothing and touch a shared function.
      out <- numerical_cdf_deriv_k(distrib, q, theta, k)
      closed <- loc_scale_cdf_deriv_k(distrib, q, theta, k)
      out[names(closed)] <- closed
      out[deriv_names(params, k)]
    })
    cdf_scale_k(distrib, distrib_cdf(distrib, q, theta), tabs, order,
                lower.tail, log)
  }
}

for (.cls in list(StudentT1Distrib, PseudoHuberDistrib, SkewNormal1Distrib,
                  SkewTDistrib)) {
  S7::method(distrib_deriv3_cdf, .cls) <- partial_loc_scale_deriv_cdf_k(3L)
  S7::method(distrib_deriv4_cdf, .cls) <- partial_loc_scale_deriv_cdf_k(4L)
}
rm(.cls)
