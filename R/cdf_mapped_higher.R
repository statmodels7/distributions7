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
#' Carries a parent's derivatives of \eqn{F} in its own parameters onto the new
#' ones, for a map given as keyed partial tables, at any order up to four. The
#' general form of [chain_cdf_deriv()], which writes orders 1 and 2 out; above
#' them the same sum comes from [chain_assemble()], the enumeration the
#' reparametrized parameter derivatives already run on, so the package carries
#' no second copy of the partition sum.
#'
#' @details
#' The parent's tables are fetched on the natural scale, with
#' `lower.tail = TRUE` and `log = FALSE`, because the chain rule applies to
#' \eqn{F} itself. The tail and the logarithm are put on afterwards by
#' [cdf_scale_k()], so one Faa di Bruno pass serves both tails and both
#' scales.
#'
#' @section Notation:
#' \eqn{F} is the parent's distribution function, \eqn{\psi} the new
#' parameters, \eqn{h} the map from them to the parent's, and \eqn{\partial^I}
#' a derivative with respect to a multi-index.
#'
#' @param parent The distribution being mapped.
#' @param q A numeric vector of quantiles, already on the parent's scale if the
#'   map transforms the response.
#' @param th_par The parent's parameters, evaluated at the new ones.
#' @param maps The map's keyed partial tables. A missing key is an exact zero.
#' @param new_params A character vector naming the new parameters, which names
#'   and orders the result.
#' @param order The derivative order, 1 to 4.
#'
#' @return A named list of numeric vectors, derivatives of \eqn{F} itself on
#'   the natural scale, keyed as
#'   [`deriv_names(new_params, order)`][deriv_names].
#'
#' @seealso [chain_cdf_deriv()] for orders 1 and 2;
#'   [chain_assemble()] for the partition sum;
#'   [mapped_cdf_deriv_k()], the caller that gates it.
#'
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

#' Higher CDF Derivatives of a Family Written as a Map of Another
#'
#' @description
#' The body every mapped family registers at orders 3 and 4, and at 2 where its
#' written-out route stops there. It asks whether the parent's cdf derivatives
#' are exact at **every** order up to the one wanted, takes
#' [chain_cdf_deriv_k()] if they are, and falls back to [cdf_tables()] on the
#' new family's own cdf if they are not.
#'
#' @details
#' # The gate
#'
#' The chain is exact only if everything it carries is exact, and the higher
#' orders read the lower ones, so the test is over all of them. It is what
#' stops a family from adding an exact transformation to a differenced parent
#' and reporting the result as closed.
#'
#' # A transformed response
#'
#' A family may be the parent's law at a transformed point as well as at a
#' mapped parameter. A lognormal is a Gaussian at \eqn{\log q}, and since the
#' transformation carries no parameter, the derivatives in \eqn{\theta} are the
#' parent's with the point substituted. `q_par` is where that substitution is
#' handed in.
#'
#' @param distrib The mapped family, whose cdf and parameter names are read.
#' @param parent The distribution it maps onto.
#' @param th_par The parent's parameters, evaluated at the new ones.
#' @param maps The map's keyed partial tables.
#' @param q A numeric vector of quantiles, on the new family's own scale.
#' @param theta A named list of the new parameters.
#' @param order The derivative order, 2 to 4.
#' @param lower.tail Is the lower tail wanted? A single logical.
#' @param log Are derivatives of the log probability wanted? A single logical.
#' @param q_par The quantiles on the parent's scale. `q` by default, and
#'   `log(q)` for a lognormal.
#'
#' @return A named list of numeric vectors of the requested order, on the
#'   requested tail and scale.
#'
#' @seealso [chain_cdf_deriv_k()] for the exact route;
#'   [cdf_tables()] for the fallback;
#'   [has_exact_cdf_deriv()] for the gate;
#'   [mapped_cdf_deriv()] for orders 1 and 2.
#'
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

#' Register the Higher CDF Orders on a Mapped Family
#'
#' @description
#' Turns a parent and a map into the two or three S7 methods a mapped family
#' needs, so that the family states its map once instead of once per order.
#' Five families are registered through it: the two further Gaussian
#' parametrizations, the second inverse Gaussian, the second Laplace and the
#' lognormal.
#'
#' @details
#' `orders` is how a family takes more than the top two. The second inverse
#' Gaussian is registered at 2 to 4 because its written-out route in
#' `cdf_derivatives_families.R` stops at the gradient; that was right while its
#' parent differenced its own second order and stopped being right when the
#' parent gained a closed one.
#'
#' `force(o)` inside the factory is what keeps the registrations from sharing
#' one order.
#'
#' @param cls The S7 class to register on.
#' @param parent_fn A function of no arguments returning the parent
#'   distribution. A function rather than the object, so that the parent is
#'   built at call time; at load time the class it names may not exist yet.
#' @param th_fn A function of `theta` returning the parent's parameters.
#' @param md_fn A function of `theta` returning the map's keyed partial tables.
#' @param q_fn The transformation of the response, for a parent that is the
#'   same law on a transformed scale. [identity()] by default, and `log` for
#'   the lognormal.
#' @param orders An integer vector of the orders to register, `3:4` by default.
#'   A family whose written-out route stops below the fourth order takes the
#'   rest here.
#'
#' @return Invisibly `NULL`. Called for the registration.
#'
#' @seealso [mapped_cdf_deriv_k()], the body it registers;
#'   [distrib_deriv3_cdf()] for the generics.
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
#'
#' @description
#' The chain rule on the parent's cdf derivatives at orders 3 and 4, through
#' [mapped_cdf_deriv_k()]. It is exact whenever the parent's are exact at every
#' order up to the one wanted, and falls to one product stencil on the
#' reparametrized cdf otherwise.
#'
#' @details
#' The map's partials come from the object itself, `reparam_tables()` reading
#' the keyed tables the reparametrization was built with. A family created by
#' [reparametrize()] over a parent with closed cdf derivatives therefore
#' reaches the fourth order with no arithmetic of its own.
#'
#' @param distrib A `ReparamContinuousDistrib`, from [reparametrize()].
#' @param q A numeric vector of quantiles.
#' @param theta A named list of the new parameters, on the new parameter scale.
#' @param lower.tail Is the lower tail wanted? A single logical, `TRUE` by
#'   default.
#' @param log Are derivatives of the log probability wanted? A single logical,
#'   `TRUE` by default.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return For `distrib_deriv3_cdf()`, a named list of third-derivative
#'   components keyed as [`deriv_names(distrib@params, 3)`][deriv_names]; for
#'   `distrib_deriv4_cdf()`, the fourth-order components. Each vector is the
#'   length of `q` recycled against `theta`.
#'
#' @seealso [chain_cdf_deriv_k()] for the identity;
#'   [distrib_grad_cdf.ReparamContinuousDistrib()] for the first order;
#'   [reparametrize()].
#'
#' @examples
#' d <- reparametrize(
#'   gaussian1_distrib(),
#'   map = function(psi) list(mu = psi$mu, sigma = sqrt(psi$sigma2)),
#'   params = c("mu", "sigma2"),
#'   bounds = list(mu = c(-Inf, Inf), sigma2 = c(0, Inf)),
#'   links = list(mu = linkfunctions7::identity_link(),
#'                sigma2 = linkfunctions7::log_link())
#' )
#' th <- list(mu = 0.3, sigma2 = 1.44)
#'
#' # It agrees with the family written out by hand, to the closed route's own
#' # accuracy rather than to the stencil's.
#' a <- distrib_deriv3_cdf(d, 1, th)
#' b <- distrib_deriv3_cdf(gaussian2_distrib(), 1, th)
#' max(abs(unlist(a) - unlist(b)))
#'
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
#' Builds the [distrib_deriv3_cdf()] or [distrib_deriv4_cdf()] body the Student
#' t, the pseudo-Huber and the skew t register: the components over the
#' location and the scale from [loc_scale_cdf_deriv_k()], and every component
#' naming a shape parameter from [numerical_cdf_deriv_k()].
#'
#' @details
#' # Why the stencil is taken over everything
#'
#' The stencil runs over every component and the closed ones then overwrite it,
#' so the location and the scale are computed twice. That is deliberate: the
#' whole cdf surface costs milliseconds at a thousand quantiles, and the
#' alternative, widening the stencil's signature to take a subset at these
#' orders, would touch a shared function for no measurable gain. The orders
#' below do use `which`, where the same components are asked for far more
#' often.
#'
#' # Who is not here
#'
#' The skew normal was among these families and is not any more. Owen's T has
#' elementary partial derivatives in both of its arguments, so its shape
#' components close too and it has a route of its own in
#' `cdf_skewnormal_higher.R`.
#'
#' @param order The derivative order, 3 or 4.
#'
#' @return A function of `(distrib, q, theta, lower.tail, log, ...)` suitable
#'   for registering as an S7 method on either generic, returning a named list
#'   of numeric vectors of that order.
#'
#' @seealso [loc_scale_cdf_deriv_k()] for the closed components;
#'   [numerical_cdf_deriv_k()] for the rest;
#'   [partial_loc_scale_hess_cdf()] for the second order.
#'
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

# the skew normal is not among them any more: Owen's T has elementary partial
# derivatives, so its shape components close too and it has a route of its own
for (.cls in list(StudentT1Distrib, PseudoHuberDistrib, SkewTDistrib)) {
  S7::method(distrib_deriv3_cdf, .cls) <- partial_loc_scale_deriv_cdf_k(3L)
  S7::method(distrib_deriv4_cdf, .cls) <- partial_loc_scale_deriv_cdf_k(4L)
}
rm(.cls)
