#' @include cdf_derivatives.R reparametrize.R reparam_maps.R
#' @include gaussian2_distrib.R gaussian3_distrib.R invgauss2_distrib.R
#' @include gumbel_distrib.R skewnormal1_distrib.R skewt_distrib.R
#' @include exponential_distrib.R weibull1_distrib.R gpd_distrib.R
NULL

# ===========================================================================
# Closed-form cdf derivatives beyond the location-scale families of
# cdf_derivatives.R, and the chain rule that carries a parent's closed forms
# into a reparametrization.
#
# The governing identity is the ordinary chain rule on the composition
# F(q; h(psi)), the region of integration being independent of the parameters
# in both parametrizations:
#
#   dF/dpsi_a           = sum_k F_k h^k_a
#   d2F/dpsi_a dpsi_b   = sum_kl F_kl h^k_a h^l_b + sum_k F_k h^k_ab
#
# which is why a family written as a reparametrization of another does not
# have to difference its own cdf: the parent's closed form and the map's own
# derivatives, which the family already declares for its likelihood
# derivatives, are between them enough.
# ===========================================================================

#' The Chain Rule on a Parent's CDF Derivatives
#'
#' @description
#' Carries the derivatives of a parent's distribution function onto a new
#' parametrization, given the map's partial derivatives. A family written as a
#' map of another therefore needs no cdf derivatives of its own: the parent's
#' closed forms and the map's partials, which the family already declares for
#' its likelihood derivatives, are between them enough.
#'
#' @details
#' # The identity
#'
#' The region of integration is independent of the parameters in both
#' parametrizations, so the ordinary chain rule applies to the composition
#' \eqn{F(q; h(\psi))}:
#' \deqn{\frac{\partial F}{\partial\psi_a} = \sum_k F_k\,h^k_a, \qquad
#'       \frac{\partial^2 F}{\partial\psi_a \partial\psi_b}
#'         = \sum_{k,l} F_{kl}\,h^k_a h^l_b + \sum_k F_k\,h^k_{ab}.}
#'
#' # When it is worth taking
#'
#' The result is exact whenever the parent's own cdf derivatives are, and no
#' better than they are otherwise: applied to a parent that differences its own
#' cdf, the chain adds an exact transformation to an approximate quantity and
#' buys nothing over differencing the new cdf directly. Checking that is the
#' caller's business, and [mapped_cdf_deriv()] is where it is done.
#'
#' @section Notation:
#' \eqn{F} is the parent's distribution function, \eqn{\psi} the new
#' parameters, \eqn{h} the map from them to the parent's, \eqn{F_k} and
#' \eqn{F_{kl}} the parent's cdf derivatives and \eqn{h^k_a}, \eqn{h^k_{ab}}
#' the map's partials.
#'
#' @param parent The distribution being reparametrized.
#' @param q A numeric vector of quantiles.
#' @param th_par The parent's parameters, evaluated at the new ones.
#' @param maps The map's keyed partial tables, in the form
#'   [reparam_map_derivs()] produces: a missing key is an exact zero, so a map
#'   with many vanishing partials costs nothing for them.
#' @param new_params A character vector naming the new parameters, which names
#'   and orders the result.
#' @param order The derivative order, 1 or 2.
#'
#' @return A named list of numeric vectors, derivatives of \eqn{F} itself on
#'   the natural scale: one per new parameter at order 1, and one per
#'   [hess_names()] component at order 2.
#'
#' @seealso [mapped_cdf_deriv()], the caller that gates it;
#'   [chain_derivatives()], the same construction for likelihood derivatives;
#'   [cdf_tail_scale()].
#'
#' @keywords internal
chain_cdf_deriv <- function(parent, q, th_par, maps, new_params, order) {
  pp <- parent@params
  p <- length(pp)
  F1 <- distrib_grad_cdf(parent, q, th_par, lower.tail = TRUE, log = FALSE)
  zero <- 0 * F1[[1L]]
  hv <- function(i, tup) {
    v <- maps[[i]][[paste(sort(tup), collapse = ",")]]
    if (is.null(v)) 0 else v
  }

  if (order == 1L) {
    out <- lapply(seq_along(new_params), function(a) {
      s <- zero
      for (k in seq_len(p)) s <- s + F1[[k]] * hv(k, a)
      s
    })
    return(stats::setNames(out, new_params))
  }

  F2 <- distrib_hess_cdf(parent, q, th_par, lower.tail = TRUE, log = FALSE)
  ppairs <- hess_pairs(pp)
  tab <- new.env(parent = emptyenv())
  for (nm in names(ppairs)) {
    pr <- ppairs[[nm]]
    key <- paste(sort(c(match(pr[1], pp), match(pr[2], pp))), collapse = ",")
    assign(key, F2[[nm]], envir = tab)
  }

  npairs <- hess_pairs(new_params)
  stats::setNames(lapply(names(npairs), function(nm) {
    pr <- npairs[[nm]]
    a <- match(pr[1], new_params)
    b <- match(pr[2], new_params)
    s <- zero
    for (k in seq_len(p)) {
      for (l in seq_len(p)) {
        s <- s + get(paste(sort(c(k, l)), collapse = ","), envir = tab) *
          hv(k, a) * hv(l, b)
      }
      s <- s + F1[[k]] * hv(k, c(a, b))
    }
    s
  }), names(npairs))
}

#' CDF Derivatives of a Family Written as a Map of Another
#'
#' @description
#' The [distrib_grad_cdf()] and [distrib_hess_cdf()] body every mapped family
#' shares. It asks whether the parent's cdf derivatives are exact at every
#' order up to the one wanted, takes [chain_cdf_deriv()] if they are, and falls
#' back to [numerical_cdf_deriv()] on the new family's own cdf if they are not.
#'
#' @details
#' # The gate, and why it is a gate
#'
#' [has_exact_cdf_deriv()] is asked of the parent at order 1 and, for a
#' Hessian, at order 2 as well. The chain is exact only if what it carries is
#' exact; taken over a differenced parent it would add rounding to rounding and
#' cost more than differencing the child's cdf once. The gate is also what
#' keeps a family's page honest as the parent improves: when the inverse
#' Gaussian gained a closed second order, its second parametrization stopped
#' differencing without any edit here.
#'
#' @param distrib The mapped family, whose cdf and parameter names are read.
#' @param parent The distribution it maps onto.
#' @param th_par The parent's parameters, evaluated at the new ones.
#' @param maps The map's keyed partial tables.
#' @param q A numeric vector of quantiles.
#' @param theta A named list of the new parameters.
#' @param order The derivative order, 1 or 2.
#' @param lower.tail Is the lower tail wanted? A single logical.
#' @param log Are derivatives of the log probability wanted? A single logical.
#'
#' @return A named list of numeric vectors on the requested tail and scale: one
#'   per new parameter at order 1, and one per [hess_names()] component at
#'   order 2.
#'
#' @seealso [chain_cdf_deriv()] for the exact route;
#'   [numerical_cdf_deriv()] for the fallback;
#'   [has_exact_cdf_deriv()] for the gate.
#'
#' @keywords internal
mapped_cdf_deriv <- function(distrib, parent, th_par, maps, q, theta, order,
                             lower.tail, log) {
  exact <- has_exact_cdf_deriv(parent, 1L) &&
    (order == 1L || has_exact_cdf_deriv(parent, 2L))
  Fq <- distrib_cdf(distrib, q, theta)
  if (!exact) {
    d1 <- numerical_cdf_deriv(distrib, q, theta, 1L)
    d2 <- if (order == 1L) NULL else numerical_cdf_deriv(distrib, q, theta, 2L)
    return(cdf_tail_scale(distrib, Fq, d1, d2, lower.tail, log))
  }
  d1 <- chain_cdf_deriv(parent, q, th_par, maps, distrib@params, 1L)
  d2 <- if (order == 1L) {
    NULL
  } else {
    chain_cdf_deriv(parent, q, th_par, maps, distrib@params, 2L)
  }
  cdf_tail_scale(distrib, Fq, d1, d2, lower.tail, log)
}


# --- the reparametrize() wrapper --------------------------------------------

#' @title Log-CDF Gradient of a Reparametrized Distribution
#' @name distrib_grad_cdf.ReparamContinuousDistrib
#'
#' @description
#' The chain rule on the parent's cdf derivatives, through
#' [mapped_cdf_deriv()]. It is exact whenever the parent's are; where the
#' parent differences its own cdf, so does this, the chain having nothing
#' closed to carry.
#'
#' @details
#' The map's partials come from the object itself, `reparam_tables()` reading
#' the keyed tables the reparametrization was built with, so a family created
#' by [reparametrize()] gets closed cdf derivatives for free as soon as its
#' parent has them.
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
#' @return A named list of numeric vectors, one per new parameter, each the
#'   length of `q` recycled against `theta`.
#'
#' @seealso [distrib_hess_cdf.ReparamContinuousDistrib()] for the second order;
#'   [chain_cdf_deriv()] for the identity; [reparametrize()].
#'
#' @examples
#' # A Gaussian reparametrized in the variance. Its parent has closed cdf
#' # derivatives, so the chain is exact.
#' d <- reparametrize(
#'   gaussian1_distrib(),
#'   map = function(psi) list(mu = psi$mu, sigma = sqrt(psi$sigma2)),
#'   params = c("mu", "sigma2"),
#'   bounds = list(mu = c(-Inf, Inf), sigma2 = c(0, Inf)),
#'   links = list(mu = linkfunctions7::identity_link(),
#'                sigma2 = linkfunctions7::log_link())
#' )
#' distrib_grad_cdf(d, c(-1, 0.5, 2), list(mu = 0.3, sigma2 = 1.44))
#'
#' @keywords internal
S7::method(distrib_grad_cdf, ReparamContinuousDistrib) <-
  function(distrib, q, theta, lower.tail = TRUE, log = TRUE, ...) {
    mapped_cdf_deriv(distrib, distrib@parent_distrib,
                     reparam_theta(distrib, theta),
                     reparam_tables(distrib, theta),
                     q, theta, 1L, lower.tail, log)
  }

#' @title Log-CDF Hessian of a Reparametrized Distribution
#' @name distrib_hess_cdf.ReparamContinuousDistrib
#'
#' @description
#' The second-order chain rule on the parent's cdf derivatives, through
#' [mapped_cdf_deriv()]:
#' \eqn{\partial^2 F/\partial\psi_a\partial\psi_b
#'      = \sum_{k,l} F_{kl} h^k_a h^l_b + \sum_k F_k h^k_{ab}}. The second sum
#' is the one a first-order chain does not have, and it is why the map's second
#' partials are needed as well as its first.
#'
#' @details
#' The gate asks the parent for exactness at **both** orders. A parent with a
#' closed gradient and a differenced Hessian sends this whole method to
#' [numerical_cdf_deriv()], because carrying an exact first order onto an
#' approximate second would not make the second exact.
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
#' @return A named list of numeric vectors keyed as [hess_names()], each the
#'   length of `q` recycled against `theta`. The gradient is not returned
#'   alongside.
#'
#' @seealso [distrib_grad_cdf.ReparamContinuousDistrib()] for the first order;
#'   [chain_cdf_deriv()] for the identity; [reparametrize()].
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
#' q <- c(-1, 0.5, 2)
#' th <- list(mu = 0.3, sigma2 = 1.44)
#'
#' # Against a central difference of the reparametrized cdf.
#' exact <- distrib_hess_cdf(d, q, th, log = FALSE)
#' fd <- numerical_cdf_deriv(d, q, th, order = 2)
#' max(abs(unlist(exact[names(fd)]) - unlist(fd)))
#'
#' @keywords internal
S7::method(distrib_hess_cdf, ReparamContinuousDistrib) <-
  function(distrib, q, theta, lower.tail = TRUE, log = TRUE, ...) {
    mapped_cdf_deriv(distrib, distrib@parent_distrib,
                     reparam_theta(distrib, theta),
                     reparam_tables(distrib, theta),
                     q, theta, 2L, lower.tail, log)
  }


# --- the hand-written second parametrizations -------------------------------
#
# These are families of their own rather than reparametrize() wrappers, but
# their cdf is still the parent's at a mapped parameter, so the same chain
# applies with the map derivatives written out.

#' @title Gaussian Log-CDF Gradient in Mean and Variance
#' @name distrib_grad_cdf.Gaussian2Distrib
#'
#' @description
#' Closed form, by the chain rule on the scale parametrization's derivatives
#' through \eqn{\sigma = \sqrt{\sigma^2}}. The mean component is unchanged,
#' \eqn{-f(q)}, and the variance component is the scale one divided by
#' \eqn{2\sigma}, which is the map's only non-zero partial.
#'
#' @section Notation:
#' \eqn{\mu} is the mean, \eqn{\sigma^2 > 0} the variance, \eqn{f} the density
#' and \eqn{F} the distribution function.
#'
#' @param distrib A `Gaussian2Distrib` object, from [gaussian2_distrib()].
#' @param q A numeric vector of quantiles.
#' @param theta A named list with components `mu` and `sigma2` (positive), each
#'   a numeric vector of length 1 or `n`.
#' @param lower.tail Is the lower tail wanted? A single logical, `TRUE` by
#'   default.
#' @param log Are derivatives of the log probability wanted? A single logical,
#'   `TRUE` by default.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A named list of two numeric vectors, `mu` and `sigma2`, each the
#'   length of `q` recycled against `theta`.
#'
#' @seealso [distrib_hess_cdf.Gaussian2Distrib()] for the second order;
#'   [distrib_grad_cdf.Gaussian1Distrib()], the parent;
#'   [mapped_cdf_deriv()]; [gaussian2_distrib()].
#'
#' @examples
#' q <- c(-1, 0.5, 2)
#' g2 <- distrib_grad_cdf(gaussian2_distrib(), q,
#'                        list(mu = 0.3, sigma2 = 1.44), log = FALSE)
#' g1 <- distrib_grad_cdf(gaussian1_distrib(), q,
#'                        list(mu = 0.3, sigma = 1.2), log = FALSE)
#'
#' # The mean component is unchanged by the map.
#' all.equal(g2$mu, g1$mu)
#'
#' # The variance component is the scale one over 2 sigma.
#' all.equal(g2$sigma2, g1$sigma / (2 * 1.2))
#'
#' @keywords internal
S7::method(distrib_grad_cdf, Gaussian2Distrib) <-
  function(distrib, q, theta, lower.tail = TRUE, log = TRUE, ...) {
    mapped_cdf_deriv(distrib, gaussian1_distrib(),
                     list(mu = theta[[1]], sigma = sqrt(theta[[2]])),
                     md_gaussian2(theta), q, theta, 1L, lower.tail, log)
  }

#' @title Gaussian Log-CDF Hessian in Mean and Variance
#' @name distrib_hess_cdf.Gaussian2Distrib
#'
#' @description
#' Closed form, by the second-order chain rule on the scale parametrization's
#' derivatives through \eqn{\sigma = \sqrt{\sigma^2}}. The map's second partial
#' \eqn{\partial^2\sigma/\partial(\sigma^2)^2 = -1/(4\sigma^3)} contributes to
#' the variance-variance component, which is the term a first-order chain does
#' not have.
#'
#' @section Notation:
#' \eqn{\mu} is the mean, \eqn{\sigma^2 > 0} the variance and \eqn{F} the
#' distribution function.
#'
#' @param distrib A `Gaussian2Distrib` object, from [gaussian2_distrib()].
#' @param q A numeric vector of quantiles.
#' @param theta A named list with components `mu` and `sigma2` (positive), each
#'   a numeric vector of length 1 or `n`.
#' @param lower.tail Is the lower tail wanted? A single logical, `TRUE` by
#'   default.
#' @param log Are derivatives of the log probability wanted? A single logical,
#'   `TRUE` by default.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A named list of three numeric vectors keyed as [hess_names()], each
#'   the length of `q` recycled against `theta`.
#'
#' @seealso [distrib_grad_cdf.Gaussian2Distrib()] for the first order;
#'   [chain_cdf_deriv()] for the identity; [gaussian2_distrib()].
#'
#' @examples
#' q <- c(-1, 0.5, 2)
#' th <- list(mu = 0.3, sigma2 = 1.44)
#'
#' # Against a central difference of the cdf, which shares no arithmetic.
#' exact <- distrib_hess_cdf(gaussian2_distrib(), q, th, log = FALSE)
#' fd <- numerical_cdf_deriv(gaussian2_distrib(), q, th, order = 2)
#' max(abs(unlist(exact[names(fd)]) / unlist(fd) - 1))
#'
#' @keywords internal
S7::method(distrib_hess_cdf, Gaussian2Distrib) <-
  function(distrib, q, theta, lower.tail = TRUE, log = TRUE, ...) {
    mapped_cdf_deriv(distrib, gaussian1_distrib(),
                     list(mu = theta[[1]], sigma = sqrt(theta[[2]])),
                     md_gaussian2(theta), q, theta, 2L, lower.tail, log)
  }

#' @title Gaussian Log-CDF Gradient in Mean and Precision
#' @name distrib_grad_cdf.Gaussian3Distrib
#'
#' @description
#' Closed form, by the chain rule on the scale parametrization's derivatives
#' through \eqn{\sigma = \tau^{-1/2}}. The mean component is unchanged,
#' \eqn{-f(q)}; the precision component is the scale one times
#' \eqn{-\tfrac12\tau^{-3/2}}, so it carries the opposite sign, a larger
#' precision being a tighter distribution.
#'
#' @section Notation:
#' \eqn{\mu} is the mean, \eqn{\tau > 0} the precision, \eqn{f} the density and
#' \eqn{F} the distribution function.
#'
#' @param distrib A `Gaussian3Distrib` object, from [gaussian3_distrib()].
#' @param q A numeric vector of quantiles.
#' @param theta A named list with components `mu` and `tau` (positive), each a
#'   numeric vector of length 1 or `n`.
#' @param lower.tail Is the lower tail wanted? A single logical, `TRUE` by
#'   default.
#' @param log Are derivatives of the log probability wanted? A single logical,
#'   `TRUE` by default.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A named list of two numeric vectors, `mu` and `tau`, each the length
#'   of `q` recycled against `theta`.
#'
#' @seealso [distrib_hess_cdf.Gaussian3Distrib()] for the second order;
#'   [distrib_grad_cdf.Gaussian1Distrib()], the parent; [gaussian3_distrib()].
#'
#' @examples
#' q <- c(-1, 0.5, 2)
#' g3 <- distrib_grad_cdf(gaussian3_distrib(), q,
#'                        list(mu = 0.3, tau = 1 / 1.44), log = FALSE)
#' g1 <- distrib_grad_cdf(gaussian1_distrib(), q,
#'                        list(mu = 0.3, sigma = 1.2), log = FALSE)
#'
#' # The precision component is the scale one times -0.5 tau^(-3/2).
#' all.equal(g3$tau, g1$sigma * (-0.5) * (1 / 1.44)^(-1.5))
#'
#' @keywords internal
S7::method(distrib_grad_cdf, Gaussian3Distrib) <-
  function(distrib, q, theta, lower.tail = TRUE, log = TRUE, ...) {
    mapped_cdf_deriv(distrib, gaussian1_distrib(),
                     list(mu = theta[[1]], sigma = 1 / sqrt(theta[[2]])),
                     md_gaussian3(theta), q, theta, 1L, lower.tail, log)
  }

#' @title Gaussian Log-CDF Hessian in Mean and Precision
#' @name distrib_hess_cdf.Gaussian3Distrib
#'
#' @description
#' Closed form, by the second-order chain rule through
#' \eqn{\sigma = \tau^{-1/2}}. The map's second partial
#' \eqn{\partial^2\sigma/\partial\tau^2 = \tfrac34\tau^{-5/2}} contributes to
#' the precision-precision component.
#'
#' @section Notation:
#' \eqn{\mu} is the mean, \eqn{\tau > 0} the precision and \eqn{F} the
#' distribution function.
#'
#' @param distrib A `Gaussian3Distrib` object, from [gaussian3_distrib()].
#' @param q A numeric vector of quantiles.
#' @param theta A named list with components `mu` and `tau` (positive), each a
#'   numeric vector of length 1 or `n`.
#' @param lower.tail Is the lower tail wanted? A single logical, `TRUE` by
#'   default.
#' @param log Are derivatives of the log probability wanted? A single logical,
#'   `TRUE` by default.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A named list of three numeric vectors keyed as [hess_names()], each
#'   the length of `q` recycled against `theta`.
#'
#' @seealso [distrib_grad_cdf.Gaussian3Distrib()] for the first order;
#'   [chain_cdf_deriv()] for the identity; [gaussian3_distrib()].
#'
#' @examples
#' q <- c(-1, 0.5, 2)
#' th <- list(mu = 0.3, tau = 1 / 1.44)
#'
#' # Against a central difference of the cdf, which shares no arithmetic.
#' exact <- distrib_hess_cdf(gaussian3_distrib(), q, th, log = FALSE)
#' fd <- numerical_cdf_deriv(gaussian3_distrib(), q, th, order = 2)
#' max(abs(unlist(exact[names(fd)]) / unlist(fd) - 1))
#'
#' @keywords internal
S7::method(distrib_hess_cdf, Gaussian3Distrib) <-
  function(distrib, q, theta, lower.tail = TRUE, log = TRUE, ...) {
    mapped_cdf_deriv(distrib, gaussian1_distrib(),
                     list(mu = theta[[1]], sigma = 1 / sqrt(theta[[2]])),
                     md_gaussian3(theta), q, theta, 2L, lower.tail, log)
  }

#' @title Inverse Gaussian Log-CDF Gradient in Mean and Rate
#' @name distrib_grad_cdf.InvGauss2Distrib
#'
#' @description
#' Closed form, by the chain rule on the dispersion parametrization's gradient
#' through \eqn{\phi = 1/\lambda}. The mean component is unchanged; the rate
#' component is the dispersion one times \eqn{-1/\lambda^2}.
#'
#' @details
#' The inverse Gaussian is unusual among positive families in having an
#' elementary distribution function,
#' \eqn{F(y) = \Phi(a) + e^{2/(\phi\mu)}\,\Phi(b)}, so the parent's derivatives
#' are closed and the chain carries them.
#'
#' The second and higher orders are registered elsewhere, by
#' `register_mapped_cdf_k()` in `cdf_mapped_higher.R`, at orders 2 to 4. That
#' file's registration supersedes anything this one would give: the Hessian was
#' left on the fallback while the parent differenced its own second order, and
#' moved onto the chain when the parent stopped.
#'
#' @section Notation:
#' \eqn{\mu > 0} is the mean, \eqn{\lambda > 0} the rate, \eqn{\phi = 1/\lambda}
#' the dispersion and \eqn{\Phi} the standard normal distribution function.
#'
#' @param distrib An `InvGauss2Distrib` object, from [invgauss2_distrib()].
#' @param q A numeric vector of quantiles, positive.
#' @param theta A named list with components `mu` (positive) and `lambda`
#'   (positive), each a numeric vector of length 1 or `n`.
#' @param lower.tail Is the lower tail wanted? A single logical, `TRUE` by
#'   default.
#' @param log Are derivatives of the log probability wanted? A single logical,
#'   `TRUE` by default.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A named list of two numeric vectors, `mu` and `lambda`, each the
#'   length of `q` recycled against `theta`.
#'
#' @seealso [distrib_grad_cdf()] for the generic;
#'   [mapped_cdf_deriv()] and [chain_cdf_deriv()];
#'   [invgauss2_distrib()].
#'
#' @examples
#' d <- invgauss2_distrib()
#' q <- c(0.5, 2, 5)
#' th <- list(mu = 2, lambda = 8)
#'
#' # Against a central difference of the cdf, which shares no arithmetic.
#' fd <- numerical_cdf_deriv(d, q, th, order = 1)
#' max(abs(unlist(distrib_grad_cdf(d, q, th, log = FALSE)) / unlist(fd) - 1))
#'
#' @keywords internal
S7::method(distrib_grad_cdf, InvGauss2Distrib) <-
  function(distrib, q, theta, lower.tail = TRUE, log = TRUE, ...) {
    mapped_cdf_deriv(distrib, invgauss1_distrib(),
                     list(mu = theta[[1]], phi = 1 / theta[[2]]),
                     md_invgauss2(theta), q, theta, 1L, lower.tail, log)
  }


# --- families that are location-scale after all -----------------------------

#' @title Gumbel Log-CDF Gradient
#' @name distrib_grad_cdf.GumbelDistrib
#'
#' @description
#' Closed form from the location-scale structure, as for the Gaussian:
#' \eqn{\partial F/\partial\mu = -f(q)} and
#' \eqn{\partial F/\partial\sigma = -z f(q)} with \eqn{z = (q-\mu)/\sigma}. The
#' method is [loc_scale_grad_cdf()] itself.
#'
#' @details
#' The family's distribution function is \eqn{F(q) = \exp(-e^{-z})}, so nothing
#' here needs a quadrature or a series, which is why a Gumbel is usable for
#' censored extreme-value data: the score of a right-censored observation is
#' \eqn{f(q)/S(q)} and comes from the density alone.
#'
#' @section Notation:
#' \eqn{\mu} is the location, \eqn{\sigma > 0} the scale,
#' \eqn{z = (q-\mu)/\sigma}, \eqn{f} the density and \eqn{F} the distribution
#' function. The mean is \eqn{\mu + \gamma\sigma} and not \eqn{\mu}.
#'
#' @param distrib A `GumbelDistrib` object, from [gumbel_distrib()].
#' @param q A numeric vector of quantiles.
#' @param theta A named list with components `mu` and `sigma` (positive), each
#'   a numeric vector of length 1 or `n`.
#' @param lower.tail Is the lower tail wanted? A single logical, `TRUE` by
#'   default.
#' @param log Are derivatives of the log probability wanted? A single logical,
#'   `TRUE` by default.
#'
#' @return A named list of two numeric vectors, `mu` and `sigma`, each the
#'   length of `q` recycled against `theta`.
#'
#' @seealso [loc_scale_grad_cdf()] for the shared body;
#'   [distrib_hess_cdf.GumbelDistrib()] for the second order;
#'   [gumbel_distrib()].
#'
#' @examples
#' d <- gumbel_distrib()
#' th <- list(mu = 0.3, sigma = 1.2)
#' q <- c(-1, 0.5, 2)
#'
#' all.equal(distrib_grad_cdf(d, q, th, log = FALSE)$mu,
#'           -distrib_pdf(d, q, th))
#'
#' @keywords internal
S7::method(distrib_grad_cdf, GumbelDistrib) <- loc_scale_grad_cdf

#' @title Gumbel Log-CDF Hessian
#' @name distrib_hess_cdf.GumbelDistrib
#'
#' @description
#' Closed form from the same location-scale structure, through
#' [loc_scale_hess_cdf()]. The response derivative it reads is
#' \eqn{\ell_y = (e^{-z} - 1)/\sigma}, which is smooth everywhere, so all three
#' components are continuous on the whole line.
#'
#' @section Notation:
#' \eqn{\mu} is the location, \eqn{\sigma > 0} the scale,
#' \eqn{z = (q-\mu)/\sigma}, \eqn{f} the density and
#' \eqn{\ell_y = \partial\log f/\partial y}.
#'
#' @param distrib A `GumbelDistrib` object, from [gumbel_distrib()].
#' @param q A numeric vector of quantiles.
#' @param theta A named list with components `mu` and `sigma` (positive), each
#'   a numeric vector of length 1 or `n`.
#' @param lower.tail Is the lower tail wanted? A single logical, `TRUE` by
#'   default.
#' @param log Are derivatives of the log probability wanted? A single logical,
#'   `TRUE` by default.
#'
#' @return A named list of three numeric vectors keyed as [hess_names()], each
#'   the length of `q` recycled against `theta`.
#'
#' @seealso [loc_scale_hess_cdf()] for the shared body;
#'   [distrib_grad_cdf.GumbelDistrib()] for the first order;
#'   [gumbel_distrib()].
#'
#' @examples
#' d <- gumbel_distrib()
#' th <- list(mu = 0.3, sigma = 1.2)
#' q <- c(-1, 0.5, 2)
#'
#' # Against a central difference of the cdf, which shares no arithmetic.
#' exact <- distrib_hess_cdf(d, q, th, log = FALSE)
#' fd <- numerical_cdf_deriv(d, q, th, order = 2)
#' max(abs(unlist(exact[names(fd)]) / unlist(fd) - 1))
#'
#' @keywords internal
S7::method(distrib_hess_cdf, GumbelDistrib) <- loc_scale_hess_cdf


# --- families location-scale in their first two parameters ------------------

#' CDF Hessian When Only Some Parameters Are Location-Scale
#'
#' @description
#' The second-order companion of [partial_loc_scale_grad_cdf()], and the
#' [distrib_hess_cdf()] body the Student t, the pseudo-Huber and the skew t
#' share. The three components in the location and the scale come from
#' [loc_scale_cdf_deriv()]'s formulas, and every component touching a shape
#' parameter is differenced.
#'
#' @details
#' The gradient is assembled here as well, on the same split, because the
#' log-scale conversion in [cdf_tail_scale()] reads it. Both differencing calls
#' pass a `which` argument, so no component that has a closed form costs a cdf
#' evaluation.
#'
#' For a family with \eqn{p} parameters of which two are the location and the
#' scale, the closed block is 3 components of the \eqn{p(p+1)/2}: 3 of 6 for
#' the Student t and the pseudo-Huber, 3 of 10 for the skew t.
#'
#' @section Notation:
#' \eqn{\mu} is the location, \eqn{\sigma > 0} the scale,
#' \eqn{z = (q-\mu)/\sigma}, \eqn{f} the density and
#' \eqn{\ell_y = \partial\log f/\partial y}.
#'
#' @param distrib An object inheriting from `distrib` whose first two
#'   parameters are a location and a scale, and which has at least one more.
#' @param q A numeric vector of quantiles.
#' @param theta A named list of parameters, the location first and the scale
#'   second.
#' @param lower.tail Is the lower tail wanted? A single logical, `TRUE` by
#'   default.
#' @param log Are derivatives of the log probability wanted? A single logical,
#'   `TRUE` by default.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A named list of numeric vectors keyed as [hess_names()], in that
#'   enumeration's order.
#'
#' @seealso [partial_loc_scale_grad_cdf()] for the first order;
#'   [loc_scale_cdf_deriv()] for the closed formulas;
#'   [numerical_cdf_deriv()] for the differenced components.
#'
#' @keywords internal
partial_loc_scale_hess_cdf <- function(distrib, q, theta, lower.tail = TRUE,
                                       log = TRUE, ...) {
  params <- distrib@params
  s <- theta[[2]]
  z <- (q - theta[[1]]) / s
  f <- distrib_pdf(distrib, q, theta)
  ly <- distrib_grad_y(distrib, q, theta)

  closed <- stats::setNames(
    list(f * ly, f * (z^2 * ly + 2 * z / s), f * (z * ly + 1 / s)),
    hess_names(params[1:2])
  )
  nms <- hess_names(params)
  rest <- setdiff(nms, names(closed))
  d2 <- c(closed, numerical_cdf_deriv(distrib, q, theta, 2L, which = rest))[nms]

  rest1 <- params[-(1:2)]
  d1 <- c(list(-f, -z * f),
          numerical_cdf_deriv(distrib, q, theta, 1L, which = rest1))
  names(d1) <- params
  cdf_tail_scale(distrib, distrib_cdf(distrib, q, theta), d1, d2,
                 lower.tail, log)
}

#' @title Student t Log-CDF Hessian
#' @name distrib_hess_cdf.StudentT1Distrib
#'
#' @description
#' Closed form in the location-scale block, `mu_mu`, `sigma_sigma` and
#' `mu_sigma`; the three components touching the degrees of freedom are
#' differenced, that derivative having no elementary form. The method is
#' [partial_loc_scale_hess_cdf()] itself, shared with the pseudo-Huber and the
#' skew t.
#'
#' @section Notation:
#' \eqn{\mu} is the location, \eqn{\sigma > 0} the scale, \eqn{\nu > 0} the
#' degrees of freedom, \eqn{z = (q-\mu)/\sigma} and \eqn{f} the density.
#'
#' @param distrib A `StudentT1Distrib` object, from [student_t1_distrib()].
#' @param q A numeric vector of quantiles.
#' @param theta A named list with components `mu`, `sigma` (positive) and `nu`
#'   (positive), each a numeric vector of length 1 or `n`.
#' @param lower.tail Is the lower tail wanted? A single logical, `TRUE` by
#'   default.
#' @param log Are derivatives of the log probability wanted? A single logical,
#'   `TRUE` by default.
#'
#' @return A named list of six numeric vectors keyed as [hess_names()], each
#'   the length of `q` recycled against `theta`.
#'
#' @seealso [partial_loc_scale_hess_cdf()] for the shared body;
#'   [distrib_grad_cdf.StudentT1Distrib()] for the first order;
#'   [student_t1_distrib()].
#'
#' @examples
#' d <- student_t1_distrib()
#' th <- list(mu = 0.3, sigma = 1.2, nu = 6)
#' q <- c(-1, 0.5, 2)
#'
#' # Six components: three closed, three differenced.
#' names(distrib_hess_cdf(d, q, th))
#'
#' @keywords internal
S7::method(distrib_hess_cdf, StudentT1Distrib) <- partial_loc_scale_hess_cdf

#' @title Pseudo-Huber Log-CDF Hessian
#' @name distrib_hess_cdf.PseudoHuberDistrib
#'
#' @description
#' Closed form in the location-scale block; the three components touching the
#' shape are differenced. The method is [partial_loc_scale_hess_cdf()] itself,
#' shared with the Student t and the skew t.
#'
#' @details
#' The saving is larger here than for the Student t, this family's distribution
#' function being a quadrature: the closed block reads the density and its
#' response derivative, where differencing it would run the quadrature four
#' times per component.
#'
#' @section Notation:
#' \eqn{\mu} is the location, \eqn{\sigma > 0} the scale, \eqn{\nu > 0} the
#' shape, \eqn{z = (q-\mu)/\sigma} and \eqn{f} the density.
#'
#' @param distrib A `PseudoHuberDistrib` object, from [pseudohuber_distrib()].
#' @param q A numeric vector of quantiles.
#' @param theta A named list with components `mu`, `sigma` (positive) and `nu`
#'   (positive), each a numeric vector of length 1 or `n`.
#' @param lower.tail Is the lower tail wanted? A single logical, `TRUE` by
#'   default.
#' @param log Are derivatives of the log probability wanted? A single logical,
#'   `TRUE` by default.
#'
#' @return A named list of six numeric vectors keyed as [hess_names()], each
#'   the length of `q` recycled against `theta`.
#'
#' @seealso [partial_loc_scale_hess_cdf()] for the shared body;
#'   [distrib_grad_cdf.PseudoHuberDistrib()] for the first order;
#'   [pseudohuber_distrib()].
#'
#' @examples
#' d <- pseudohuber_distrib()
#' th <- list(mu = 0.3, sigma = 1.2, nu = 4)
#'
#' distrib_hess_cdf(d, c(-1, 2), th)$mu_mu
#'
#' @keywords internal
S7::method(distrib_hess_cdf, PseudoHuberDistrib) <- partial_loc_scale_hess_cdf

#' @title Skew t Log-CDF Gradient
#' @name distrib_grad_cdf.SkewTDistrib
#'
#' @description
#' Closed form in the location and the scale, \eqn{-f(q)} and \eqn{-z f(q)}
#' with \eqn{z = (q-\mu)/\sigma}; the shape and the degrees of freedom are
#' differenced. The method is [partial_loc_scale_grad_cdf()] itself, shared
#' with the Student t and the pseudo-Huber.
#'
#' @details
#' Two of the four components are closed. The shape and the degrees of freedom
#' enter the distribution function through a Student t distribution function at
#' \eqn{\nu+1} degrees of freedom, whose derivatives in either have no
#' elementary form, so both are differenced.
#'
#' @section Notation:
#' \eqn{\mu} is the location, \eqn{\sigma > 0} the scale, \eqn{\alpha} the
#' shape, \eqn{\nu > 0} the degrees of freedom, \eqn{z = (q-\mu)/\sigma} and
#' \eqn{f} the density.
#'
#' @param distrib A `SkewTDistrib` object, from [skewt_distrib()].
#' @param q A numeric vector of quantiles.
#' @param theta A named list with components `mu`, `sigma` (positive), `alpha`
#'   (any sign) and `nu` (positive), each a numeric vector of length 1 or `n`.
#' @param lower.tail Is the lower tail wanted? A single logical, `TRUE` by
#'   default.
#' @param log Are derivatives of the log probability wanted? A single logical,
#'   `TRUE` by default.
#'
#' @return A named list of four numeric vectors, `mu`, `sigma`, `alpha` and
#'   `nu`, each the length of `q` recycled against `theta`.
#'
#' @seealso [partial_loc_scale_grad_cdf()] for the shared body;
#'   [distrib_hess_cdf.SkewTDistrib()] for the second order;
#'   [skewt_distrib()].
#'
#' @examples
#' d <- skewt_distrib()
#' th <- list(mu = 0.3, sigma = 1.2, alpha = 2, nu = 6)
#' q <- c(-1, 0.5, 2)
#'
#' # The location component is exact, the density itself.
#' all.equal(distrib_grad_cdf(d, q, th, log = FALSE)$mu,
#'           -distrib_pdf(d, q, th))
#'
#' @keywords internal
S7::method(distrib_grad_cdf, SkewTDistrib) <- partial_loc_scale_grad_cdf

#' @title Skew t Log-CDF Hessian
#' @name distrib_hess_cdf.SkewTDistrib
#'
#' @description
#' Closed form in the location-scale block, three of the ten components; the
#' seven touching the shape or the degrees of freedom are differenced. The
#' method is [partial_loc_scale_hess_cdf()] itself, shared with the Student t
#' and the pseudo-Huber.
#'
#' @section Notation:
#' \eqn{\mu} is the location, \eqn{\sigma > 0} the scale, \eqn{\alpha} the
#' shape, \eqn{\nu > 0} the degrees of freedom, \eqn{z = (q-\mu)/\sigma} and
#' \eqn{f} the density.
#'
#' @param distrib A `SkewTDistrib` object, from [skewt_distrib()].
#' @param q A numeric vector of quantiles.
#' @param theta A named list with components `mu`, `sigma` (positive), `alpha`
#'   (any sign) and `nu` (positive), each a numeric vector of length 1 or `n`.
#' @param lower.tail Is the lower tail wanted? A single logical, `TRUE` by
#'   default.
#' @param log Are derivatives of the log probability wanted? A single logical,
#'   `TRUE` by default.
#'
#' @return A named list of ten numeric vectors keyed as [hess_names()], each
#'   the length of `q` recycled against `theta`.
#'
#' @seealso [partial_loc_scale_hess_cdf()] for the shared body;
#'   [distrib_grad_cdf.SkewTDistrib()] for the first order;
#'   [skewt_distrib()].
#'
#' @examples
#' d <- skewt_distrib()
#' th <- list(mu = 0.3, sigma = 1.2, alpha = 2, nu = 6)
#'
#' # Ten components, of which mu_mu, sigma_sigma and mu_sigma are closed.
#' names(distrib_hess_cdf(d, c(-1, 2), th))
#'
#' @keywords internal
S7::method(distrib_hess_cdf, SkewTDistrib) <- partial_loc_scale_hess_cdf


# --- positive families with an elementary distribution function -------------

#' @title Lognormal Log-CDF Hessian
#' @name distrib_hess_cdf.Lognormal1Distrib
#'
#' @description
#' Closed form. On the log scale the family is location-scale, so with
#' \eqn{z = (\log q - \mu)/\sigma} and \eqn{\varphi} the standard normal
#' density,
#' \deqn{\frac{\partial^2 F}{\partial\mu^2} = -\frac{z\varphi}{\sigma^2},
#'       \qquad
#'       \frac{\partial^2 F}{\partial\mu\,\partial\sigma^2}
#'         = \frac{\varphi(1-z^2)}{2\sigma^3},
#'       \qquad
#'       \frac{\partial^2 F}{\partial(\sigma^2)^2}
#'         = \frac{\varphi z(3-z^2)}{4\sigma^4}.}
#'
#' @details
#' The formulas are written in \eqn{\varphi(z)}, so the Jacobian factor \eqn{q}
#' of the gradient does not appear: at second order it would have to be
#' differentiated too, and expressing everything on the log scale avoids that.
#'
#' The three components vanish at the points where their polynomial in \eqn{z}
#' does, so `mu_mu` is exactly zero at the median \eqn{q = e^{\mu}}. A relative
#' comparison against a numerical derivative is meaningless there and an
#' absolute one is what to use.
#'
#' @section Notation:
#' \eqn{\mu} is the mean of \eqn{\log Y}, \eqn{\sigma^2 > 0} its variance,
#' \eqn{z = (\log q - \mu)/\sigma}, \eqn{\varphi} the standard normal density
#' and \eqn{F} the distribution function of \eqn{Y}.
#'
#' @param distrib A `Lognormal1Distrib` object, from [lognormal1_distrib()].
#' @param q A numeric vector of quantiles. At or below zero every component is
#'   zero, the distribution function being flat there.
#' @param theta A named list with components `mu` (any real value) and `sigma2`
#'   (positive), each a numeric vector of length 1 or `n`.
#' @param lower.tail Is the lower tail wanted? A single logical, `TRUE` by
#'   default.
#' @param log Are derivatives of the log probability wanted? A single logical,
#'   `TRUE` by default.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A named list of three numeric vectors keyed as [hess_names()], each
#'   the length of `q` recycled against `theta`.
#'
#' @seealso [distrib_grad_cdf.Lognormal1Distrib()] for the first order;
#'   [distrib_hess_cdf.Gaussian1Distrib()], the family on the log scale;
#'   [lognormal1_distrib()].
#'
#' @examples
#' d <- lognormal1_distrib()
#' th <- list(mu = 0, sigma2 = 1)
#' q <- c(0.5, 1, 3)
#'
#' # Against a central difference of the cdf, on an absolute scale: mu_mu is
#' # exactly zero at the median q = exp(mu) = 1.
#' exact <- distrib_hess_cdf(d, q, th, log = FALSE)
#' fd <- numerical_cdf_deriv(d, q, th, order = 2)
#' max(abs(unlist(exact[names(fd)]) - unlist(fd)))
#'
#' @keywords internal
S7::method(distrib_hess_cdf, Lognormal1Distrib) <-
  function(distrib, q, theta, lower.tail = TRUE, log = TRUE, ...) {
    v <- theta[[2]]
    s <- sqrt(v)
    pos <- q > 0
    z <- ifelse(pos, (base::log(ifelse(pos, q, 1)) - theta[[1]]) / s, 0)
    ph <- stats::dnorm(z) * pos
    d1 <- list(mu = -ph / s, sigma2 = -ph * z / (2 * v))
    d2 <- list(mu_mu = -z * ph / v,
               sigma2_sigma2 = ph * z * (3 - z^2) / (4 * v^2),
               mu_sigma2 = ph * (1 - z^2) / (2 * s^3))
    names(d2) <- hess_names(distrib@params)
    cdf_tail_scale(distrib, distrib_cdf(distrib, q, theta), d1, d2,
                   lower.tail, log)
  }
