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
#' Carries the derivatives of the parent's distribution function onto a new
#' parametrization, given the map's partial derivatives.
#'
#' @details
#' The map derivatives arrive in the keyed form
#' \code{\link{reparam_map_derivs}} produces, so a missing key is an exact
#' zero, and the result is exact whenever the parent's own cdf derivatives
#' are. It is the caller's business to check that they are: applied to a
#' parent whose derivatives are themselves differenced, the chain adds an
#' exact transformation to an approximate quantity and buys nothing over
#' differencing the new cdf directly.
#'
#' @param parent The distribution being reparametrized.
#' @param q A numeric vector of quantiles.
#' @param th_par The parent's parameters at the new ones.
#' @param maps The map's keyed partial tables.
#' @param new_params The new parameter names.
#' @param order The derivative order, 1 or 2.
#'
#' @return A named list of derivative components of \eqn{F}.
#'
#' @seealso \code{\link{cdf_tail_scale}}, \code{\link{chain_derivatives}}
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

#' CDF Derivatives Through a Map, With the Fallback as a Guard
#'
#' @description
#' The body shared by the families that obtain their cdf derivatives from a
#' parent's: the chain rule of \code{\link{chain_cdf_deriv}} when the parent
#' has a closed form at that order, and the finite-difference fallback
#' otherwise.
#'
#' @param distrib The distribution in the new parametrization.
#' @param parent The parent distribution.
#' @param th_par The parent's parameters at the new ones.
#' @param maps The map's keyed partial tables.
#' @param q A numeric vector of quantiles.
#' @param theta A named list of the new parameters.
#' @param order The derivative order, 1 or 2.
#' @param lower.tail Logical; whether the lower tail is wanted.
#' @param log Logical; whether derivatives of the log probability are wanted.
#'
#' @return A named list of derivative component vectors.
#'
#' @seealso \code{\link{chain_cdf_deriv}}
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
#' @description
#' The chain rule on the parent's cdf derivatives, which is exact whenever the
#' parent's are; when the parent differences its own cdf, so does this, the
#' chain having nothing closed to carry.
#' @param distrib A reparametrized distribution.
#' @param q A numeric vector of quantiles.
#' @param theta A named list of the new parameters.
#' @param lower.tail Logical; if \code{TRUE} (default), the lower tail.
#' @param log Logical; if \code{TRUE} (default), derivatives of the log probability.
#' @return A named list, one vector per parameter.
#' @seealso \code{\link{reparametrize}}
#' @keywords internal
S7::method(distrib_grad_cdf, ReparamContinuousDistrib) <-
  function(distrib, q, theta, lower.tail = TRUE, log = TRUE) {
    mapped_cdf_deriv(distrib, distrib@parent_distrib,
                     reparam_theta(distrib, theta),
                     reparam_tables(distrib, theta),
                     q, theta, 1L, lower.tail, log)
  }

#' @rdname distrib_grad_cdf.ReparamContinuousDistrib
#' @name distrib_hess_cdf.ReparamContinuousDistrib
#' @keywords internal
S7::method(distrib_hess_cdf, ReparamContinuousDistrib) <-
  function(distrib, q, theta, lower.tail = TRUE, log = TRUE) {
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

#' @title Gaussian Log-CDF Derivatives in Mean and Variance
#' @name distrib_grad_cdf.Gaussian2Distrib
#' @description
#' Closed form, by the chain rule on the scale parametrization's derivatives
#' through \eqn{\sigma = \sqrt{\sigma^2}}.
#' @param distrib A \code{Gaussian2Distrib} object.
#' @param q A numeric vector of quantiles.
#' @param theta A list containing \code{mu} and \code{sigma2}.
#' @param lower.tail Logical; if \code{TRUE} (default), the lower tail.
#' @param log Logical; if \code{TRUE} (default), derivatives of the log probability.
#' @return A named list, one vector per parameter.
#' @seealso \code{\link{gaussian2_distrib}}
#' @keywords internal
S7::method(distrib_grad_cdf, Gaussian2Distrib) <-
  function(distrib, q, theta, lower.tail = TRUE, log = TRUE) {
    mapped_cdf_deriv(distrib, gaussian1_distrib(),
                     list(mu = theta[[1]], sigma = sqrt(theta[[2]])),
                     md_gaussian2(theta), q, theta, 1L, lower.tail, log)
  }

#' @rdname distrib_grad_cdf.Gaussian2Distrib
#' @name distrib_hess_cdf.Gaussian2Distrib
#' @keywords internal
S7::method(distrib_hess_cdf, Gaussian2Distrib) <-
  function(distrib, q, theta, lower.tail = TRUE, log = TRUE) {
    mapped_cdf_deriv(distrib, gaussian1_distrib(),
                     list(mu = theta[[1]], sigma = sqrt(theta[[2]])),
                     md_gaussian2(theta), q, theta, 2L, lower.tail, log)
  }

#' @title Gaussian Log-CDF Derivatives in Mean and Precision
#' @name distrib_grad_cdf.Gaussian3Distrib
#' @description
#' Closed form, by the chain rule on the scale parametrization's derivatives
#' through \eqn{\sigma = \tau^{-1/2}}.
#' @param distrib A \code{Gaussian3Distrib} object.
#' @param q A numeric vector of quantiles.
#' @param theta A list containing \code{mu} and \code{tau}.
#' @param lower.tail Logical; if \code{TRUE} (default), the lower tail.
#' @param log Logical; if \code{TRUE} (default), derivatives of the log probability.
#' @return A named list, one vector per parameter.
#' @seealso \code{\link{gaussian3_distrib}}
#' @keywords internal
S7::method(distrib_grad_cdf, Gaussian3Distrib) <-
  function(distrib, q, theta, lower.tail = TRUE, log = TRUE) {
    mapped_cdf_deriv(distrib, gaussian1_distrib(),
                     list(mu = theta[[1]], sigma = 1 / sqrt(theta[[2]])),
                     md_gaussian3(theta), q, theta, 1L, lower.tail, log)
  }

#' @rdname distrib_grad_cdf.Gaussian3Distrib
#' @name distrib_hess_cdf.Gaussian3Distrib
#' @keywords internal
S7::method(distrib_hess_cdf, Gaussian3Distrib) <-
  function(distrib, q, theta, lower.tail = TRUE, log = TRUE) {
    mapped_cdf_deriv(distrib, gaussian1_distrib(),
                     list(mu = theta[[1]], sigma = 1 / sqrt(theta[[2]])),
                     md_gaussian3(theta), q, theta, 2L, lower.tail, log)
  }

#' @title Inverse Gaussian Log-CDF Gradient in Mean and Shape
#' @name distrib_grad_cdf.InvGauss2Distrib
#' @description
#' Closed form, by the chain rule on the dispersion parametrization's
#' gradient through \eqn{\phi = 1/\lambda}. The second order stays on the
#' fallback, the parent having no closed form there.
#' @param distrib An \code{InvGauss2Distrib} object.
#' @param q A numeric vector of quantiles.
#' @param theta A list containing \code{mu} and \code{lambda}.
#' @param lower.tail Logical; if \code{TRUE} (default), the lower tail.
#' @param log Logical; if \code{TRUE} (default), derivatives of the log probability.
#' @return A named list, one vector per parameter.
#' @seealso \code{\link{invgauss2_distrib}}
#' @keywords internal
S7::method(distrib_grad_cdf, InvGauss2Distrib) <-
  function(distrib, q, theta, lower.tail = TRUE, log = TRUE) {
    mapped_cdf_deriv(distrib, invgauss1_distrib(),
                     list(mu = theta[[1]], phi = 1 / theta[[2]]),
                     md_invgauss2(theta), q, theta, 1L, lower.tail, log)
  }


# --- families that are location-scale after all -----------------------------

#' @title Gumbel Log-CDF Derivatives
#' @name distrib_grad_cdf.GumbelDistrib
#' @description
#' Closed form from the location-scale structure, as for the Gaussian:
#' \eqn{\partial F/\partial\mu = -f} and \eqn{\partial F/\partial\sigma = -zf}.
#' @param distrib A \code{GumbelDistrib} object.
#' @param q A numeric vector of quantiles.
#' @param theta A list containing \code{mu} and \code{sigma}.
#' @param lower.tail Logical; if \code{TRUE} (default), the lower tail.
#' @param log Logical; if \code{TRUE} (default), derivatives of the log probability.
#' @return A named list, one vector per parameter.
#' @seealso \code{\link{gumbel_distrib}}
#' @keywords internal
S7::method(distrib_grad_cdf, GumbelDistrib) <- loc_scale_grad_cdf

#' @rdname distrib_grad_cdf.GumbelDistrib
#' @name distrib_hess_cdf.GumbelDistrib
#' @keywords internal
S7::method(distrib_hess_cdf, GumbelDistrib) <- loc_scale_hess_cdf


# --- families location-scale in their first two parameters ------------------

#' CDF Hessian When Only Some Parameters Are Location-Scale
#'
#' @description
#' The second-order companion of \code{\link{partial_loc_scale_grad_cdf}}: the
#' three components in the location and scale in closed form, the components
#' involving a shape parameter by finite differences.
#'
#' @param distrib An object inheriting from class \code{"distrib"}.
#' @param q A numeric vector of quantiles.
#' @param theta A named list of parameters.
#' @param lower.tail Logical; whether the lower tail is wanted.
#' @param log Logical; whether derivatives of the log probability are wanted.
#'
#' @return A named list of Hessian component vectors.
#'
#' @seealso \code{\link{loc_scale_cdf_deriv}}
#' @keywords internal
partial_loc_scale_hess_cdf <- function(distrib, q, theta, lower.tail = TRUE,
                                       log = TRUE) {
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
#' @description
#' Closed form in the location and scale block; the components involving the
#' degrees of freedom are differenced, having no elementary form.
#' @param distrib A \code{StudentT1Distrib} object.
#' @param q A numeric vector of quantiles.
#' @param theta A list containing \code{mu}, \code{sigma} and \code{nu}.
#' @param lower.tail Logical; if \code{TRUE} (default), the lower tail.
#' @param log Logical; if \code{TRUE} (default), derivatives of the log probability.
#' @return A named list keyed as \code{\link{hess_names}}.
#' @seealso \code{\link{student_t1_distrib}}
#' @keywords internal
S7::method(distrib_hess_cdf, StudentT1Distrib) <- partial_loc_scale_hess_cdf

#' @title Pseudo-Huber Log-CDF Hessian
#' @name distrib_hess_cdf.PseudoHuberDistrib
#' @description
#' Closed form in the location and scale block; the shape is differenced.
#' Differencing the cdf itself would be poor here, it being a quadrature.
#' @param distrib A \code{PseudoHuberDistrib} object.
#' @param q A numeric vector of quantiles.
#' @param theta A list containing \code{mu}, \code{sigma} and \code{nu}.
#' @param lower.tail Logical; if \code{TRUE} (default), the lower tail.
#' @param log Logical; if \code{TRUE} (default), derivatives of the log probability.
#' @return A named list keyed as \code{\link{hess_names}}.
#' @seealso \code{\link{pseudohuber_distrib}}
#' @keywords internal
S7::method(distrib_hess_cdf, PseudoHuberDistrib) <- partial_loc_scale_hess_cdf

#' @title Skew t Log-CDF Derivatives
#' @name distrib_grad_cdf.SkewTDistrib
#' @description
#' Closed form in the location and scale; the shape and the degrees of freedom
#' are differenced.
#' @param distrib A \code{SkewTDistrib} object.
#' @param q A numeric vector of quantiles.
#' @param theta A list containing \code{mu}, \code{sigma}, \code{alpha} and
#'   \code{nu}.
#' @param lower.tail Logical; if \code{TRUE} (default), the lower tail.
#' @param log Logical; if \code{TRUE} (default), derivatives of the log probability.
#' @return A named list, one vector per parameter.
#' @seealso \code{\link{skewt_distrib}}
#' @keywords internal
S7::method(distrib_grad_cdf, SkewTDistrib) <- partial_loc_scale_grad_cdf

#' @rdname distrib_grad_cdf.SkewTDistrib
#' @name distrib_hess_cdf.SkewTDistrib
#' @keywords internal
S7::method(distrib_hess_cdf, SkewTDistrib) <- partial_loc_scale_hess_cdf


# --- positive families with an elementary distribution function -------------

#' @title Lognormal Log-CDF Hessian
#' @name distrib_hess_cdf.Lognormal1Distrib
#' @description
#' Closed form. On the log scale the family is location-scale, so with
#' \eqn{z = (\log q - \mu)/\sigma} and \eqn{\varphi} the standard normal
#' density, \eqn{\partial^2 F/\partial\mu^2 = -z\varphi/\sigma^2},
#' \eqn{\partial^2 F/\partial\mu\partial\sigma^2 = \varphi(1-z^2)/(2\sigma^3)}
#' and
#' \eqn{\partial^2 F/\partial(\sigma^2)^2 = \varphi z(3-z^2)/(4\sigma^4)}.
#' @param distrib A \code{Lognormal1Distrib} object.
#' @param q A numeric vector of quantiles.
#' @param theta A list containing \code{mu} and \code{sigma2}.
#' @param lower.tail Logical; if \code{TRUE} (default), the lower tail.
#' @param log Logical; if \code{TRUE} (default), derivatives of the log probability.
#' @return A named list keyed as \code{\link{hess_names}}.
#' @seealso \code{\link{lognormal1_distrib}}
#' @keywords internal
S7::method(distrib_hess_cdf, Lognormal1Distrib) <-
  function(distrib, q, theta, lower.tail = TRUE, log = TRUE) {
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
