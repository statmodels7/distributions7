#' @include distrib.R generics.R utility_functions.R y_derivatives.R cross_derivatives.R cross2_derivatives.R gaussian1_distrib.R
NULL

# One or two derivatives with respect to the response, TWO with respect to the
# parameters. These complete the mixed grid:
#
#            d^0/dtheta^0        d^1/dtheta^1        d^2/dtheta^2
#   d/dy     distrib_grad_y      distrib_cross_y     distrib_grad_y_hess
#   d2/dy2   distrib_hess_y      distrib_cross2_y    distrib_hess_y_hess
#
# The right-hand column is what the SECOND derivative of a marginal criterion
# needs of a penalty: a penalty is a negative log-density at the coefficients,
# so d2rho/dbeta2 carries the density's response curvature and d3rho/dbeta
# dtheta2 and d4rho/dbeta2 dtheta2 carry these.
#
# Continuous distributions only, like the other response derivatives.

#' @title Hyperparameter Hessian of the Response Gradient
#'
#' @description
#' Computes \eqn{\partial^3 \ell / \partial y\, \partial\theta_i
#' \partial\theta_j}, with one component per unordered pair of parameters, and
#' each component is a vector along `y`. It gives the second derivative, in the
#' parameters, of the response gradient, and it is the third-order entry of the
#' mixed grid below.
#'
#' @details
#' These are the second-order column of the mixed grid whose first-order column
#' is [distrib_cross_y()] and [distrib_cross2_y()]:
#'
#' | | \eqn{\partial^0/\partial\theta^0} | \eqn{\partial^1/\partial\theta^1} | \eqn{\partial^2/\partial\theta^2} |
#' |---|---|---|---|
#' | \eqn{\partial/\partial y} | [distrib_grad_y()] | [distrib_cross_y()] | [distrib_grad_y_hess()] |
#' | \eqn{\partial^2/\partial y^2} | [distrib_hess_y()] | [distrib_cross2_y()] | [distrib_hess_y_hess()] |
#'
#' The second derivative of a marginal criterion requires the right-hand column
#' from a penalty. A penalty is a negative log-density read
#' at the
#' coefficients, so \eqn{\partial^2\rho/\partial\beta^2} carries the density's
#' response curvature and the two mixed derivatives above carry
#' \eqn{\partial^3\rho/\partial\beta\,\partial\theta^2} and
#' \eqn{\partial^4\rho/\partial\beta^2\partial\theta^2}.
#'
#' The components are keyed by [hess_names()], the same enumeration that
#' [distrib_hessian()] uses, so a consumer that looks up a pair finds it under
#' the name that [distrib_hessian()] uses.
#'
#' # The link scale
#'
#' Each component is multiplied by \eqn{h_i'(\eta_i) h_j'(\eta_j)} and, on a
#' diagonal pair, gains \eqn{h_i''(\eta_i)} times the corresponding
#' first-order component from [distrib_cross_y()]. The response derivative is
#' untouched by a reparametrization of \eqn{\theta}, so this is the ordinary
#' diagonal chain rule at second order, exactly as for [distrib_hessian()].
#'
#' # Where the numbers come from
#'
#' A distribution with a closed form provides it. A continuous distribution
#' without one takes [distrib_grad_y_hess.continuous_distrib()], which applies
#' one tensor stencil, through [mixed_tensor_derivatives()], to the highest
#' quantity that the family implements itself, so a single variable is never
#' differenced twice.
#'
#' Methods exist for the continuous families and for the multivariate Gaussian
#' and Student t families. A discrete family has no method, and the call signals
#' an error that names its class.
#'
#' @param distrib A distribution object inheriting from `continuous_distrib`,
#'   or a multivariate Gaussian or Student t distribution.
#' @param y A numeric vector of observations.
#' @param theta A named list, or named numeric vector, of distribution
#'   parameters. The generic aligns it before dispatch.
#' @inheritParams distrib_gradient
#' @param ... Passed to the method.
#'
#' @return A named list with one numeric vector per unordered pair of
#'   parameters. Each vector is as long as `y`, and the list is keyed as
#'   [`hess_names(distrib@params)`][hess_names].
#'
#' @section Notation:
#' \eqn{\ell} is the log-density of one observation, \eqn{y} the response,
#' \eqn{\theta_i} a distribution parameter, \eqn{\eta_i} its value on the
#' unconstrained scale and \eqn{h_i = g_i^{-1}} the inverse link carrying one to
#' the other.
#'
#' @seealso [distrib_hess_y_hess()] for the fourth-order twin,
#'   [distrib_cross_y()] for the first-order column, and
#'   [mixed_tensor_derivatives()] for the fallback.
#'
#' @examples
#' d <- gaussian1_distrib()
#' y <- c(-1, 0, 2)
#' theta <- list(mu = 0.4, sigma = 1.3)
#' distrib_grad_y_hess(d, y, theta)
#'
#' # In closed form, the response gradient is -r / sigma^2, so the location pair
#' # vanishes and the other two are elementary.
#' r <- y - 0.4
#' c(mu_sigma = -2 / 1.3^3, sigma_sigma = -6 * r[1] / 1.3^4)
#'
#' # The result is compared with a numerical Hessian of the response gradient
#' # in the parameters.
#' g <- function(v) distrib_grad_y(d, y[1], list(mu = v[1], sigma = v[2]))
#' numDeriv::hessian(g, c(0.4, 1.3))
#'
#' # On the link scale, with sigma on a log link, both h' and h'' equal sigma.
#' cy <- distrib_cross_y(d, y, theta)
#' distrib_grad_y_hess(d, y, theta, scale = "link")$sigma_sigma
#' distrib_grad_y_hess(d, y, theta)$sigma_sigma * 1.3^2 + cy$sigma * 1.3
#'
#' # A discrete family has no method at all.
#' try(distrib_grad_y_hess(poisson_distrib(), 1:3, list(mu = 2)))
#'
#' @export
distrib_grad_y_hess <- S7::new_generic("distrib_grad_y_hess", "distrib", function(distrib, y, theta, scale = c("parameter", "link"), ...) {
  args <- check_derivative_args(distrib, y, theta)
  y <- args$y
  theta <- args$theta
  scale <- match.arg(scale)
  res <- S7::S7_dispatch()
  if (scale == "link") {
    theta2_link_scale(distrib, y, theta, res, distrib_cross_y(distrib, y, theta))
  } else {
    res
  }
})

#' @title Hyperparameter Hessian of the Response Curvature
#'
#' @description
#' Computes \eqn{\partial^4 \ell / \partial y^2\, \partial\theta_i
#' \partial\theta_j}, with one component per unordered pair of parameters, and
#' each component is a vector along `y`. It gives the second derivative, in the
#' parameters, of the response curvature, and it is the fourth-order entry of
#' the mixed grid below.
#'
#' @details
#' These are the second-order column of the mixed grid whose first-order column
#' is [distrib_cross_y()] and [distrib_cross2_y()]:
#'
#' | | \eqn{\partial^0/\partial\theta^0} | \eqn{\partial^1/\partial\theta^1} | \eqn{\partial^2/\partial\theta^2} |
#' |---|---|---|---|
#' | \eqn{\partial/\partial y} | [distrib_grad_y()] | [distrib_cross_y()] | [distrib_grad_y_hess()] |
#' | \eqn{\partial^2/\partial y^2} | [distrib_hess_y()] | [distrib_cross2_y()] | [distrib_hess_y_hess()] |
#'
#' The second derivative of a marginal criterion requires the right-hand column
#' from a penalty. A penalty is a negative log-density read
#' at the
#' coefficients, so \eqn{\partial^2\rho/\partial\beta^2} carries the density's
#' response curvature and the two mixed derivatives above carry
#' \eqn{\partial^3\rho/\partial\beta\,\partial\theta^2} and
#' \eqn{\partial^4\rho/\partial\beta^2\partial\theta^2}.
#'
#' The components are keyed by [hess_names()], the same enumeration that
#' [distrib_hessian()] uses.
#'
#' # The link scale
#'
#' Each component is multiplied by \eqn{h_i'(\eta_i) h_j'(\eta_j)} and, on a
#' diagonal pair, gains \eqn{h_i''(\eta_i)} times the corresponding
#' first-order component from [distrib_cross2_y()]. That is the same diagonal
#' chain rule that [distrib_grad_y_hess()] obeys, except that the first-order
#' quantity is differentiated one order higher in \eqn{y}.
#'
#' # Where the numbers come from
#'
#' A distribution with a closed form provides it. A continuous distribution
#' without one takes [distrib_hess_y_hess.continuous_distrib()], which applies
#' one tensor stencil, through [mixed_tensor_derivatives()], to the highest
#' quantity that the family implements itself. The multivariate Gaussian and
#' Student t families have methods of their own; a discrete family has no
#' method, and the call signals an error that names its class.
#'
#' @param distrib A distribution object inheriting from `continuous_distrib`,
#'   or a multivariate Gaussian or Student t distribution.
#' @param y A numeric vector of observations.
#' @param theta A named list, or named numeric vector, of distribution
#'   parameters. The generic aligns it before dispatch.
#' @inheritParams distrib_gradient
#' @param ... Passed to the method.
#'
#' @return A named list with one numeric vector per unordered pair of
#'   parameters. Each vector is as long as `y`, and the list is keyed as
#'   [`hess_names(distrib@params)`][hess_names].
#'
#' @section Notation:
#' \eqn{\ell} is the log-density of one observation, \eqn{y} the response,
#' \eqn{\theta_i} a distribution parameter, \eqn{\eta_i} its value on the
#' unconstrained scale and \eqn{h_i = g_i^{-1}} the inverse link carrying one to
#' the other.
#'
#' @seealso [distrib_grad_y_hess()] for the third-order twin,
#'   [distrib_cross2_y()] for the first-order column, and
#'   [mixed_tensor_derivatives()] for the fallback.
#'
#' @examples
#' d <- gaussian1_distrib()
#' y <- c(-1, 0, 2)
#' theta <- list(mu = 0.4, sigma = 1.3)
#' distrib_hess_y_hess(d, y, theta)
#'
#' # The gaussian's response curvature is -1 / sigma^2, carrying no location,
#' # so only the scale pair survives and it does not vary with y.
#' -6 / 1.3^4
#'
#' # The result is compared with a numerical Hessian of the response curvature
#' # in the parameters.
#' h <- function(v) distrib_hess_y(d, y[1], list(mu = v[1], sigma = v[2]))
#' numDeriv::hessian(h, c(0.4, 1.3))
#'
#' # A discrete family has no method at all.
#' try(distrib_hess_y_hess(poisson_distrib(), 1:3, list(mu = 2)))
#'
#' @export
distrib_hess_y_hess <- S7::new_generic("distrib_hess_y_hess", "distrib", function(distrib, y, theta, scale = c("parameter", "link"), ...) {
  args <- check_derivative_args(distrib, y, theta)
  y <- args$y
  theta <- args$theta
  scale <- match.arg(scale)
  res <- S7::S7_dispatch()
  if (scale == "link") {
    theta2_link_scale(distrib, y, theta, res,
                      distrib_cross2_y(distrib, y, theta))
  } else {
    res
  }
})


#' @title The Link Scale of a Second-Order Parameter Derivative
#'
#' @description
#' Carries a quantity keyed by parameter pair from the parameter scale onto the
#' unconstrained one. Both [distrib_grad_y_hess()] and [distrib_hess_y_hess()]
#' call it after dispatch, so a method returns the parameter scale and contains
#' no code for links.
#'
#' @details
#' The chain rule is DIAGONAL because each parameter carries its own link, so a
#' pair \eqn{(i, j)} becomes
#' \deqn{\frac{\partial^2}{\partial\eta_i \partial\eta_j}
#'   = h_i'(\eta_i)\, h_j'(\eta_j) \frac{\partial^2}{\partial\theta_i
#'     \partial\theta_j}
#'   + \delta_{ij}\, h_i''(\eta_i) \frac{\partial}{\partial\theta_i}.}
#' The response derivatives do not enter, because a reparametrization of
#' \eqn{\theta} leaves them alone. For that reason `first` and `second` are the
#' same quantity at two orders in \eqn{\theta}, and the function takes no third
#' argument.
#'
#' @param distrib A distribution object. Its `link_params` supply \eqn{h'} and
#'   \eqn{h''} through [inverse_link_derivs()].
#' @param y The response. It is unused and is present so that the signature
#'   matches that of the callers.
#' @param theta A named list of parameters, on the parameter scale.
#' @param second The parameter-scale second-order components, keyed by
#'   parameter pair.
#' @param first The parameter-scale first-order components, keyed by parameter:
#'   [distrib_cross_y()] for the third-order quantity and [distrib_cross2_y()]
#'   for the fourth.
#'
#' @return A named list keyed exactly as `second`.
#'
#' @section Notation:
#' \eqn{\ell} is the log-density of one observation, \eqn{y} the response,
#' \eqn{\theta_i} a distribution parameter, \eqn{\eta_i} its value on the
#' unconstrained scale and \eqn{h_i = g_i^{-1}} the inverse link carrying one to
#' the other.
#'
#' @seealso [distrib_grad_y_hess()] and [distrib_hess_y_hess()], its two
#'   callers, and [inverse_link_derivs()] for the link derivatives.
#'
#' @keywords internal
#'
#' @examples
#' d <- gaussian1_distrib()
#' y <- c(-1, 0, 2)
#' theta <- list(mu = 0.4, sigma = 1.3)
#'
#' # The generic applies the same computation when scale = "link".
#' second <- distrib_grad_y_hess(d, y, theta)
#' first <- distrib_cross_y(d, y, theta)
#' linked <- distributions7:::theta2_link_scale(d, y, theta, second, first)
#' identical(linked, distrib_grad_y_hess(d, y, theta, scale = "link"))
#'
#' # sigma rides a log link, so h' = h'' = sigma and the diagonal pair gains
#' # the first-order term while the off-diagonal one does not.
#' c(diagonal = second$sigma_sigma[1] * 1.3^2 + first$sigma[1] * 1.3,
#'   off = second$mu_sigma[1] * 1 * 1.3)
#' c(linked$sigma_sigma[1], linked$mu_sigma[1])
theta2_link_scale <- function(distrib, y, theta, second, first) {
  params <- distrib@params
  # indexed by PARAMETER first and by ORDER second: d[[i]][[1]] is h'_i and
  # d[[i]][[2]] is h''_i
  d <- inverse_link_derivs(distrib, theta, 2L)
  prs <- hess_pairs(params)
  stats::setNames(lapply(names(second), function(nm) {
    pr <- prs[[nm]]
    i <- match(pr[1L], params)
    j <- match(pr[2L], params)
    out <- second[[nm]] * d[[i]][[1L]] * d[[j]][[1L]]
    if (i == j) out <- out + first[[pr[1L]]] * d[[i]][[2L]]
    out
  }), names(second))
}


#' @title Numerical Hyperparameter Hessians of the Response Derivatives
#'
#' @description
#' Computes a second-order mixed component by one central difference, in each
#' parameter, of a first-order component. It supplies the pairs that have no
#' closed form in the methods for location-scale families with a shape parameter
#' and for families with a scale and no location. The default methods for a
#' continuous family use [mixed_tensor_derivatives()] instead.
#'
#' @details
#' The reference is [distrib_cross_y()] or [distrib_cross2_y()]. If a
#' distribution has a closed form for those, the function applies exactly one
#' difference to that closed form.
#'
#' A mixed pair is differenced in both orders and the two results are averaged.
#' They agree in exact arithmetic and differ in floating point because the two
#' steps have different sizes. The average is symmetric in the pair.
#'
#' The step is chosen by [fd_stable_step()] from the two candidates of
#' [fd_steps()], so that both evaluation points of the central difference of a
#' parameter near a bound stay inside its domain.
#'
#' @param distrib A distribution object. Its `params` and `params_bounds` set
#'   the enumeration and the steps.
#' @param y A numeric vector of observations. The function does not use it,
#'   because `inner` has already closed over it.
#' @param theta A named list of parameters, the point at which to differentiate.
#' @param inner A function that takes `theta` as its only argument and returns
#'   the first-order components as a list with one entry per parameter,
#'   typically `function(th) distrib_cross_y(distrib, y, th)`.
#' @param rel_step The relative step, defaulting to `.Machine$double.eps^(1/3)`,
#'   which is the optimal exponent for a central first difference.
#'
#' @return A named list with one numeric vector per unordered pair of
#'   parameters, keyed as [`hess_names(distrib@params)`][hess_names].
#'
#' @section Notation:
#' \eqn{\ell} is the log-density of one observation, \eqn{y} the response,
#' \eqn{\theta_i} a distribution parameter, \eqn{\eta_i} its value on the
#' unconstrained scale and \eqn{h_i = g_i^{-1}} the inverse link carrying one to
#' the other.
#'
#' @seealso [distrib_grad_y_hess()] and [distrib_hess_y_hess()], the two
#'   generics that it serves, and [fd_steps()] for the step rule.
#'
#' @examples
#' d <- gaussian1_distrib()
#' y <- c(-1, 0, 2)
#' theta <- list(mu = 0.4, sigma = 1.3)
#'
#' num <- numerical_theta2_y(d, y, theta,
#'                           function(th) distrib_cross_y(d, y, th))
#' num$sigma_sigma
#'
#' # The result is compared with the gaussian's own closed form, for which this
#' # family does not need the fallback.
#' max(abs(unlist(num) - unlist(distrib_grad_y_hess(d, y, theta))))
#'
#' # A family that takes the default continuous method is checked against
#' # numDeriv.
#' dg <- gamma2_distrib()
#' yg <- c(0.5, 1, 2)
#' thg <- list(mu = 2, sigma2 = 1)
#' g <- distrib_grad_y_hess(dg, yg, thg)
#' f <- function(v) distrib_grad_y(dg, yg[1], list(mu = v[1], sigma2 = v[2]))
#' numDeriv::hessian(f, c(2, 1))
#' c(g$mu_mu[1], g$mu_sigma2[1], g$sigma2_sigma2[1])
#'
#' @export
numerical_theta2_y <- function(distrib, y, theta, inner,
                               rel_step = .Machine$double.eps^(1 / 3)) {
  params <- distrib@params
  diffs <- stats::setNames(vector("list", length(params)), params)
  for (j in seq_along(params)) {
    p <- params[j]
    quotient <- function(h) {
      th_up <- th_dn <- theta
      th_up[[j]] <- theta[[j]] + h
      th_dn[[j]] <- theta[[j]] - h
      up <- inner(th_up)
      dn <- inner(th_dn)
      stats::setNames(lapply(params, function(q) (up[[q]] - dn[[q]]) / (2 * h)), params)
    }
    diffs[[p]] <- fd_stable_step(quotient, theta[[j]],
                                 distrib@params_bounds[[p]], rel_step)$value
  }
  prs <- hess_pairs(params)
  stats::setNames(lapply(names(prs), function(nm) {
    pr <- prs[[nm]]
    (diffs[[pr[1L]]][[pr[2L]]] + diffs[[pr[2L]]][[pr[1L]]]) / 2
  }), names(prs))
}


#' @title Default Hyperparameter Hessian of the Response Gradient
#' @name distrib_grad_y_hess.continuous_distrib
#'
#' @description
#' Falls back to one tensor stencil on the highest quantity that the family
#' implements itself, through [mixed_tensor_derivatives()]: a central
#' difference of an analytic [distrib_cross_y()] in each parameter, a
#' difference in the response of an analytic Hessian, or one mixed stencil of
#' the log-density, whichever leaves the smallest order. Registering the
#' fallback on `continuous_distrib` gives the third-order mixed derivative to
#' every continuous family, whether or not it writes one out.
#'
#' @details
#' The stencil is never applied to a quantity that is itself a numerical
#' fallback; only quantities that the family's class registers are used as the
#' base.
#'
#' @param distrib A `continuous_distrib` object with no closed form of its own.
#' @param y A numeric vector of observations.
#' @param theta A named list of parameters.
#' @param scale One of `"parameter"` or `"link"`. The generic applies the link
#'   scale to the value that the method returns.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A named list with one numeric vector per unordered pair of
#'   parameters, keyed as [`hess_names(distrib@params)`][hess_names].
#'
#' @seealso [distrib_grad_y_hess()] for the generic,
#'   [distrib_hess_y_hess.continuous_distrib()] for the fourth-order twin, and
#'   [mixed_tensor_derivatives()], which does the work.
#'
#' @keywords internal
#'
#' @examples
#' # The gamma has no closed form of its own, so this method is used for it.
#' d <- gamma2_distrib()
#' y <- c(0.5, 1, 2)
#' theta <- list(mu = 2, sigma2 = 1)
#' attr(S7::method(distrib_grad_y_hess, S7::S7_class(d)), "signature")[[1]]
#'
#' g <- distrib_grad_y_hess(d, y, theta)
#' c(g$mu_mu[1], g$mu_sigma2[1], g$sigma2_sigma2[1])
#'
#' # The result is compared with a numerical Hessian of the response gradient.
#' f <- function(v) distrib_grad_y(d, y[1], list(mu = v[1], sigma2 = v[2]))
#' numDeriv::hessian(f, c(2, 1))
S7::method(distrib_grad_y_hess, continuous_distrib) <- function(distrib, y, theta,
                                                                scale = c("parameter", "link"),
                                                                ...) {
  mixed_tensor_derivatives(distrib, y, theta, 1L, 2L)
}

#' @title Default Hyperparameter Hessian of the Response Curvature
#' @name distrib_hess_y_hess.continuous_distrib
#'
#' @description
#' Falls back to one tensor stencil on the highest quantity that the family
#' implements itself, through [mixed_tensor_derivatives()]. It is
#' [distrib_grad_y_hess.continuous_distrib()] read one order higher in the
#' response, and it makes the fourth-order mixed derivative available for
#' every continuous family. For [transformation()], which has an analytic
#' Hessian in the parameters and no response derivative, it is a second
#' difference in the response of that Hessian.
#'
#' @param distrib A `continuous_distrib` object with no closed form of its own.
#' @param y A numeric vector of observations.
#' @param theta A named list of parameters.
#' @param scale One of `"parameter"` or `"link"`. The generic applies the link
#'   scale to the value that the method returns.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A named list with one numeric vector per unordered pair of
#'   parameters, keyed as [`hess_names(distrib@params)`][hess_names].
#'
#' @seealso [distrib_hess_y_hess()] for the generic,
#'   [distrib_grad_y_hess.continuous_distrib()] for the third-order twin, and
#'   [mixed_tensor_derivatives()], which does the work.
#'
#' @keywords internal
#'
#' @examples
#' d <- gamma2_distrib()
#' y <- c(0.5, 1, 2)
#' theta <- list(mu = 2, sigma2 = 1)
#' h <- distrib_hess_y_hess(d, y, theta)
#' c(h$mu_mu[1], h$mu_sigma2[1], h$sigma2_sigma2[1])
#'
#' # The result is compared with a numerical Hessian of the response curvature.
#' f <- function(v) distrib_hess_y(d, y[1], list(mu = v[1], sigma2 = v[2]))
#' numDeriv::hessian(f, c(2, 1))
S7::method(distrib_hess_y_hess, continuous_distrib) <- function(distrib, y, theta,
                                                                scale = c("parameter", "link"),
                                                                ...) {
  mixed_tensor_derivatives(distrib, y, theta, 2L, 2L)
}


#' @title Gaussian Hyperparameter Hessian of the Response Gradient
#' @name distrib_grad_y_hess.Gaussian1Distrib
#'
#' @description
#' Computes the derivative in closed form. With \eqn{r = y - \mu} the response
#' gradient is
#' \eqn{-r/\sigma^2}, so differentiating it twice in the parameters gives
#' \deqn{\frac{\partial^3\ell}{\partial y\,\partial\mu^2} = 0, \qquad
#'   \frac{\partial^3\ell}{\partial y\,\partial\mu\,\partial\sigma}
#'     = -\frac{2}{\sigma^3}, \qquad
#'   \frac{\partial^3\ell}{\partial y\,\partial\sigma^2}
#'     = -\frac{6r}{\sigma^4}.}
#' The location pair vanishes because the response gradient is LINEAR in
#' \eqn{\mu}, and only the scale pair varies with the data.
#'
#' @param distrib A `Gaussian1Distrib` object.
#' @param y A numeric vector of observations.
#' @param theta A named list with `mu` and `sigma`.
#' @param scale One of `"parameter"` or `"link"`. The generic applies the link
#'   scale to the value that the method returns.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A named list of three numeric vectors of length `length(y)`, keyed
#'   `mu_mu`, `sigma_sigma` and `mu_sigma`.
#'
#' @section Notation:
#' \eqn{\ell} is the log-density of one observation, \eqn{y} the response,
#' \eqn{\theta_i} a distribution parameter, \eqn{\eta_i} its value on the
#' unconstrained scale and \eqn{h_i = g_i^{-1}} the inverse link carrying one to
#' the other.
#'
#' @seealso [distrib_grad_y_hess()] for the generic,
#'   [distrib_hess_y_hess.Gaussian1Distrib()] for the fourth-order twin, and
#'   [distrib_cross_y()] for the first-order column.
#'
#' @keywords internal
#'
#' @examples
#' d <- gaussian1_distrib()
#' y <- c(-1, 0, 2)
#' theta <- list(mu = 0.4, sigma = 1.3)
#' g <- distrib_grad_y_hess(d, y, theta)
#' g
#'
#' # The three formulas are written out.
#' r <- y - 0.4
#' c(mu_mu = 0, mu_sigma = -2 / 1.3^3, sigma_sigma = -6 * r[1] / 1.3^4)
#'
#' # The result is compared with a numerical Hessian of the response gradient.
#' f <- function(v) distrib_grad_y(d, y[1], list(mu = v[1], sigma = v[2]))
#' numDeriv::hessian(f, c(0.4, 1.3))
S7::method(distrib_grad_y_hess, Gaussian1Distrib) <- function(distrib, y, theta,
                                                              scale = c("parameter", "link"),
                                                              ...) {
  mu <- theta[[1]]
  sigma <- theta[[2]]
  r <- y - mu
  n <- length(y)
  list(mu_mu = rep(0, length.out = n),
       sigma_sigma = rep_len(-6 * r / sigma^4, n),
       mu_sigma = rep(-2 / sigma^3, length.out = n))
}

#' @title Gaussian Hyperparameter Hessian of the Response Curvature
#' @name distrib_hess_y_hess.Gaussian1Distrib
#'
#' @description
#' Computes the closed form, in which nearly every component is zero. The
#' gaussian's response curvature is
#' \eqn{-1/\sigma^2}, which carries no location at all, so
#' \deqn{\frac{\partial^4\ell}{\partial y^2\partial\sigma^2}
#'     = -\frac{6}{\sigma^4}}
#' is the only component that does not vanish, and it does not vary with the
#' data either.
#'
#' @param distrib A `Gaussian1Distrib` object.
#' @param y A numeric vector of observations. Only its length is used.
#' @param theta A named list with `mu` and `sigma`.
#' @param scale One of `"parameter"` or `"link"`. The generic applies the link
#'   scale to the value that the method returns.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A named list of three numeric vectors of length `length(y)`, keyed
#'   `mu_mu`, `sigma_sigma` and `mu_sigma`, of which the first and last are
#'   exactly zero.
#'
#' @section Notation:
#' \eqn{\ell} is the log-density of one observation, \eqn{y} the response,
#' \eqn{\theta_i} a distribution parameter, \eqn{\eta_i} its value on the
#' unconstrained scale and \eqn{h_i = g_i^{-1}} the inverse link carrying one to
#' the other.
#'
#' @seealso [distrib_hess_y_hess()] for the generic,
#'   [distrib_grad_y_hess.Gaussian1Distrib()] for the third-order twin, and
#'   [distrib_cross2_y()] for the first-order column.
#'
#' @keywords internal
#'
#' @examples
#' d <- gaussian1_distrib()
#' y <- c(-1, 0, 2)
#' theta <- list(mu = 0.4, sigma = 1.3)
#' h <- distrib_hess_y_hess(d, y, theta)
#' h
#'
#' # One component survives, and it is constant along y.
#' c(reported = h$sigma_sigma[1], formula = -6 / 1.3^4)
#'
#' # The result is compared with a numerical Hessian of the response curvature.
#' f <- function(v) distrib_hess_y(d, y[1], list(mu = v[1], sigma = v[2]))
#' numDeriv::hessian(f, c(0.4, 1.3))
S7::method(distrib_hess_y_hess, Gaussian1Distrib) <- function(distrib, y, theta,
                                                              scale = c("parameter", "link"),
                                                              ...) {
  sigma <- theta[[2]]
  n <- length(y)
  list(mu_mu = rep(0, length.out = n),
       sigma_sigma = rep(-6 / sigma^4, length.out = n),
       mu_sigma = rep(0, length.out = n))
}
