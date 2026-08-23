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

#' Hyperparameter Hessians of the Response Derivatives
#'
#' @description
#' `distrib_grad_y_hess()` computes
#' \eqn{\partial^3 \ell / \partial y\, \partial\theta_i \partial\theta_j} and
#' `distrib_hess_y_hess()` computes
#' \eqn{\partial^4 \ell / \partial y^2\, \partial\theta_i \partial\theta_j},
#' one component per unordered pair of parameters, each a vector along
#' `y`.
#'
#' @details
#' These are the second-order column of the mixed grid whose first-order one
#' is [distrib_cross_y()] and [distrib_cross2_y()]: how the
#' response gradient and the response curvature CURVE in the parameters. A
#' marginal likelihood needs them to differentiate its Laplace approximation
#' twice, the penalty being a negative log-density evaluated at the
#' coefficients.
#'
#' The components are keyed by [hess_names()], the same enumeration
#' [distrib_hessian()] uses, so a consumer that looks a pair up finds
#' it under the name it already knows.
#'
#' On the link scale each component is multiplied by \eqn{h_i'(\eta_i)
#' h_j'(\eta_j)} and, on a diagonal pair, gains the second-order term
#' \eqn{h_i''(\eta_i)} times the corresponding first-order component: the
#' response derivatives are untouched by a reparametrization of \eqn{\theta},
#' so this is the ordinary diagonal chain rule at second order, exactly as for
#' [distrib_hessian()].
#'
#' Distributions with closed forms provide them; the others take one central
#' difference of the analytic first-order quantity in each parameter (see
#' [numerical_theta2_y()]). The two differences act on different
#' parameters off the diagonal, so they compose into a single mixed stencil
#' rather than the nested differencing of one variable the package forbids.
#'
#' @param distrib A distribution object inheriting from the `distrib`
#'   class.
#' @param y A numeric vector of observations.
#' @param theta A named list (or named numeric vector) of distribution
#'   parameters.
#' @inheritParams distrib_gradient
#' @param ... Additional arguments passed to the specific method.
#'
#' @return A named list with one numeric vector per parameter pair, keyed by
#'   `hess_names(distrib@params)`.
#'
#' @examples
#' d <- gaussian1_distrib()
#' distrib_grad_y_hess(d, c(-1, 0, 2), list(mu = 0, sigma = 1))
#' distrib_hess_y_hess(d, c(-1, 0, 2), list(mu = 0, sigma = 1))
#'
#' @seealso [distrib_cross_y()], [distrib_cross2_y()]
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

#' @rdname distrib_grad_y_hess
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


#' The Link Scale of a Second-Order Parameter Derivative
#'
#' @description
#' Carries a component keyed by parameter pair onto the unconstrained scale.
#'
#' @details
#' The chain rule is diagonal, each parameter having its own link, so a pair
#' \eqn{(i, j)} is multiplied by \eqn{h_i' h_j'} and a diagonal pair gains
#' \eqn{h_i''} times the first-order component. The response derivatives do not
#' enter it: a reparametrization of \eqn{\theta} leaves them alone.
#'
#' @param distrib A distribution object.
#' @param y The response.
#' @param theta The parameters.
#' @param second The parameter-scale second-order components.
#' @param first The parameter-scale first-order components.
#'
#' @return A named list keyed as `second`.
#'
#' @keywords internal
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


#' Numerical Hyperparameter Hessians of the Response Derivatives
#'
#' @description
#' Computes the second-order mixed components by one central difference of an
#' analytic first-order component in each parameter.
#'
#' @details
#' The reference is [distrib_cross_y()] or
#' [distrib_cross2_y()], so a distribution with a closed form for
#' those pays for exactly one difference. A mixed pair is differenced both ways
#' and averaged: the two agree in exact arithmetic and not quite in floating
#' point, the steps differing, and a second derivative of a scalar has to come
#' out symmetric.
#'
#' @param distrib A distribution object.
#' @param y A numeric vector of observations.
#' @param theta A named list of parameters.
#' @param inner A function of `theta` returning the first-order
#'   components, one per parameter.
#' @param h_rel The relative step.
#'
#' @return A named list keyed by `hess_names(distrib@params)`.
#'
#' @examples
#' d <- gaussian1_distrib()
#' numerical_theta2_y(d, c(-1, 0, 1), list(mu = 0, sigma = 1),
#'                    function(th) distrib_cross_y(d, c(-1, 0, 1), th))
#'
#' @export
numerical_theta2_y <- function(distrib, y, theta, inner,
                               h_rel = .Machine$double.eps^(1 / 3)) {
  params <- distrib@params
  diffs <- stats::setNames(vector("list", length(params)), params)
  for (j in seq_along(params)) {
    p <- params[j]
    h <- fd_steps(theta[[j]], distrib@params_bounds[[p]], h_rel)
    th_up <- theta
    th_dn <- theta
    th_up[[j]] <- theta[[j]] + h
    th_dn[[j]] <- theta[[j]] - h
    up <- inner(th_up)
    dn <- inner(th_dn)
    diffs[[p]] <- stats::setNames(
      lapply(params, function(q) (up[[q]] - dn[[q]]) / (2 * h)), params)
  }
  prs <- hess_pairs(params)
  stats::setNames(lapply(names(prs), function(nm) {
    pr <- prs[[nm]]
    (diffs[[pr[1L]]][[pr[2L]]] + diffs[[pr[2L]]][[pr[1L]]]) / 2
  }), names(prs))
}


#' @title Default Second-Order Mixed Derivatives for Continuous Distributions
#' @name distrib_grad_y_hess.continuous_distrib
#' @description Fallback: one central difference of the analytic
#'   [distrib_cross_y()] or [distrib_cross2_y()] in each
#'   parameter (see [numerical_theta2_y()]).
#' @param distrib A `continuous_distrib` object.
#' @param y A numeric vector of observations.
#' @param theta A named list of parameters.
#' @param scale Handled by the generic before dispatch.
#' @param ... Unused.
#' @return A named list keyed by parameter pair.
#' @keywords internal
S7::method(distrib_grad_y_hess, continuous_distrib) <- function(distrib, y, theta,
                                                                scale = c("parameter", "link"),
                                                                ...) {
  numerical_theta2_y(distrib, y, theta,
                     function(th) distrib_cross_y(distrib, y, th))
}

#' @rdname distrib_grad_y_hess.continuous_distrib
#' @name distrib_hess_y_hess.continuous_distrib
#' @keywords internal
S7::method(distrib_hess_y_hess, continuous_distrib) <- function(distrib, y, theta,
                                                                scale = c("parameter", "link"),
                                                                ...) {
  numerical_theta2_y(distrib, y, theta,
                     function(th) distrib_cross2_y(distrib, y, th))
}


#' @title Gaussian Second-Order Mixed Derivatives
#' @name distrib_grad_y_hess.Gaussian1Distrib
#' @description
#' Closed forms. With \eqn{r = y - \mu}, the response gradient is
#' \eqn{-r/\sigma^2} and the response curvature \eqn{-1/\sigma^2}, so
#' \deqn{\partial^3\ell/\partial y\,\partial\mu^2 = 0, \quad
#'   \partial^3\ell/\partial y\,\partial\mu\partial\sigma = -2/\sigma^3, \quad
#'   \partial^3\ell/\partial y\,\partial\sigma^2 = -6r/\sigma^4,}
#' and the only component of the fourth derivative that does not vanish is
#' \eqn{\partial^4\ell/\partial y^2\partial\sigma^2 = -6/\sigma^4}, the
#' curvature carrying no location at all.
#' @param distrib A `Gaussian1Distrib` object.
#' @param y A numeric vector of observations.
#' @param theta A list containing `mu` and `sigma`.
#' @param scale Handled by the generic before dispatch.
#' @param ... Unused.
#' @return A named list keyed by parameter pair.
#' @keywords internal
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

#' @rdname distrib_grad_y_hess.Gaussian1Distrib
#' @name distrib_hess_y_hess.Gaussian1Distrib
#' @keywords internal
S7::method(distrib_hess_y_hess, Gaussian1Distrib) <- function(distrib, y, theta,
                                                              scale = c("parameter", "link"),
                                                              ...) {
  sigma <- theta[[2]]
  n <- length(y)
  list(mu_mu = rep(0, length.out = n),
       sigma_sigma = rep(-6 / sigma^4, length.out = n),
       mu_sigma = rep(0, length.out = n))
}
