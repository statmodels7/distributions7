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
#' \partial\theta_j}, one component per unordered pair of parameters, each a
#' vector along `y`. It says how the response GRADIENT curves in the
#' parameters, and it is the third-order entry of the mixed grid below.
#'
#' @details
#' These are the second-order column of the mixed grid whose first-order column
#' is [distrib_cross_y()] and [distrib_cross2_y()]:
#'
#' |          | \eqn{\partial^0/\partial\theta^0} | \eqn{\partial^1/\partial\theta^1} | \eqn{\partial^2/\partial\theta^2} |
#' |---|---|---|---|
#' | \eqn{\partial/\partial y}   | [distrib_grad_y()] | [distrib_cross_y()]  | [distrib_grad_y_hess()] |
#' | \eqn{\partial^2/\partial y^2} | [distrib_hess_y()] | [distrib_cross2_y()] | [distrib_hess_y_hess()] |
#'
#' The right-hand column is the one a marginal criterion's SECOND derivative
#' asks of a penalty. A penalty is a negative log-density read at the
#' coefficients, so \eqn{\partial^2\rho/\partial\beta^2} carries the density's
#' response curvature and the two mixed derivatives above carry
#' \eqn{\partial^3\rho/\partial\beta\,\partial\theta^2} and
#' \eqn{\partial^4\rho/\partial\beta^2\partial\theta^2}.
#'
#' The components are keyed by [hess_names()], the same enumeration
#' [distrib_hessian()] uses, so a consumer looking a pair up finds it under the
#' name it already knows.
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
#' A distribution with a closed form provides it. The rest take one central
#' difference of the analytic [distrib_cross_y()] in each parameter, through
#' [numerical_theta2_y()]. Off the diagonal the two differences act on
#' DIFFERENT parameters, so they compose into one mixed stencil rather than the
#' nested differencing of a single variable the package forbids.
#'
#' Continuous distributions only, as with every other response derivative: a
#' discrete family has no method and the call raises with its class named.
#'
#' @param distrib A distribution object inheriting from `continuous_distrib`.
#' @param y A numeric vector of observations.
#' @param theta A named list, or named numeric vector, of distribution
#'   parameters. Aligned by the generic before dispatch.
#' @inheritParams distrib_gradient
#' @param ... Passed to the method.
#'
#' @return A named list with one numeric vector per unordered pair of
#'   parameters, each as long as `y`, keyed as
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
#'   [numerical_theta2_y()] for the fallback.
#'
#' @examples
#' d <- gaussian1_distrib()
#' y <- c(-1, 0, 2)
#' theta <- list(mu = 0.4, sigma = 1.3)
#' distrib_grad_y_hess(d, y, theta)
#'
#' # Closed form: the response gradient is -r / sigma^2, so the location pair
#' # vanishes and the other two are elementary.
#' r <- y - 0.4
#' c(mu_sigma = -2 / 1.3^3, sigma_sigma = -6 * r[1] / 1.3^4)
#'
#' # Against a numerical Hessian of the response gradient in the parameters.
#' g <- function(v) distrib_grad_y(d, y[1], list(mu = v[1], sigma = v[2]))
#' numDeriv::hessian(g, c(0.4, 1.3))
#'
#' # On the link scale, with sigma on a log link: h' = h'' = sigma.
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
#' \partial\theta_j}, one component per unordered pair of parameters, each a
#' vector along `y`. It says how the response CURVATURE curves in the
#' parameters, and it is the fourth-order entry of the mixed grid below.
#'
#' @details
#' These are the second-order column of the mixed grid whose first-order column
#' is [distrib_cross_y()] and [distrib_cross2_y()]:
#'
#' |          | \eqn{\partial^0/\partial\theta^0} | \eqn{\partial^1/\partial\theta^1} | \eqn{\partial^2/\partial\theta^2} |
#' |---|---|---|---|
#' | \eqn{\partial/\partial y}   | [distrib_grad_y()] | [distrib_cross_y()]  | [distrib_grad_y_hess()] |
#' | \eqn{\partial^2/\partial y^2} | [distrib_hess_y()] | [distrib_cross2_y()] | [distrib_hess_y_hess()] |
#'
#' The right-hand column is the one a marginal criterion's SECOND derivative
#' asks of a penalty. A penalty is a negative log-density read at the
#' coefficients, so \eqn{\partial^2\rho/\partial\beta^2} carries the density's
#' response curvature and the two mixed derivatives above carry
#' \eqn{\partial^3\rho/\partial\beta\,\partial\theta^2} and
#' \eqn{\partial^4\rho/\partial\beta^2\partial\theta^2}.
#'
#' The components are keyed by [hess_names()], the same enumeration
#' [distrib_hessian()] uses.
#'
#' # The link scale
#'
#' Each component is multiplied by \eqn{h_i'(\eta_i) h_j'(\eta_j)} and, on a
#' diagonal pair, gains \eqn{h_i''(\eta_i)} times the corresponding
#' first-order component from [distrib_cross2_y()]. That is the same diagonal
#' chain rule [distrib_grad_y_hess()] obeys, with the first-order quantity
#' taken one order higher in \eqn{y}.
#'
#' # Where the numbers come from
#'
#' A distribution with a closed form provides it. The rest take one central
#' difference of the analytic [distrib_cross2_y()] in each parameter, through
#' [numerical_theta2_y()]. Continuous distributions only: a discrete family has
#' no method and the call raises with its class named.
#'
#' @param distrib A distribution object inheriting from `continuous_distrib`.
#' @param y A numeric vector of observations.
#' @param theta A named list, or named numeric vector, of distribution
#'   parameters. Aligned by the generic before dispatch.
#' @inheritParams distrib_gradient
#' @param ... Passed to the method.
#'
#' @return A named list with one numeric vector per unordered pair of
#'   parameters, each as long as `y`, keyed as
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
#'   [numerical_theta2_y()] for the fallback.
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
#' # Against a numerical Hessian of the response curvature in the parameters.
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
#' call it after dispatch, so a method returns the parameter scale and never
#' has to know about links.
#'
#' @details
#' The chain rule is DIAGONAL, each parameter carrying its own link, so a pair
#' \eqn{(i, j)} becomes
#' \deqn{\frac{\partial^2}{\partial\eta_i \partial\eta_j}
#'   = h_i'(\eta_i)\, h_j'(\eta_j) \frac{\partial^2}{\partial\theta_i
#'     \partial\theta_j}
#'   + \delta_{ij}\, h_i''(\eta_i) \frac{\partial}{\partial\theta_i}.}
#' The response derivatives do not enter: a reparametrization of \eqn{\theta}
#' leaves them alone, which is why `first` and `second` are the SAME quantity
#' at two orders in \eqn{\theta} and no third argument is needed.
#'
#' @param distrib A distribution object. Its `link_params` supply \eqn{h'} and
#'   \eqn{h''} through [inverse_link_derivs()].
#' @param y The response. Unused, and present so that the signature reads like
#'   its callers'.
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
#' # What the generic does when scale = "link" is asked for.
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
#' parameter, of an ANALYTIC first-order component. It is the fallback behind
#' [distrib_grad_y_hess()] and [distrib_hess_y_hess()] for a family that
#' provides no closed form.
#'
#' @details
#' The reference is [distrib_cross_y()] or [distrib_cross2_y()], so a
#' distribution with a closed form for those pays for exactly one difference
#' and never differences a difference.
#'
#' A mixed pair is differenced BOTH WAYS and averaged. The two agree in exact
#' arithmetic, and in floating point they differ, the two steps being
#' different sizes, while a second derivative of a scalar has to come out
#' symmetric. On a gamma the two orders differ by about
#' \eqn{3\times 10^{-10}}, which is the size of the answer's own error.
#'
#' The step is [fd_steps()]'s, so a parameter near a bound is differenced
#' inward rather than across it.
#'
#' @param distrib A distribution object. Its `params` and `params_bounds` set
#'   the enumeration and the steps.
#' @param y A numeric vector of observations. Passed to nothing here; `inner`
#'   has already closed over it.
#' @param theta A named list of parameters, the point to differentiate at.
#' @param inner A function of `theta` alone returning the first-order
#'   components as a list with one entry per parameter, typically
#'   `function(th) distrib_cross_y(distrib, y, th)`.
#' @param h_rel The relative step, defaulting to `.Machine$double.eps^(1/3)`,
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
#'   generics it serves, and [fd_steps()] for the step rule.
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
#' # Against the gaussian's own closed form, which this family does not need
#' # the fallback for.
#' max(abs(unlist(num) - unlist(distrib_grad_y_hess(d, y, theta))))
#'
#' # A family that does take the fallback, checked against numDeriv.
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


#' @title Default Hyperparameter Hessian of the Response Gradient
#' @name distrib_grad_y_hess.continuous_distrib
#'
#' @description
#' Falls back to one central difference of the analytic [distrib_cross_y()] in
#' each parameter, through [numerical_theta2_y()]. Registering the fallback on
#' `continuous_distrib` gives the third-order mixed derivative to every
#' continuous family, whether or not it writes one out.
#'
#' @details
#' The difference lands on an ANALYTIC quantity wherever the family provides
#' [distrib_cross_y()] in closed form, so the answer carries the error of one
#' stencil rather than two. Where that first-order quantity is itself a
#' fallback the two differences act on different variables and still compose
#' into a single mixed stencil.
#'
#' @param distrib A `continuous_distrib` object with no closed form of its own.
#' @param y A numeric vector of observations.
#' @param theta A named list of parameters.
#' @param scale One of `"parameter"` or `"link"`, applied by the generic before
#'   dispatch.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A named list with one numeric vector per unordered pair of
#'   parameters, keyed as [`hess_names(distrib@params)`][hess_names].
#'
#' @seealso [distrib_grad_y_hess()] for the generic,
#'   [distrib_hess_y_hess.continuous_distrib()] for the fourth-order twin, and
#'   [numerical_theta2_y()], which does the work.
#'
#' @keywords internal
#'
#' @examples
#' # The gamma writes no closed form, so this method answers for it.
#' d <- gamma2_distrib()
#' y <- c(0.5, 1, 2)
#' theta <- list(mu = 2, sigma2 = 1)
#' attr(S7::method(distrib_grad_y_hess, S7::S7_class(d)), "signature")[[1]]
#'
#' g <- distrib_grad_y_hess(d, y, theta)
#' c(g$mu_mu[1], g$mu_sigma2[1], g$sigma2_sigma2[1])
#'
#' # Against a numerical Hessian of the response gradient.
#' f <- function(v) distrib_grad_y(d, y[1], list(mu = v[1], sigma2 = v[2]))
#' numDeriv::hessian(f, c(2, 1))
S7::method(distrib_grad_y_hess, continuous_distrib) <- function(distrib, y, theta,
                                                                scale = c("parameter", "link"),
                                                                ...) {
  numerical_theta2_y(distrib, y, theta,
                     function(th) distrib_cross_y(distrib, y, th))
}

#' @title Default Hyperparameter Hessian of the Response Curvature
#' @name distrib_hess_y_hess.continuous_distrib
#'
#' @description
#' Falls back to one central difference of the analytic [distrib_cross2_y()] in
#' each parameter, through [numerical_theta2_y()]. It is
#' [distrib_grad_y_hess.continuous_distrib()] read one order higher in the
#' response, and it makes the fourth-order mixed derivative available for every
#' continuous family.
#'
#' @param distrib A `continuous_distrib` object with no closed form of its own.
#' @param y A numeric vector of observations.
#' @param theta A named list of parameters.
#' @param scale One of `"parameter"` or `"link"`, applied by the generic before
#'   dispatch.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A named list with one numeric vector per unordered pair of
#'   parameters, keyed as [`hess_names(distrib@params)`][hess_names].
#'
#' @seealso [distrib_hess_y_hess()] for the generic,
#'   [distrib_grad_y_hess.continuous_distrib()] for the third-order twin, and
#'   [numerical_theta2_y()], which does the work.
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
#' # Against a numerical Hessian of the response curvature.
#' f <- function(v) distrib_hess_y(d, y[1], list(mu = v[1], sigma2 = v[2]))
#' numDeriv::hessian(f, c(2, 1))
S7::method(distrib_hess_y_hess, continuous_distrib) <- function(distrib, y, theta,
                                                                scale = c("parameter", "link"),
                                                                ...) {
  numerical_theta2_y(distrib, y, theta,
                     function(th) distrib_cross2_y(distrib, y, th))
}


#' @title Gaussian Hyperparameter Hessian of the Response Gradient
#' @name distrib_grad_y_hess.Gaussian1Distrib
#'
#' @description
#' Closed form. With \eqn{r = y - \mu} the response gradient is
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
#' @param scale One of `"parameter"` or `"link"`, applied by the generic before
#'   dispatch.
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
#' # The three formulas written out.
#' r <- y - 0.4
#' c(mu_mu = 0, mu_sigma = -2 / 1.3^3, sigma_sigma = -6 * r[1] / 1.3^4)
#'
#' # Against a numerical Hessian of the response gradient.
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
#' Closed form, and nearly all of it zero. The gaussian's response curvature is
#' \eqn{-1/\sigma^2}, which carries no location at all, so
#' \deqn{\frac{\partial^4\ell}{\partial y^2\partial\sigma^2}
#'     = -\frac{6}{\sigma^4}}
#' is the only component that does not vanish, and it does not vary with the
#' data either.
#'
#' @param distrib A `Gaussian1Distrib` object.
#' @param y A numeric vector of observations. Only its length is used.
#' @param theta A named list with `mu` and `sigma`.
#' @param scale One of `"parameter"` or `"link"`, applied by the generic before
#'   dispatch.
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
#' # One surviving component, constant along y.
#' c(reported = h$sigma_sigma[1], formula = -6 / 1.3^4)
#'
#' # Against a numerical Hessian of the response curvature.
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
