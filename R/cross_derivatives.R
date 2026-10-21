#' @include distrib.R generics.R utility_functions.R y_derivatives.R gaussian1_distrib.R student_t1_distrib.R truncated.R
NULL

# Mixed second derivatives of the log-density: one derivative with respect to
# the response, one with respect to each parameter. Together with
# distrib_hess_y and distrib_hessian these complete the joint Hessian of a
# quantity that depends on both, which is what estimating the coefficients and
# the hyperparameters of a penalty together requires: the penalty is a negative
# log-density evaluated at the coefficients, and this block is its off-diagonal.
#
# Defined for continuous distributions, like the other response derivatives;
# a discrete distribution has no derivative in y to cross.

#' @title Mixed Response-Parameter Derivatives of the Log-Density
#'
#' @description
#' Computes the mixed second derivatives
#' \eqn{\partial^2 \ell / \partial y\, \partial \theta_i}, with one component
#' per parameter, and each component is a vector along `y`. Together with
#' [distrib_hess_y()] and [distrib_hessian()] this completes the joint Hessian
#' of the log-density in \eqn{(y, \theta)}: those two are the diagonal blocks
#' and this is the off-diagonal one.
#'
#' The function is defined only for continuous families. A discrete family has
#' no derivative in \eqn{y} to cross, so the package registers no method for it
#' and the call signals an error.
#'
#' @details
#' # What consumes it
#' A penalty in \pkg{penalties7} is a negative log-density read at the
#' coefficients, so estimating coefficients and hyperparameters together needs
#' exactly this block, and so does the gradient of a profiled objective through
#' the implicit function theorem.
#'
#' # The link scale
#' The component for \eqn{\eta_i} is the parameter-scale component multiplied
#' by \eqn{h_i'(\eta_i)}, and nothing else. The response derivative is
#' untouched by a reparametrization of \eqn{\theta}, so only the first-order
#' diagonal chain rule enters, exactly as for [distrib_gradient()]: on a
#' Gaussian whose scale carries a log link the `sigma` component is multiplied
#' by \eqn{\sigma} and the `mu` component by 1.
#'
#' # Where the closed forms are
#' Every continuous family in the package registers one. The fallback in
#' [distrib_cross_y.continuous_distrib()] therefore serves families defined
#' outside the package. It applies one tensor stencil to the highest quantity
#' that the family implements itself, through [mixed_tensor_derivatives()].
#'
#' @param distrib An object inheriting from `continuous_distrib`.
#' @param y A numeric vector of observations.
#' @param theta A named list, or a named numeric vector, of parameters. Each
#'   component must have length 1 or `length(y)`; a component of length 1 is
#'   recycled. The generic reorders `theta` by name and validates it against
#'   `params_bounds`, which are treated as open, before dispatch.
#' @param scale Either `"parameter"` (the default) or `"link"`, matched by
#'   [base::match.arg()]. The generic applies it after dispatch, so a method
#'   always returns the parameter scale and never reads this argument.
#' @param ... Passed to the method. The shipped methods do not read it.
#'
#' @return A named list with one numeric vector per parameter, keyed and
#'   ordered by `distrib@params`. Each vector has length `length(y)`.
#'
#' @examples
#' d <- gaussian1_distrib()
#' y <- c(-1, 0, 2)
#' th <- list(mu = 0.3, sigma = 1.4)
#' distrib_cross_y(d, y, th)
#'
#' # The mu component is 1/sigma^2 at every observation, and the sigma
#' # component is 2r/sigma^3, so it changes sign with the residual.
#' all.equal(distrib_cross_y(d, y, th)$sigma, 2 * (y - 0.3) / 1.4^3)
#'
#' # On the link scale only the diagonal chain-rule factor enters: sigma
#' # carries a log link, so its component is multiplied by sigma itself.
#' distrib_cross_y(d, y, th, scale = "link")$sigma /
#'   distrib_cross_y(d, y, th)$sigma
#'
#' @seealso [distrib_grad_y()] and [distrib_hess_y()], the two diagonal blocks;
#'   [numerical_cross_y()] for the fallback;
#'   [distrib_cross2_y()] for the third-order block.
#' @export
distrib_cross_y <- S7::new_generic("distrib_cross_y", "distrib", function(distrib, y, theta, scale = c("parameter", "link"), ...) {
  args <- check_derivative_args(distrib, y, theta)
  y <- args$y
  theta <- args$theta
  scale <- match.arg(scale)
  res <- S7::S7_dispatch()
  if (scale == "link") to_link_scale(distrib, theta, list(res), 1L) else res
})

#' @title Numerical Mixed Response-Parameter Derivatives
#'
#' @description
#' Computes \eqn{\partial^2 \ell / \partial y\, \partial \theta_i} by one
#' central difference of [distrib_grad_y()] in each parameter, at a step of
#' `h_rel` relative to the parameter's own magnitude. A family with a closed
#' form for some of its parameters calls it for the remaining ones through
#' `which`. The default [distrib_cross_y()] method for a continuous family
#' uses [mixed_tensor_derivatives()] instead.
#'
#' @details
#' # The quantity that is differenced
#' The quantity differenced is the response gradient, not the log-density,
#' so a family with an analytical `distrib_grad_y` requires exactly one
#' finite-difference layer. Where the response gradient is itself the
#' finite-difference fallback, the composition is the four-point mixed stencil
#' on the log-density: the two differences act on different variables, so they
#' commute into a single stencil instead of compounding as nested differences
#' in one variable do.
#'
#' # The step
#' \eqn{\varepsilon^{1/3} \approx 6.1\times10^{-6}} balances the \eqn{O(h^2)}
#' truncation of a central difference against a rounding term growing as
#' \eqn{1/h}, which is the usual optimum for a first derivative.
#'
#' @param distrib An object inheriting from `continuous_distrib`.
#' @param y A numeric vector of observations.
#' @param theta A named list of parameters, aligned to `distrib@params`.
#' @param h_rel The relative finite-difference step, a single positive number.
#'   Defaults to `.Machine$double.eps^(1/3)`. The step is taken by
#'   [fd_stable_step()], which keeps both evaluation points inside the bounds
#'   of the parameter.
#' @param which A character vector of parameter names to differentiate, or
#'   `NULL` (the default) for all of them. A family with a closed form for some
#'   of its parameters uses this argument, so that only the remaining ones cost
#'   a pair of evaluations.
#'
#' @return A named list with one numeric vector per requested parameter, in the
#'   order that `which` gives, or in the order of `distrib@params` when `which`
#'   is `NULL`. Each vector has length `length(y)`.
#'
#' @examples
#' d <- gaussian1_distrib()
#' y <- c(-1, 0, 1)
#' th <- list(mu = 0, sigma = 1)
#' numerical_cross_y(d, y, th)
#'
#' # It agrees with the family's own closed form to about 1e-10.
#' all.equal(numerical_cross_y(d, y, th), distrib_cross_y(d, y, th),
#'           tolerance = 1e-8)
#'
#' # 'which' costs a pair of evaluations only for the named parameter.
#' numerical_cross_y(d, y, th, which = "sigma")
#'
#' @seealso [distrib_cross_y()], the generic that it serves;
#'   [numerical_grad_y()] for the quantity that it differences;
#'   [numerical_cross2_y()] for the third-order analogue.
#' @export
numerical_cross_y <- function(distrib, y, theta, h_rel = .Machine$double.eps^(1 / 3),
                              which = NULL) {
  params <- distrib@params
  keep <- if (is.null(which)) params else which
  out <- vector("list", length(keep))
  names(out) <- keep
  for (j in match(keep, params)) {
    p <- params[j]
    quotient <- function(h) {
      th_up <- th_dn <- theta
      th_up[[j]] <- theta[[j]] + h
      th_dn[[j]] <- theta[[j]] - h
      (distrib_grad_y(distrib, y, th_up) -
        distrib_grad_y(distrib, y, th_dn)) / (2 * h)
    }
    out[[p]] <- fd_stable_step(quotient, theta[[j]],
                               distrib@params_bounds[[p]], h_rel)$value
  }
  out
}

#' @title Default Mixed Derivatives for Continuous Distributions
#' @name distrib_cross_y.continuous_distrib
#'
#' @description
#' The fallback applies one tensor stencil to the highest quantity that the
#' family implements itself, through [mixed_tensor_derivatives()]. A family
#' with an analytic response gradient or score differences it once, in the
#' parameter or in the response; a family with only a density takes one
#' four-point mixed stencil of the log-density.
#'
#' Every continuous family shipped in this package registers a closed form,
#' so this method serves a family defined outside the package, which obtains
#' the mixed block only from its density.
#'
#' @param distrib An object inheriting from `continuous_distrib` that
#'   registers no method of its own.
#' @param y A numeric vector of observations.
#' @param theta A named list of parameters, aligned by the generic.
#' @param scale Handled by the generic after dispatch; this method always
#'   returns the parameter scale.
#' @param ... Unused.
#'
#' @return A named list with one numeric vector per parameter, keyed by
#'   `distrib@params`. Each vector has length `length(y)`.
#'
#' @seealso [mixed_tensor_derivatives()], which does the work;
#'   [numerical_cross_y()] for the differenced part of a family's own method;
#'   [distrib_cross_y()] for the generic;
#'   [distrib_cross_y.Gaussian1Distrib()] for a closed form to compare against.
#' @keywords internal
S7::method(distrib_cross_y, continuous_distrib) <- function(distrib, y, theta,
                                                            scale = c("parameter", "link"),
                                                            ...) {
  mixed_tensor_derivatives(distrib, y, theta, 1L, 1L)
}

#' @title Gaussian Mixed Derivatives
#' @name distrib_cross_y.Gaussian1Distrib
#'
#' @description
#' Computes the closed form. With \eqn{r = y - \mu},
#' \deqn{\frac{\partial^2 \ell}{\partial y\, \partial \mu} = \frac{1}{\sigma^2},
#'       \qquad
#'       \frac{\partial^2 \ell}{\partial y\, \partial \sigma} =
#'       \frac{2r}{\sigma^3}.}
#'
#' The first component is constant in `y`, because the log-density is
#' quadratic in the response and in the location. The second is zero at the
#' mean and grows linearly with the residual, so it changes sign there.
#'
#' @param distrib A `Gaussian1Distrib` object.
#' @param y A numeric vector of observations.
#' @param theta A named list with components `mu` and `sigma`, aligned by the
#'   generic. Each may be length 1 or `length(y)`.
#' @param scale Handled by the generic after dispatch; this method always
#'   returns the parameter scale.
#' @param ... Unused.
#'
#' @return A named list with components `mu` and `sigma`. Each component is a
#'   numeric vector of length `length(y)`.
#'
#' @examples
#' d <- gaussian1_distrib()
#' y <- c(-1, 0, 2)
#' distrib_cross_y(d, y, list(mu = 0.3, sigma = 1.4))
#'
#' # The sigma component vanishes exactly at the mean.
#' distrib_cross_y(d, 0.3, list(mu = 0.3, sigma = 1.4))$sigma
#'
#' @seealso [distrib_cross_y()] for the generic;
#'   [distrib_cross_y.StudentT1Distrib()], whose mu component is not constant;
#'   [distrib_grad_y.Gaussian1Distrib()] for the quantity differentiated.
#' @keywords internal
S7::method(distrib_cross_y, Gaussian1Distrib) <- function(distrib, y, theta,
                                                         scale = c("parameter", "link"),
                                                         ...) {
  mu <- theta[[1]]
  sigma <- theta[[2]]
  r <- y - mu
  list(
    mu = rep(1 / sigma^2, length.out = length(y)),
    sigma = 2 * r / sigma^3
  )
}

#' @title Student's t Mixed Derivatives
#' @name distrib_cross_y.StudentT1Distrib
#'
#' @description
#' Computes the closed form in all three parameters. With \eqn{r = y - \mu} and
#' \eqn{D = \nu\sigma^2 + r^2},
#' \deqn{\frac{\partial^2 \ell}{\partial y\, \partial \mu} =
#'       \frac{(\nu+1)(\nu\sigma^2 - r^2)}{D^2}, \qquad
#'       \frac{\partial^2 \ell}{\partial y\, \partial \sigma} =
#'       \frac{2\nu\sigma(\nu+1) r}{D^2}, \qquad
#'       \frac{\partial^2 \ell}{\partial y\, \partial \nu} =
#'       -\frac{r(r^2 - \sigma^2)}{D^2}.}
#'
#' The \eqn{\nu} component is elementary here, unlike the \eqn{\nu} components
#' of the cdf derivatives: the log-density carries `lgamma` and a logarithm of
#' \eqn{D}, both differentiable in \eqn{\nu} in closed form, and the
#' expression does not involve a distribution function.
#'
#' The \eqn{\mu} component changes sign at \eqn{|r| = \sigma\sqrt\nu} and
#' decays as \eqn{r^{-2}}, where the Gaussian's is constant. It is the
#' derivative in \eqn{\mu} of the response score, which redescends.
#'
#' @param distrib A `StudentT1Distrib` object.
#' @param y A numeric vector of observations.
#' @param theta A named list with components `mu`, `sigma` and `nu`, aligned by
#'   the generic. Each may be length 1 or `length(y)`.
#' @param scale Handled by the generic after dispatch; this method always
#'   returns the parameter scale.
#' @param ... Unused.
#'
#' @return A named list with components `mu`, `sigma` and `nu`. Each component
#'   is a numeric vector of length `length(y)`.
#'
#' @examples
#' d <- student_t1_distrib()
#' th <- list(mu = 0.2, sigma = 1.1, nu = 6)
#' distrib_cross_y(d, c(-1, 0, 2), th)
#'
#' # The mu component changes sign where the residual passes sigma*sqrt(nu),
#' # which is where the score of a t stops growing and starts to redescend.
#' r0 <- 1.1 * sqrt(6)
#' distrib_cross_y(d, 0.2 + c(0.5, 1, 1.5) * r0, th)$mu
#'
#' @seealso [distrib_cross_y()] for the generic;
#'   [distrib_cross_y.Gaussian1Distrib()], the limit as \eqn{\nu} grows;
#'   [distrib_grad_y.StudentT1Distrib()] for the redescending score itself.
#' @keywords internal
S7::method(distrib_cross_y, StudentT1Distrib) <- function(distrib, y, theta,
                                                         scale = c("parameter", "link"),
                                                         ...) {
  student_t1_cross_y_cpp(y, theta[[1]], theta[[2]], theta[[3]])
}

#' @title Mixed Derivatives of a Truncated Distribution
#' @name distrib_cross_y.TruncatedContinuousDistrib
#'
#' @description
#' Returns the parent's mixed derivatives unchanged. Truncation replaces the
#' parent's log-density by \eqn{\ell(y;\theta) - \log Z(\theta)}, and the
#' normalizing constant does not depend on \eqn{y}, so it vanishes from any
#' derivative taking one derivative in the response. The equality is exact, not
#' approximate, and holds at every observation inside the support.
#'
#' Like the response gradient and the response Hessian of a truncated family,
#' this block costs nothing beyond the parent's. The gradient and the Hessian
#' in the parameters both need \eqn{\log Z} and its derivatives, which cost one
#' quadrature or one pair of cdf evaluations per component.
#'
#' @param distrib A `TruncatedContinuousDistrib` object.
#' @param y A numeric vector of observations inside the truncation interval.
#'   Outside it the parent's value is returned, because the method performs
#'   no test.
#' @param theta A named list of parameters of the parent, aligned by the
#'   generic. Truncation adds no parameter, because the endpoints are constants.
#' @param scale Handled by the generic after dispatch; this method always
#'   returns the parameter scale, and passes `"parameter"` down so the chain
#'   rule is applied once rather than twice.
#' @param ... Passed to the parent's method.
#'
#' @return A named list with one numeric vector per parameter, keyed by the
#'   parent's `params`. Each vector has length `length(y)`.
#'
#' @examples
#' d <- gaussian1_distrib()
#' tr <- truncated(d, lower = -2, upper = 3)
#' y <- c(-1, 0, 2)
#' th <- list(mu = 0.3, sigma = 1.4)
#'
#' # The result is identical to the parent's, component for component.
#' all.equal(distrib_cross_y(tr, y, th), distrib_cross_y(d, y, th))
#'
#' # The gradient is not: there log Z does depend on theta.
#' all.equal(distrib_gradient(tr, y, th), distrib_gradient(d, y, th))
#'
#' @seealso [truncated()] for the wrapper;
#'   [distrib_cross_y()] for the generic;
#'   [distrib_gradient()], which does pay for the normalizing constant.
#' @keywords internal
S7::method(distrib_cross_y, TruncatedContinuousDistrib) <- function(distrib, y, theta,
                                                                    scale = c("parameter", "link"),
                                                                    ...) {
  distrib_cross_y(distrib@parent_distrib, y, theta, scale = "parameter", ...)
}
