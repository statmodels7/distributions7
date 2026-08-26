#' @include distrib.R generics.R
NULL

# Derivatives of the log-density with respect to the response y (not the
# parameters). Registered on continuous_distrib; distributions with a closed
# form override these. For discrete distributions the derivative w.r.t. y is not
# defined, so no default method is provided.

#' @title Finite-Difference Steps That Respect the Support
#'
#' @description
#' Returns the step for a central difference in the RESPONSE: `h_rel` scaled by
#' \eqn{\max(1, |y|)}, then shrunk so that \eqn{y \pm h} stays strictly inside
#' the distribution's support. It is the response counterpart of [fd_steps()],
#' and both numerical response derivatives take their step from it.
#'
#' @details
#' The scaling by \eqn{\max(1, |y|)} makes the step relative where the response
#' is large and absolute where it is small, so a value near zero is not
#' differenced with a step below the resolution of a double.
#'
#' The clamp is what the support requires. A gamma observation at
#' \eqn{y = 10^{-3}} differenced with the default step of
#' \eqn{6 \times 10^{-6}} needs no help, but one at \eqn{y = 10^{-8}} would be
#' evaluated at a negative point, where the density is not defined and
#' `distrib_pdf()` returns `-Inf`. The factor 0.49 leaves the step under half
#' the distance to the bound, so both evaluation points stay inside with room
#' to spare: at \eqn{y = 10^{-8}} the step becomes \eqn{4.9\times 10^{-9}} and
#' the left point \eqn{5.1\times 10^{-9}}.
#'
#' @param y A numeric vector of observations.
#' @param bounds A numeric vector of length two, the distribution's support.
#'   An infinite endpoint imposes no clamp on that side.
#' @param h_rel The relative step size, a single positive number. The callers
#'   pass \eqn{\varepsilon^{1/3}} at first order and \eqn{\varepsilon^{1/4}} at
#'   second.
#'
#' @return A numeric vector of steps, as long as `y`, every entry positive.
#'
#' @section Notation:
#' \eqn{\ell} is the log-density of one observation, \eqn{y} the response,
#' \eqn{h} the finite-difference step and \eqn{\varepsilon} the machine
#' epsilon, `.Machine$double.eps`.
#'
#' @seealso [fd_steps()] for the parameter counterpart, and
#'   [numerical_grad_y()] and [numerical_hess_y()], its two callers.
#'
#' @keywords internal
#'
#' @examples
#' # On the whole line the step is h_rel times max(1, |y|).
#' h_rel <- .Machine$double.eps^(1 / 3)
#' distributions7:::fd_steps_y(c(-100, 0, 0.5, 100), c(-Inf, Inf), h_rel)
#' h_rel * pmax(1, abs(c(-100, 0, 0.5, 100)))
#'
#' # Near a bound the step is cut to 0.49 of the distance to it, so the left
#' # evaluation point stays inside the support.
#' y <- c(1e-8, 1e-3, 1)
#' h <- distributions7:::fd_steps_y(y, c(0, Inf), h_rel)
#' rbind(step = h, left_point = y - h)
#'
#' # The two steps the callers use differ by a factor of twenty.
#' c(first_order = .Machine$double.eps^(1 / 3),
#'   second_order = .Machine$double.eps^(1 / 4))
fd_steps_y <- function(y, bounds, h_rel) {
  h <- h_rel * pmax(1, abs(y))
  if (is.finite(bounds[1])) h <- pmin(h, 0.49 * (y - bounds[1]))
  if (is.finite(bounds[2])) h <- pmin(h, 0.49 * (bounds[2] - y))
  h
}

#' @title Numerical Gradient of the Log-Density with Respect to the Response
#'
#' @description
#' Computes \eqn{\partial \ell / \partial y} by the two-point central
#' difference
#' \deqn{\frac{\partial\ell}{\partial y} \approx
#'   \frac{\ell(y+h) - \ell(y-h)}{2h},}
#' evaluated on `distrib_pdf(..., log = TRUE)`. It is what the default
#' [distrib_grad_y()] method runs for a continuous family with no closed form.
#'
#' @details
#' # The step and the accuracy it buys
#'
#' The stencil's truncation error is \eqn{O(h^2)} and its rounding error
#' \eqn{O(\varepsilon/h)}, so the total is smallest at
#' \eqn{h \sim \varepsilon^{1/3}}, which is the default. Measured against the
#' gamma's own closed form at `h_rel` from \eqn{10^{-1}} to \eqn{10^{-8}}, the
#' error is 8.2e-02, 8.0e-04, 8.0e-06, 8.0e-08, 8.0e-10, 1.4e-10, 8.2e-10 and
#' 2.0e-08: it divides by a hundred per decade while the truncation dominates,
#' turns near the default and rises again as the rounding takes over. A step
#' ten times shorter than the default is therefore worse, not better.
#'
#' The step itself comes from [fd_steps_y()], which scales it by
#' \eqn{\max(1, |y|)} and shrinks it near a finite bound.
#'
#' # The cost
#'
#' Two evaluations of the log-density, both vectorized over `y`, so one call
#' costs two `distrib_pdf()` calls whatever the sample size. A family with a
#' closed form pays neither.
#'
#' @param distrib An object inheriting from `continuous_distrib`. A discrete
#'   family has no derivative in the response and registers no method.
#' @param y A numeric vector of observations.
#' @param theta A named list of parameters.
#' @param h_rel The relative step, defaulting to `.Machine$double.eps^(1/3)`,
#'   the exponent that minimizes the total error of a central first difference.
#'
#' @return A numeric vector as long as `y`.
#'
#' @section Notation:
#' \eqn{\ell} is the log-density of one observation, \eqn{y} the response,
#' \eqn{h} the finite-difference step and \eqn{\varepsilon} the machine
#' epsilon, `.Machine$double.eps`.
#'
#' @seealso [distrib_grad_y()] for the generic, [numerical_hess_y()] for the
#'   second order, and [fd_steps_y()] for the step rule.
#'
#' @examples
#' d <- gamma2_distrib()
#' y <- c(0.5, 1, 2)
#' theta <- list(mu = 2, sigma2 = 1)
#'
#' numerical_grad_y(d, y, theta)
#'
#' # Against the gamma's own closed form: the difference is good to about
#' # 3e-10, not to machine precision.
#' distrib_grad_y(d, y, theta)
#' max(abs(numerical_grad_y(d, y, theta) - distrib_grad_y(d, y, theta)))
#'
#' # Too short a step is worse than too long. The error falls as h^2 and then
#' # rises as eps / h, and the default sits at the turn.
#' err <- function(p) max(abs(numerical_grad_y(d, y, theta, h_rel = 10^(-p)) -
#'                            distrib_grad_y(d, y, theta)))
#' setNames(vapply(1:8, err, numeric(1)), paste0("1e-", 1:8))
#'
#' @export
numerical_grad_y <- function(distrib, y, theta, h_rel = .Machine$double.eps^(1 / 3)) {
  h <- fd_steps_y(y, distrib@bounds, h_rel)
  (distrib_pdf(distrib, y + h, theta, log = TRUE) -
    distrib_pdf(distrib, y - h, theta, log = TRUE)) / (2 * h)
}

#' @title Numerical Second Response Derivative of the Log-Density
#'
#' @description
#' Computes \eqn{\partial^2 \ell / \partial y^2} by the three-point central
#' stencil
#' \deqn{\frac{\partial^2\ell}{\partial y^2} \approx
#'   \frac{\ell(y+h) - 2\ell(y) + \ell(y-h)}{h^2},}
#' evaluated on `distrib_pdf(..., log = TRUE)`. It is what the default
#' [distrib_hess_y()] method runs for a continuous family with no closed form.
#'
#' @details
#' # Why the step is longer and the answer coarser
#'
#' The truncation error is again \eqn{O(h^2)}, but the rounding error is
#' \eqn{O(\varepsilon/h^2)} where the first order's is
#' \eqn{O(\varepsilon/h)}, the numerator being a difference of quantities of
#' the same size divided by \eqn{h^2}. The
#' total is smallest at \eqn{h \sim \varepsilon^{1/4}}, which is the default
#' and is twenty times the first order's step, and the accuracy attainable is
#' correspondingly coarser. Measured against the gamma's own closed form at
#' `h_rel` from \eqn{10^{-1}} to \eqn{10^{-8}}: 2.5e-01, 2.4e-03, 2.4e-05,
#' 2.8e-07, 5.4e-06, 6.2e-04, 3.5e-02 and 3.1e+00. The rise past the optimum is
#' far steeper than the first order's, which is the \eqn{h^{-2}} at work.
#'
#' # The cost
#'
#' Three evaluations of the log-density, one of them at `y` itself, all
#' vectorized over `y`.
#'
#' @param distrib An object inheriting from `continuous_distrib`. A discrete
#'   family has no derivative in the response and registers no method.
#' @param y A numeric vector of observations.
#' @param theta A named list of parameters.
#' @param h_rel The relative step, defaulting to `.Machine$double.eps^(1/4)`,
#'   the exponent that minimizes the total error of a central second
#'   difference.
#'
#' @return A numeric vector as long as `y`.
#'
#' @section Notation:
#' \eqn{\ell} is the log-density of one observation, \eqn{y} the response,
#' \eqn{h} the finite-difference step and \eqn{\varepsilon} the machine
#' epsilon, `.Machine$double.eps`.
#'
#' @seealso [distrib_hess_y()] for the generic, [numerical_grad_y()] for the
#'   first order, and [fd_steps_y()] for the step rule.
#'
#' @examples
#' d <- gamma2_distrib()
#' y <- c(0.5, 1, 2)
#' theta <- list(mu = 2, sigma2 = 1)
#'
#' numerical_hess_y(d, y, theta)
#' distrib_hess_y(d, y, theta)
#'
#' # The second difference is the cruder of the two, by about three orders of
#' # magnitude on the same family and the same points.
#' c(order2 = max(abs(numerical_hess_y(d, y, theta) -
#'                    distrib_hess_y(d, y, theta))),
#'   order1 = max(abs(numerical_grad_y(d, y, theta) -
#'                    distrib_grad_y(d, y, theta))))
#'
#' # And it degrades far faster below its optimal step.
#' err <- function(p) max(abs(numerical_hess_y(d, y, theta, h_rel = 10^(-p)) -
#'                            distrib_hess_y(d, y, theta)))
#' setNames(vapply(3:8, err, numeric(1)), paste0("1e-", 3:8))
#'
#' @export
numerical_hess_y <- function(distrib, y, theta, h_rel = .Machine$double.eps^(1 / 4)) {
  h <- fd_steps_y(y, distrib@bounds, h_rel)
  (distrib_pdf(distrib, y + h, theta, log = TRUE) -
    2 * distrib_pdf(distrib, y, theta, log = TRUE) +
    distrib_pdf(distrib, y - h, theta, log = TRUE)) / (h^2)
}

#' @title Default Response Gradient for Continuous Distributions
#' @name distrib_grad_y.continuous_distrib
#'
#' @description
#' Computes \eqn{\partial \ell / \partial y} by the two-point central
#' difference \eqn{\{\ell(y+h) - \ell(y-h)\}/(2h)} on
#' `distrib_pdf(..., log = TRUE)`, through [numerical_grad_y()]. Registering
#' the fallback on `continuous_distrib` gives a response gradient to a family
#' that supplies a density and nothing else.
#'
#' @details
#' The step is \eqn{h = \varepsilon^{1/3}\max(1, |y|)}, shrunk near a finite
#' bound by [fd_steps_y()] so that both evaluation points stay inside the
#' support. The stencil's truncation error is \eqn{O(h^2)} and its rounding
#' error \eqn{O(\varepsilon/h)}, so that exponent minimizes their sum;
#' measured against the gamma's closed form it delivers about
#' \eqn{3\times 10^{-10}}. The cost is two vectorized `distrib_pdf()` calls.
#'
#' EVERY family the package ships writes its own response gradient, so this
#' method answers for user-defined families alone. A DISCRETE family has no
#' derivative in the response, so no default is registered there and the call
#' raises with the class named.
#'
#' @param distrib A `continuous_distrib` object with no closed form of its own.
#' @param y A numeric vector of observations.
#' @param theta A named list of parameters.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector as long as `y`.
#'
#' @section Notation:
#' \eqn{\ell} is the log-density of one observation, \eqn{y} the response,
#' \eqn{h} the finite-difference step and \eqn{\varepsilon} the machine
#' epsilon, `.Machine$double.eps`.
#'
#' @seealso [distrib_grad_y()] for the generic, [numerical_grad_y()], which
#'   does the work, and [distrib_hess_y.continuous_distrib()] for the second
#'   order.
#'
#' @keywords internal
#'
#' @examples
#' # A family that supplies a density and nothing else, which is what the
#' # fallback exists for: every family the package ships writes its own
#' # response derivatives.
#' Toy <- S7::new_class("Toy", parent = continuous_distrib)
#' S7::method(distrib_pdf, Toy) <- function(distrib, y, theta, ...) {
#'   z <- (y - theta$mu) / theta$sigma
#'   v <- -log(2 * theta$sigma) - abs(z) - 0.1 * z^2
#'   if (isTRUE(list(...)$log)) v else exp(v)
#' }
#' d <- Toy(distrib_name = "toy", dimension = "univariate",
#'          bounds = c(-Inf, Inf), params = c("mu", "sigma"),
#'          params_interpretation = c("location", "scale"), n_params = 2L,
#'          params_bounds = list(mu = c(-Inf, Inf), sigma = c(0, Inf)),
#'          link_params = list(mu = linkfunctions7::identity_link(),
#'                             sigma = linkfunctions7::log_link()))
#' y <- c(-1, 0.4, 1)
#' theta <- list(mu = 0.3, sigma = 1.2)
#'
#' attr(S7::method(distrib_grad_y, Toy), "signature")[[1]]
#' distrib_grad_y(d, y, theta)
#'
#' # Against a Richardson-extrapolated derivative of the same log-density,
#' # which shares no arithmetic with the stencil.
#' f <- function(v) distrib_pdf(d, v, theta, log = TRUE)
#' vapply(y, function(v) numDeriv::grad(f, v), numeric(1))
#'
#' # A discrete family has no method at all.
#' try(distrib_grad_y(poisson_distrib(), 1:3, list(mu = 2)))
S7::method(distrib_grad_y, continuous_distrib) <- function(distrib, y, theta, ...) {
  numerical_grad_y(distrib, y, theta)
}

#' @title Default Response Hessian for Continuous Distributions
#' @name distrib_hess_y.continuous_distrib
#'
#' @description
#' Computes \eqn{\partial^2 \ell / \partial y^2} by the three-point central
#' stencil \eqn{\{\ell(y+h) - 2\ell(y) + \ell(y-h)\}/h^2} on
#' `distrib_pdf(..., log = TRUE)`, through [numerical_hess_y()]. It is one
#' stencil on the log-density and never a difference of
#' [distrib_grad_y.continuous_distrib()], which would be two layers of rounding
#' on the same variable.
#'
#' @details
#' The step is \eqn{h = \varepsilon^{1/4}\max(1, |y|)}, shrunk near a finite
#' bound by [fd_steps_y()]. The longer exponent is what the second difference
#' needs: its rounding error is \eqn{O(\varepsilon/h^2)} where the first
#' order's is \eqn{O(\varepsilon/h)}, so the optimum sits at a step twenty
#' times longer and the accuracy reached is coarser, measured at about
#' \eqn{4\times 10^{-7}} against the gamma's closed form where the first order
#' reaches \eqn{3\times 10^{-10}}. The cost is three vectorized
#' `distrib_pdf()` calls.
#'
#' As at first order, every family the package ships writes its own, and a
#' DISCRETE family registers no default.
#'
#' @param distrib A `continuous_distrib` object with no closed form of its own.
#' @param y A numeric vector of observations.
#' @param theta A named list of parameters.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector as long as `y`.
#'
#' @section Notation:
#' \eqn{\ell} is the log-density of one observation, \eqn{y} the response,
#' \eqn{h} the finite-difference step and \eqn{\varepsilon} the machine
#' epsilon, `.Machine$double.eps`.
#'
#' @seealso [distrib_hess_y()] for the generic, [numerical_hess_y()], which
#'   does the work, and [distrib_grad_y.continuous_distrib()] for the first
#'   order.
#'
#' @keywords internal
#'
#' @examples
#' # A family that supplies a density and nothing else, which is what the
#' # fallback exists for: every family the package ships writes its own
#' # response derivatives.
#' Toy <- S7::new_class("Toy", parent = continuous_distrib)
#' S7::method(distrib_pdf, Toy) <- function(distrib, y, theta, ...) {
#'   z <- (y - theta$mu) / theta$sigma
#'   v <- -log(2 * theta$sigma) - abs(z) - 0.1 * z^2
#'   if (isTRUE(list(...)$log)) v else exp(v)
#' }
#' d <- Toy(distrib_name = "toy", dimension = "univariate",
#'          bounds = c(-Inf, Inf), params = c("mu", "sigma"),
#'          params_interpretation = c("location", "scale"), n_params = 2L,
#'          params_bounds = list(mu = c(-Inf, Inf), sigma = c(0, Inf)),
#'          link_params = list(mu = linkfunctions7::identity_link(),
#'                             sigma = linkfunctions7::log_link()))
#' y <- c(-1, 0.4, 1)
#' theta <- list(mu = 0.3, sigma = 1.2)
#'
#' distrib_hess_y(d, y, theta)
#'
#' # Against a Richardson-extrapolated second derivative of the log-density.
#' f <- function(v) distrib_pdf(d, v, theta, log = TRUE)
#' vapply(y, function(v) numDeriv::hessian(f, v)[1, 1], numeric(1))
#'
#' # A discrete family has no method at all.
#' try(distrib_hess_y(poisson_distrib(), 1:3, list(mu = 2)))
S7::method(distrib_hess_y, continuous_distrib) <- function(distrib, y, theta, ...) {
  numerical_hess_y(distrib, y, theta)
}
