#' @include negbin1_distrib.R
NULL

# NB1's third and fourth derivatives, written out.
#
# In the size r = mu/theta the log-likelihood is sparse:
#
#   l = G(r) + r B(theta) + C(theta),
#   G(r) = lgamma(y + r) - lgamma(r),  B = -log(1 + theta),
#   C = y log theta - y log(1 + theta),
#
# so the only composite piece is G(mu/theta). Its mixed derivatives follow a
# recursion that closes on itself. Writing
#
#   d^a/dmu^a d^b/dtheta^b G(mu/theta)
#       = theta^{-(a+b)} sum_j c_j r^j G^(a+j)(r),
#
# one more theta-derivative sends
#
#   c_j r^j G^(a+j)  ->  -(a + b + j) c_j r^j G^(a+j) - c_j r^{j+1} G^(a+j+1),
#
# because dr/dtheta = -r/theta contributes to both the power of r and the
# order of G. The coefficients are integers and the recursion is run rather
# than solved, which is what keeps every order exact with nothing transcribed
# beyond this one step.

#' Derivative Components of NB1
#'
#' @description
#' Returns the components of
#' \eqn{\partial^{a+b}\ell/\partial\mu^a\partial\theta^b} at any order from one
#' to four, from the sparse form the NB1 log-likelihood takes in the size
#' \eqn{r = \mu/\theta}.
#'
#' @details
#' # The sparse form, and the one composite piece
#'
#' Writing \eqn{r = \mu/\theta} the log-likelihood is
#' \deqn{\ell = G(r) + r B(\theta) + C(\theta), \qquad
#'       G(r) = \log\Gamma(y+r) - \log\Gamma(r),}
#' with \eqn{B = -\log(1+\theta)} and
#' \eqn{C = y\log\theta - y\log(1+\theta)}. The term \eqn{rB(\theta)} is
#' \eqn{\mu B(\theta)/\theta}, linear in \eqn{\mu}, so it contributes to
#' components carrying at most one \eqn{\mu}; \eqn{C} carries none. The only
#' composite piece is \eqn{G(\mu/\theta)}.
#'
#' # The recursion that closes on itself
#'
#' Its mixed derivatives take the form
#' \deqn{\frac{\partial^{a+b}}{\partial\mu^a\partial\theta^b}G(\mu/\theta)
#'   = \theta^{-(a+b)}\sum_j c_j\, r^j\, G^{(a+j)}(r),}
#' and one further \eqn{\theta}-derivative sends
#' \deqn{c_j r^j G^{(a+j)} \;\longrightarrow\;
#'       -(a+b+j)\,c_j r^j G^{(a+j)} - c_j r^{j+1} G^{(a+j+1)},}
#' because \eqn{\mathrm{d}r/\mathrm{d}\theta = -r/\theta} contributes to both
#' the power of \eqn{r} and the order of \eqn{G}. The coefficients are
#' integers, and the recursion is **run rather than solved**, so every order is
#' exact with nothing transcribed beyond this one step.
#'
#' # The cancellation the polygamma differences carry
#'
#' Each \eqn{G^{(m)}(r)} is a polygamma differenced at the shift \eqn{y}, which
#' is a count. As \eqn{\theta \to 0} the family tends to the Poisson,
#' \eqn{r = \mu/\theta} runs away, and the two terms of the difference agree to
#' leading order while the consumers above divide by \eqn{\theta^{a+b}}. The
#' differences therefore go through [psi_shift_diff()], which forms them as an
#' exact sum of reciprocals rather than as a subtraction. What that does not
#' repair is the cancellation among the powers of \eqn{r} in the recursion
#' itself: at orders three and four those terms are of size
#' \eqn{8\times10^6} at \eqn{\theta = 5\times10^{-4}} and sum to a value of
#' order one, so neither this form nor the one it replaced is reliable there.
#'
#' @param y A numeric vector of counts.
#' @param theta A named list with components `mu` and `theta`, each a numeric
#'   vector of length 1 or of the length of `y`, both strictly positive.
#'   Shorter components are recycled. Note that `theta` names both the list and
#'   its second component, the dispersion.
#' @param order The derivative order, an integer from 1 to 4.
#'
#' @return A named list of component vectors, one per distinct multi-index of
#'   the given order and keyed as [deriv_names()] keys them: two at order 1,
#'   three at order 2, four at order 3 and five at order 4. Each has the
#'   recycled length of the inputs.
#'
#' @seealso [distrib_deriv3.NegBin1Distrib()] and
#'   [distrib_deriv4.NegBin1Distrib()], which call this;
#'   [psi_shift_diff()] for the polygamma differences; and
#'   [negbin1_distrib()] for the family.
#' @keywords internal
negbin1_components <- function(y, theta, order) {
  mu <- theta[[1]]
  th <- theta[[2]]
  n <- max(length(y), lengths(theta[1:2]))
  y <- rep_len(y, n)
  mu <- rep_len(mu, n)
  th <- rep_len(th, n)
  r <- mu / th
  om <- 1 + th

  # G^(m)(r), m = 1..4.  Each is a polygamma differenced at the shift y,
  # which is a COUNT: as theta goes to zero the family tends to the Poisson,
  # r = mu/theta runs away and the two terms agree to leading order, while
  # the consumers below divide by theta^(a+b).  See psi_shift_diff().
  Gd <- lapply(1:4, function(m) psi_shift_diff(m - 1L, y, r))

  # B(theta) = -log(1 + theta) and its derivatives
  Bd <- function(m) {
    if (m == 0L) return(-base::log(om))
    (-1)^m * factorial(m - 1L) / om^m
  }
  # u(theta) = 1/theta and its derivatives
  ud <- function(m) (-1)^m * factorial(m) / th^(m + 1L)
  # M = B/theta by Leibniz, so that the term r B = mu M is linear in mu
  Md <- function(b) {
    s <- numeric(n)
    for (i in 0:b) s <- s + choose(b, i) * ud(i) * Bd(b - i)
    s
  }
  # C(theta) = y log theta - y log(1 + theta)
  Cd <- function(b) {
    if (b == 0L) return(y * base::log(th) - y * base::log(om))
    y * (-1)^(b - 1L) * factorial(b - 1L) * (1 / th^b - 1 / om^b)
  }

  # the G part, by the coefficient recursion
  gpart <- function(a, b) {
    coef <- c(1)                      # c_j indexed from j = 0
    for (step in seq_len(b)) {
      bb <- step - 1L                 # the b already applied
      new <- numeric(length(coef) + 1L)
      for (j in seq_along(coef)) {
        jj <- j - 1L
        new[j] <- new[j] - (a + bb + jj) * coef[j]
        new[j + 1L] <- new[j + 1L] - coef[j]
      }
      coef <- new
    }
    out <- numeric(n)
    for (j in seq_along(coef)) {
      if (coef[j] == 0) next
      jj <- j - 1L
      m <- a + jj
      gm <- if (m == 0L) lgamma(y + r) - lgamma(r) else Gd[[m]]
      out <- out + coef[j] * r^jj * gm
    }
    out / th^(a + b)
  }

  comp <- function(a, b) {
    out <- gpart(a, b)
    if (a == 0L) out <- out + mu * Md(b) + Cd(b)
    if (a == 1L) out <- out + Md(b)
    out
  }

  nms <- deriv_names(c("mu", "theta"), order)
  stats::setNames(lapply(nms, function(nm) {
    parts <- strsplit(nm, "_")[[1]]
    comp(sum(parts == "mu"), sum(parts == "theta"))
  }), nms)
}


#' @title NB1 Third-Order Derivatives
#' @name distrib_deriv3.NegBin1Distrib
#' @description
#' Computes the four distinct third derivatives of the NB1 log-likelihood in
#' the mean \eqn{\mu} and the dispersion \eqn{\theta}, **in closed form**,
#' through [negbin1_components()]. In the size \eqn{r = \mu/\theta} the
#' log-likelihood is \eqn{G(r) + rB(\theta) + C(\theta)}, so the only composite
#' piece is \eqn{G(\mu/\theta)}, and its mixed derivatives follow a recursion
#' in the powers of \eqn{r} and the order of \eqn{G} that is run rather than
#' solved.
#'
#' With `expected = TRUE` the method calls [expected_derivative()] instead: the
#' expected third derivatives have no closed form. That is the one place on
#' this page where `approx` and `nsim` are read.
#'
#' @param distrib A `NegBin1Distrib` object, from [negbin1_distrib()].
#' @param y A numeric vector of counts. With `expected = TRUE` only its length
#'   is read.
#' @param theta A named list with components `mu` and `theta`, each a numeric
#'   vector of length 1 or of the length of `y`, both strictly positive. A
#'   component of length 1 is recycled.
#' @param expected Logical of length 1. When `TRUE` the expectation under the
#'   model is returned in place of the value at the data, computed numerically.
#'   Defaults to `FALSE`.
#' @param scale One of `"parameter"` (the default) or `"link"`, matched by
#'   [base::match.arg()]. Read by the generic, not by this method.
#' @param approx One of `"integrate"` (the default here), `"bartlett"`, `"mc"`
#'   or `"opg"`. Read only when `expected = TRUE`.
#' @param nsim A single positive integer, the sample size when
#'   `approx = "mc"`. Read only when `expected = TRUE`. Defaults to `10000`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A named list of four numeric vectors, `mu_mu_mu`, `mu_mu_theta`,
#'   `mu_theta_theta` and `theta_theta_theta`, each of length
#'   `max(length(y), lengths(theta))`.
#'
#' @section Notation:
#' \eqn{\ell} is the log-likelihood of one observation, \eqn{\mu > 0} the mean,
#' \eqn{\theta > 0} the dispersion, \eqn{r = \mu/\theta} the negative binomial
#' size and \eqn{G(r) = \log\Gamma(y+r) - \log\Gamma(r)}.
#'
#' @section The Poisson boundary:
#' As \eqn{\theta \to 0} the family tends to the Poisson and the recursion's
#' terms in the powers of \eqn{r} grow while their sum stays of order one. The
#' polygamma differences go through [psi_shift_diff()] and are exact, but the
#' cancellation among those powers is not repaired: at
#' \eqn{\theta = 5\times10^{-4}} the terms reach \eqn{8\times10^{6}}, and this
#' order is not reliable below about \eqn{\theta = 0.05}. The score itself,
#' which does not carry the recursion, reaches the Poisson limit to five
#' figures.
#'
#' @seealso [distrib_hessian.NegBin1Distrib()] for the order below,
#'   [distrib_deriv4.NegBin1Distrib()] for the order above,
#'   [negbin1_components()] for the recursion,
#'   [distrib_deriv3.NegBin2Distrib()] for the other negative binomial, and
#'   [distrib_deriv3()] for the generic.
#'
#' @examples
#' d <- negbin1_distrib()
#' y <- c(0, 3, 7)
#' th <- list(mu = 4, theta = 1.2)
#' d3 <- distrib_deriv3(d, y, th)
#' names(d3)
#'
#' # A central difference of the Hessian reproduces the pure-mean component.
#' eps <- 1e-5
#' up <- distrib_hessian(d, y, list(mu = 4 + eps, theta = 1.2))$mu_mu
#' dn <- distrib_hessian(d, y, list(mu = 4 - eps, theta = 1.2))$mu_mu
#' all.equal((up - dn) / (2 * eps), d3$mu_mu_mu, tolerance = 1e-6)
#'
#' # And a mixed component, which is where the recursion does its work.
#' up <- distrib_hessian(d, y, list(mu = 4, theta = 1.2 + eps))$mu_theta
#' dn <- distrib_hessian(d, y, list(mu = 4, theta = 1.2 - eps))$mu_theta
#' all.equal((up - dn) / (2 * eps), d3$mu_theta_theta, tolerance = 1e-6)
S7::method(distrib_deriv3, NegBin1Distrib) <- function(distrib, y, theta,
                                                        expected = FALSE,
                                                        scale = c("parameter", "link"),
                                                        approx = c("integrate", "bartlett", "mc", "opg"),
                                                        nsim = 10000, ...) {
  if (expected) {
    return(expected_derivative(distrib, y, theta, order = 3L,
                               approx = match.arg(approx), nsim = nsim))
  }
  negbin1_components(y, theta, 3L)
}

#' @title NB1 Fourth-Order Derivatives
#' @name distrib_deriv4.NegBin1Distrib
#' @description
#' Computes the five distinct fourth derivatives of the NB1 log-likelihood in
#' \eqn{\mu} and \eqn{\theta}, **in closed form**, by the construction
#' [distrib_deriv3.NegBin1Distrib()] describes carried one order further: the
#' same coefficient recursion over the powers of \eqn{r = \mu/\theta} and the
#' order of \eqn{G(r) = \log\Gamma(y+r) - \log\Gamma(r)}.
#'
#' With `expected = TRUE` the method calls [expected_derivative()] instead: the
#' expected fourth derivatives have no closed form.
#'
#' @param distrib A `NegBin1Distrib` object, from [negbin1_distrib()].
#' @param y A numeric vector of counts. With `expected = TRUE` only its length
#'   is read.
#' @param theta A named list with components `mu` and `theta`, each a numeric
#'   vector of length 1 or of the length of `y`, both strictly positive. A
#'   component of length 1 is recycled.
#' @param expected Logical of length 1. When `TRUE` the expectation under the
#'   model is returned in place of the value at the data, computed numerically.
#'   Defaults to `FALSE`.
#' @param scale One of `"parameter"` (the default) or `"link"`, matched by
#'   [base::match.arg()]. Read by the generic, not by this method.
#' @param approx One of `"integrate"` (the default here), `"bartlett"`, `"mc"`
#'   or `"opg"`. Read only when `expected = TRUE`.
#' @param nsim A single positive integer, the sample size when
#'   `approx = "mc"`. Read only when `expected = TRUE`. Defaults to `10000`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A named list of five numeric vectors, `mu_mu_mu_mu`,
#'   `mu_mu_mu_theta`, `mu_mu_theta_theta`, `mu_theta_theta_theta` and
#'   `theta_theta_theta_theta`, each of length
#'   `max(length(y), lengths(theta))`.
#'
#' @section Notation:
#' \eqn{\ell} is the log-likelihood of one observation, \eqn{\mu > 0} the mean,
#' \eqn{\theta > 0} the dispersion, \eqn{r = \mu/\theta} the negative binomial
#' size and \eqn{G(r) = \log\Gamma(y+r) - \log\Gamma(r)}.
#'
#' @section The Poisson boundary:
#' The caveat of [distrib_deriv3.NegBin1Distrib()] applies here and more
#' strongly: the recursion divides by \eqn{\theta^{4}}, so the cancellation
#' among the powers of \eqn{r} is worse at this order than at the one below.
#' This order is not reliable below about \eqn{\theta = 0.05}.
#'
#' @seealso [distrib_deriv3.NegBin1Distrib()] for the order below and the
#'   recursion, [negbin1_components()] for the assembly, and
#'   [distrib_deriv4()] for the generic.
#'
#' @examples
#' d <- negbin1_distrib()
#' y <- c(0, 3, 7)
#' th <- list(mu = 4, theta = 1.2)
#' d4 <- distrib_deriv4(d, y, th)
#' names(d4)
#'
#' # A central difference of the third order reproduces a mixed component.
#' eps <- 1e-4
#' up <- distrib_deriv3(d, y, list(mu = 4, theta = 1.2 + eps))$mu_mu_theta
#' dn <- distrib_deriv3(d, y, list(mu = 4, theta = 1.2 - eps))$mu_mu_theta
#' all.equal((up - dn) / (2 * eps), d4$mu_mu_theta_theta, tolerance = 1e-5)
S7::method(distrib_deriv4, NegBin1Distrib) <- function(distrib, y, theta,
                                                        expected = FALSE,
                                                        scale = c("parameter", "link"),
                                                        approx = c("integrate", "bartlett", "mc", "opg"),
                                                        nsim = 10000, ...) {
  if (expected) {
    return(expected_derivative(distrib, y, theta, order = 4L,
                               approx = match.arg(approx), nsim = nsim))
  }
  negbin1_components(y, theta, 4L)
}
