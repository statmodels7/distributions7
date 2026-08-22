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
#' The components of \eqn{\partial^{a+b}\ell/\partial\mu^a\partial\theta^b} at
#' any order from one to four, from the sparse form of the log-likelihood in
#' the size \eqn{r = \mu/\theta}.
#'
#' @param y A numeric vector of counts.
#' @param theta A list containing \code{mu} and \code{theta}.
#' @param order The derivative order, 1 to 4.
#'
#' @return A named list of component vectors, keyed as
#'   \code{\link{deriv_names}}.
#'
#' @seealso \code{\link{negbin1_distrib}}
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


#' @title NB1 Third and Fourth Derivatives
#' @name distrib_deriv3.NegBin1Distrib
#' @description
#' Closed form at both orders, from \code{\link{negbin1_components}}. In the
#' size \eqn{r = \mu/\theta} the log-likelihood is
#' \eqn{G(r) + rB(\theta) + C(\theta)}, so the only composite piece is
#' \eqn{G(\mu/\theta)} and its mixed derivatives follow a recursion in the
#' powers of \eqn{r} and the order of \eqn{G}.
#' @param distrib A \code{NegBin1Distrib} object.
#' @param y A numeric vector of counts.
#' @param theta A list containing \code{mu} and \code{theta}.
#' @param expected Logical; if \code{TRUE}, the expected derivatives.
#' @param scale Either \code{"parameter"} or \code{"link"}; handled by the generic.
#' @param approx The approximation used when \code{expected} is \code{TRUE}.
#' @param nsim Monte Carlo draws when \code{approx = "mc"}.
#' @param ... Unused.
#' @return A named list of third-derivative components.
#' @seealso \code{\link{negbin1_distrib}}
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

#' @rdname distrib_deriv3.NegBin1Distrib
#' @name distrib_deriv4.NegBin1Distrib
#' @return A named list of fourth-derivative components.
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
