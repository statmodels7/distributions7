#' @title Differences of a Polygamma at an Integer Shift
#' @name psi_shift_diff
#' @description
#' Returns \eqn{\psi^{(n)}(x+k) - \psi^{(n)}(x)} without the cancellation the
#' direct difference carries. It is the R twin of the \code{psi_A_rest()} and
#' \code{psi_T_rest()} of \code{src/psi_diff.h}, written for any order rather
#' than for the two the compiled kernels need.
#'
#' @details
#' Several families here carry a boundary limit at which a shape runs to
#' infinity and the family tends to a simpler one. Every derivative in that
#' parameter vanishes there, so it is written as a difference of terms that
#' agree to leading order: both \eqn{\psi^{(n)}(x+k)} and \eqn{\psi^{(n)}(x)}
#' are of size \eqn{x^{-n}} while their difference is of size \eqn{k x^{-n-1}},
#' and the direct subtraction therefore loses a digit for every factor of ten
#' in \eqn{x/k}. Consumers that divide the result by a power of the parameter
#' amplify what is left.
#'
#' Above \eqn{x = 100} the value comes from the asymptotic expansion
#' \deqn{\psi^{(n)}(z) \sim (-1)^{n-1}\Big[\frac{(n-1)!}{z^{n}}
#'   + \frac{n!}{2 z^{n+1}} + \frac{(n+1)!}{12 z^{n+2}}
#'   - \frac{(n+3)!}{720 z^{n+4}} + \frac{(n+5)!}{30240 z^{n+6}}\Big],}
#' whose \eqn{n = 0} case is \eqn{\log z - 1/(2z) - 1/(12 z^{2}) + \ldots}
#' with the logarithm differenced as \eqn{\log(1 + k/x)}. Each power is
#' differenced as
#' \deqn{\frac{1}{(x+k)^{p}} - \frac{1}{x^{p}}
#'   = \frac{1}{x^{p}}\,\big(e^{-p\log(1+k/x)} - 1\big),}
#' which \code{expm1} and \code{log1p} evaluate to the last bit however small
#' \eqn{k/x} is. Below the crossover the direct difference has all its digits
#' and is used as it stands.
#'
#' @param n The order of the polygamma: \code{0} for \code{digamma},
#'   \code{1} for \code{trigamma}, and so on.
#' @param k The shift, which for every caller here is a count or a size and
#'   therefore a non-negative integer. Recycled against \code{x}.
#' @param x The argument, strictly positive. Recycled against \code{k}.
#' @return A numeric vector of the same length as the recycled arguments.
#'   Exactly zero wherever \code{k} is zero, by either branch.
#' @keywords internal
psi_shift_diff <- function(n, k, x) {
  n <- as.integer(n)
  m <- max(length(k), length(x))
  k <- rep_len(k, m)
  x <- rep_len(x, m)
  out <- numeric(m)
  big <- x >= 100

  if (any(!big)) {
    i <- which(!big)
    out[i] <- if (n == 0L) {
      digamma(x[i] + k[i]) - digamma(x[i])
    } else {
      psigamma(x[i] + k[i], n) - psigamma(x[i], n)
    }
  }

  if (any(big)) {
    i <- which(big)
    a <- x[i]
    l1 <- log1p(k[i] / a)
    ## 1/(x+k)^p - 1/x^p, with no power of x formed on its own where it
    ## would overflow: past that point the quotient is zero, which is the
    ## limit of a term that has long stopped contributing.
    dp <- function(p) expm1(-p * l1) / a^p
    out[i] <- if (n == 0L) {
      l1 - 0.5 * dp(1) - dp(2) / 12 + dp(4) / 120 - dp(6) / 252
    } else {
      (-1)^(n - 1) * (
        factorial(n - 1) * dp(n) +
          factorial(n) / 2 * dp(n + 1) +
          factorial(n + 1) / 12 * dp(n + 2) -
          factorial(n + 3) / 720 * dp(n + 4) +
          factorial(n + 5) / 30240 * dp(n + 6)
      )
    }
  }
  out
}
