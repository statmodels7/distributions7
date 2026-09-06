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
#' exact sum of reciprocals rather than as a subtraction.
#'
#' # Where the recursion cedes, and the form that does not
#'
#' What that does not repair is the cancellation among the powers of \eqn{r} in
#' the recursion itself: at orders three and four those terms are of size
#' \eqn{8\times10^6} at \eqn{\theta = 5\times10^{-4}} and sum to a value of
#' order one. Measured at \eqn{\mu = 4}, \eqn{y = 3}, the third derivative in
#' \eqn{\theta} reads \eqn{-1.97\times10^{9}} at \eqn{\theta = 10^{-6}} where
#' the value is \eqn{9/32}, and the fourth \eqn{7.90\times10^{15}}.
#'
#' The size need not appear at all, which is what removes it. Since
#' \eqn{\log\Gamma(y+r) - \log\Gamma(r) = \sum_{i<y}\log(r+i)} and
#' \eqn{r = \mu/\theta}, the \eqn{y\log\theta} each such term carries cancels
#' EXACTLY against \eqn{C(\theta)}, leaving
#' \deqn{\ell = \sum_{i<y}\log(\mu + i\theta) - \log(y!)
#'   - \frac{\mu}{\theta}\log(1+\theta) - y\log(1+\theta),}
#' whose variable does not run away and whose every derivative is elementary.
#' Below [nb1_exact_cut()] that is the route taken, through
#' [nb1_components_exact()]; above it the recursion is within 4e-11 and costs
#' \eqn{O(1)} where the sum costs \eqn{y} terms an observation.
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

  # Below the crossover the route in the size loses the answer among the
  # powers of r; the form in mu + i theta does not, and costs y terms an
  # observation. Measured at mu = 4, y = 3, the two agree to between 4e-16 and
  # 5.5e-09 down to theta = 0.1 and part company below it: at theta = 1e-3 the
  # fourth derivative reads 8.97 against -1.59, and at 1e-8 5.8e+25 against
  # -1.598437. The threshold is 1, where the agreement is 4e-11 or better.
  esatta <- th < nb1_exact_cut()
  if (all(esatta)) return(nb1_components_exact(y, mu, th, order))
  if (any(esatta)) {
    dentro <- nb1_components_exact(y[esatta], mu[esatta], th[esatta], order)
    fuori <- negbin1_components(y[!esatta],
                                list(mu = mu[!esatta], theta = th[!esatta]),
                                order)
    return(stats::setNames(lapply(names(fuori), function(nm) {
      v <- numeric(n)
      v[esatta] <- dentro[[nm]]
      v[!esatta] <- fuori[[nm]]
      v
    }), names(fuori)))
  }

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

#' @rdname negbin1_components
#' @description
#' `nb1_exact_cut()` is the dispersion below which the cancellation-free
#' assembly is used in place of the recursion in the size. It sits where the
#' two agree and each is still comfortable: measured at \eqn{\mu = 4},
#' \eqn{y = 3}, they agree to 4e-11 or better at \eqn{\theta = 1} and to
#' 5.5e-09 at \eqn{\theta = 0.1}, and part company below that.
#' @keywords internal
nb1_exact_cut <- function() 1

#' @rdname negbin1_components
#' @description
#' `nb1_M_derivs()` returns \eqn{M(\theta) = \log(1+\theta)/\theta} and its
#' derivatives to the order asked for. It is the one composite piece of the
#' cancellation-free form below, and it has a removable singularity at
#' \eqn{\theta = 0}: the recursion \eqn{\theta M^{(b+1)} + (b+1)M^{(b)} =
#' (-1)^b b!/(1+\theta)^{b+1}} divides by \eqn{\theta} and loses its digits
#' there, while the series \eqn{M^{(b)} = (-1)^b \sum_{n\ge 0} (-1)^n
#' \left[\prod_{i\le b}(n+i)\right]\theta^n/(n+b+1)} converges only below one.
#' The crossover is MEASURED and not chosen: the two agree to between 1.9e-16
#' and 3.4e-13 over \eqn{\theta} from 0.01 to 0.5, and each fails on its own
#' side -- the recursion by 1.4e-04 at \eqn{\theta = 10^{-6}} and by 2.3e+08 at
#' order four, the series by 1.1e-02 at \eqn{\theta = 0.8}. It is the shape the
#' generalized Pareto's `Lambda` already carries one family over.
#' @keywords internal
nb1_M_derivs <- function(th, order, cut = 0.5, nterm = 80L) {
  out <- vector("list", order + 1L)
  low <- th < cut
  nn <- 0:nterm
  for (b in 0:order) {
    v <- numeric(length(th))
    if (any(low)) {
      co <- rep(1, length(nn))
      for (i in seq_len(b)) co <- co * (nn + i)
      w <- (-1)^b * ((-1)^nn * co / (nn + b + 1L))
      # una potenza di theta per termine, accumulata per righe
      acc <- numeric(sum(low)); tl <- th[low]; p <- rep(1, sum(low))
      for (k in seq_along(nn)) { acc <- acc + w[k] * p; p <- p * tl }
      v[low] <- acc
    }
    if (any(!low)) {
      hi <- th[!low]
      M <- log1p(hi) / hi
      f <- 1
      if (b >= 1L) for (m in 0:(b - 1L)) {
        M <- ((-1)^m * f / (1 + hi)^(m + 1L) - (m + 1L) * M) / hi
        f <- f * (m + 1L)
      }
      v[!low] <- M
    }
    out[[b + 1L]] <- v
  }
  out
}

#' @rdname negbin1_components
#' @description
#' `nb1_components_exact()` is the cancellation-free assembly, from the form in
#' \eqn{\mu + i\theta} that `negbin1_psums_cpp()` documents.
#' @keywords internal
nb1_components_exact <- function(y, mu, th, order) {
  S <- negbin1_psums_cpp(y, mu, th, as.integer(order))
  Md <- nb1_M_derivs(th, order)
  nms <- deriv_names(c("mu", "theta"), order)
  stats::setNames(lapply(nms, function(nm) {
    parts <- strsplit(nm, "_")[[1]]
    a <- sum(parts == "mu")
    b <- sum(parts == "theta")
    out <- (-1)^(order - 1L) * factorial(order - 1L) * S[, b + 1L]
    if (a == 0L) {
      out <- out - mu * Md[[b + 1L]] -
        y * (-1)^(b - 1L) * factorial(b - 1L) / (1 + th)^b
    } else if (a == 1L) {
      out <- out - Md[[b + 1L]]
    }
    out
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
