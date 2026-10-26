#' @include distrib.R generics.R
NULL

#' @title Beta-Binomial Distribution Class, Mean Proportion and Dispersion
#' @name BetaBinom1Distrib
#'
#' @description
#' This is the S7 class of the beta-binomial family parametrized by a mean
#' proportion \eqn{\mu \in (0, 1)} and a dispersion \eqn{\sigma > 0}, on the
#' finite support \eqn{\{0, 1, \dots, n\}}. It inherits from `discrete_distrib`.
#'
#' The class carries an extra property beyond the parent's, `size`: the number
#' of trials \eqn{n}, fixed at construction as it is for
#' [BinomialDistrib()]. Build one with [betabinom1_distrib()], which
#' validates `size`, supplies the two link functions and fills the properties
#' in. This page documents the raw S7 constructor, which takes the parent's
#' properties plus `size`; the validator of `distrib` checks the parent's
#' properties, and the constructor adds no check of its own.
#'
#' @inheritParams distrib
#' @param size The number of trials \eqn{n}, a single positive integer stored
#'   as a numeric. It belongs to the object, so an object cannot be reused
#'   across data sets whose group sizes differ.
#'
#' @return An S7 object of class `BetaBinom1Distrib`, inheriting from
#'   `discrete_distrib` and from `distrib`. Beyond `size` its properties are
#'   the parent's: `distrib_name`, `dimension`, `bounds`, `params`,
#'   `params_interpretation`, `n_params`, `params_bounds`, `link_params` and
#'   `params_smooth`. For an object built by [betabinom1_distrib()] they hold
#'   `"beta-binomial [size=n]"`, `"univariate"`, `c(0, size)`,
#'   `c("mu", "sigma")`, the interpretations
#'   `c(mu = "mean proportion", sigma = "dispersion")`, `2`, the domains
#'   \eqn{(0, 1)} and \eqn{(0, \infty)}, and the two links.
#'
#' @seealso [betabinom1_distrib()] to build one;
#'   [betabinom2_distrib()] for the same law in its two beta shapes;
#'   [binomial_distrib()] for the limit at \eqn{\sigma \to 0};
#'   [distrib_pdf.BetaBinom1Distrib()] for the mass function.
#'
#' @section Methods:
#' The following methods are registered on this class, and each is linked to
#' the page that documents it:
#' - the distribution: [`distrib_pdf()`][distrib_pdf.BetaBinom1Distrib],
#'   [`distrib_cdf()`][distrib_cdf.BetaBinom1Distrib],
#'   [`distrib_quantile()`][distrib_quantile.BetaBinom1Distrib] and
#'   [`distrib_rng()`][distrib_rng.BetaBinom1Distrib];
#' - the derivatives in the parameters:
#'   [`distrib_gradient()`][distrib_gradient.BetaBinom1Distrib],
#'   [`distrib_hessian()`][distrib_hessian.BetaBinom1Distrib],
#'   [`distrib_expected_hessian()`][distrib_expected_hessian.BetaBinom1Distrib],
#'   [`distrib_deriv3()`][distrib_deriv3.BetaBinom1Distrib] and
#'   [`distrib_deriv4()`][distrib_deriv4.BetaBinom1Distrib];
#' - the derivatives of the expected information:
#'   [`distrib_dexpected_hessian()`][distrib_dexpected_hessian.BetaBinom1Distrib]
#'   and
#'   [`distrib_d2expected_hessian()`][distrib_d2expected_hessian.BetaBinom1Distrib];
#' - the moments: [`mean()`][mean.BetaBinom1Distrib],
#'   [`variance()`][variance.BetaBinom1Distrib],
#'   [`skewness()`][skewness.BetaBinom1Distrib] and
#'   [`kurtosis()`][kurtosis.BetaBinom1Distrib];
#' - the fitting helpers: [distrib_scalar_route()].
#'
#' Every other generic is inherited from [discrete_distrib()].
#'
#' @examples
#' d <- betabinom1_distrib(size = 10)
#' S7::S7_inherits(d, discrete_distrib)
#'
#' # The trial count is a property of the object, not an entry of theta.
#' d@size
#' d@bounds
#'
#' # These properties let a consumer drive the family with no code specific to
#' # it.
#' d@params
#' d@params_interpretation
#'
#' # The mean proportion rides a logit and the dispersion a log.
#' vapply(d@link_params, function(l) l@link_name, character(1))
BetaBinom1Distrib <- S7::new_class("BetaBinom1Distrib",
  parent = discrete_distrib,
  properties = list(size = S7::class_numeric)
)

# --- S7 METHODS IMPLEMENTATION ---

#' @title Beta-Binomial Probability Mass Function
#' @name distrib_pdf.BetaBinom1Distrib
#' @description
#' Computes the beta-binomial mass
#' \deqn{P(Y = y) = \binom{n}{y}
#'       \dfrac{B(y + \alpha,\; n - y + \beta)}{B(\alpha, \beta)}, \qquad
#'       \alpha = \dfrac{\mu}{\sigma}, \quad \beta = \dfrac{1-\mu}{\sigma},}
#' with \eqn{B} the beta function and \eqn{n} the object's `size`. Off the
#' support \eqn{\{0, \dots, n\}} the mass is 0, and a non-integer `y` is off
#' the support too.
#'
#' The two beta functions are of magnitude \eqn{(\alpha+\beta)\log(\alpha+\beta)}
#' while their difference is of order one, so a mass computed as that
#' difference loses every digit at a small \eqn{\sigma}. The compiled kernel
#' uses a sum of logarithms instead, because the shifts are integers, and
#' stays exact to the binomial limit.
#'
#' @param distrib A `BetaBinom1Distrib` object, from [betabinom1_distrib()].
#' @param y A numeric vector of counts. A value outside
#'   \eqn{\{0, \dots, n\}} or not an integer gives a mass of 0, or `-Inf` with
#'   `log = TRUE`. A missing value gives `NA` (`NaN` for `NaN`).
#' @param theta A named list with components `mu` and `sigma`, each a numeric
#'   vector of length 1 or of the length of `y`. A component of length 1 is
#'   recycled. `mu` must lie in \eqn{(0, 1)} and `sigma` be strictly positive.
#' @param log Logical of length 1. When `TRUE` the log-mass is returned.
#'   Defaults to `FALSE`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of probabilities in \eqn{[0, 1]}, of length
#'   `max(length(y), length(mu), length(sigma))`, one value per observation.
#'
#' @seealso [distrib_cdf.BetaBinom1Distrib()] for the cumulative sum,
#'   [distrib_gradient.BetaBinom1Distrib()] for the derivatives of the
#'   log-mass, [distrib_pdf()] for the generic and [betabinom1_distrib()] for
#'   the family.
#'
#' @examples
#' d <- betabinom1_distrib(size = 10)
#' th <- list(mu = 0.3, sigma = 0.5)
#'
#' # The mass over the whole support sums to one.
#' m <- distrib_pdf(d, 0:10, th)
#' round(m, 4)
#' sum(m)
#'
#' # It is the beta-binomial mass at the two implied shapes.
#' a <- 0.3 / 0.5; b <- 0.7 / 0.5
#' all.equal(m, choose(10, 0:10) * beta(0:10 + a, 10 - 0:10 + b) / beta(a, b))
#'
#' # Off the support, and at a non-integer count, the mass is zero.
#' distrib_pdf(d, c(-1, 2.5, 11), th)
#'
#' # As the dispersion goes to zero the family becomes the binomial, and the
#' # mass stays accurate there where the two beta functions would cancel.
#' max(abs(distrib_pdf(d, 0:10, list(mu = 0.3, sigma = 1e-6)) -
#'         dbinom(0:10, 10, 0.3)))
S7::method(distrib_pdf, BetaBinom1Distrib) <- function(distrib, y, theta, log = FALSE, ...) {
  out <- betabinom_logpmf_cpp(y, theta[[1]], theta[[2]], distrib@size)
  if (log) out else exp(out)
}

#' @title Beta-Binomial Cumulative Distribution Function
#' @name distrib_cdf.BetaBinom1Distrib
#' @description
#' Computes \eqn{F(q) = P(Y \le q) = \sum_{k \le \lfloor q \rfloor} P(Y = k)}
#' as the cumulative sum of the mass over the whole support
#' \eqn{\{0, \dots, n\}}. Because the support is finite, that sum is exact
#' rather than an approximation, and it costs \eqn{n+1} evaluations of the mass for
#' each distinct parameter setting. Below the support the probability is 0
#' and at or above \eqn{n} it is 1. With `lower.tail = FALSE` the upper tail
#' \eqn{P(Y > q)} is the sum of the mass above \eqn{q}, accumulated from the
#' top of the support, so an upper-tail probability far below machine epsilon
#' is not formed as \eqn{1 - F} and keeps its accuracy, on the log scale
#' too. The helper is [betabinom1_cum()].
#'
#' @param distrib A `BetaBinom1Distrib` object, from [betabinom1_distrib()].
#' @param q A numeric vector of quantiles. A non-integer value is floored, so
#'   `F(2.7)` is `F(2)`.
#' @param theta A named list with components `mu` and `sigma`, each a numeric
#'   vector of length 1 or of the length of `q`. `mu` must lie in
#'   \eqn{(0, 1)} and `sigma` be strictly positive. With scalar parameters one
#'   cumulative table serves the whole call; parameters that vary by
#'   observation give one table per observation.
#' @param lower.tail Logical of length 1. When `TRUE`, the default,
#'   probabilities are \eqn{P(Y \le q)}; when `FALSE` they are \eqn{P(Y > q)}.
#' @param log.p Logical of length 1. When `TRUE` the logarithm of the
#'   probability is returned. Defaults to `FALSE`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of probabilities in \eqn{[0, 1]}, of length
#'   `max(length(q), length(mu), length(sigma))`. With `log.p = TRUE` the
#'   values are logarithms and are non-positive.
#'
#' @seealso [distrib_quantile.BetaBinom1Distrib()] for the generalized
#'   inverse, [distrib_pdf.BetaBinom1Distrib()] for the mass, and
#'   [distrib_cdf()] for the generic.
#'
#' @examples
#' d <- betabinom1_distrib(size = 10)
#' th <- list(mu = 0.3, sigma = 0.5)
#'
#' # The distribution function is the cumulative sum of the mass, exactly.
#' all.equal(distrib_cdf(d, 0:10, th), cumsum(distrib_pdf(d, 0:10, th)))
#'
#' # Below the support it is zero and at the top of it one.
#' distrib_cdf(d, c(-1, 0, 3, 10, 11), th)
#'
#' # A non-integer quantile is floored.
#' c(distrib_cdf(d, 2.7, th), distrib_cdf(d, 2, th))
#'
#' # The two tails sum to one.
#' distrib_cdf(d, 4, th) + distrib_cdf(d, 4, th, lower.tail = FALSE)
S7::method(distrib_cdf, BetaBinom1Distrib) <- function(distrib, q, theta,
                                                      lower.tail = TRUE,
                                                      log.p = FALSE, ...) {
  n <- distrib@size
  k <- floor(q)
  cum <- betabinom1_cum(n, theta, length(q), upper = !lower.tail)
  if (is.matrix(cum)) k <- rep_len(k, nrow(cum))
  kk <- pmin(pmax(k, 0), n) + 1L
  at <- if (is.matrix(cum)) cum[cbind(seq_along(k), kk)] else cum[kk]
  below <- if (lower.tail) 0 else 1
  p <- ifelse(k < 0, below, ifelse(k >= n, 1 - below, at))
  p <- pmin(pmax(p, 0), 1)
  if (log.p) log(p) else p
}

#' @title Beta-Binomial Quantile Function
#' @name distrib_quantile.BetaBinom1Distrib
#' @description
#' Computes the generalized inverse \eqn{Q(p) = \min\{y : F(y) \ge p\}} by
#' walking the exact cumulative sum over the support. The result is an integer
#' count because the distribution function is a step function, so `Q(F(y))`
#' returns `y` while `F(Q(p))` is at least `p` and generally exceeds it. A
#' relative tolerance of 64 machine epsilons is allowed on the comparison so
#' that a `p` obtained from [distrib_cdf.BetaBinom1Distrib()] maps back to the
#' count from which it came. With `lower.tail = FALSE` the upper tail is
#' inverted against \eqn{P(Y > k)} summed from the top of the support, so a
#' small upper-tail probability gives the count to which it belongs.
#'
#' @param distrib A `BetaBinom1Distrib` object, from [betabinom1_distrib()].
#' @param p A numeric vector of probabilities in \eqn{[0, 1]}, or of their
#'   logarithms when `log.p = TRUE`. `NA` is returned for `NA` and `NaN` for a
#'   value outside \eqn{[0, 1]}; \eqn{p = 0} gives 0 and \eqn{p = 1} gives
#'   \eqn{n} (the reverse with `lower.tail = FALSE`).
#' @param theta A named list with components `mu` and `sigma`, each a numeric
#'   vector of length 1 or of the length of `p`. `mu` must lie in
#'   \eqn{(0, 1)} and `sigma` be strictly positive. With scalar parameters one
#'   cumulative table serves the whole call; parameters that vary by
#'   observation give one table per observation.
#' @param lower.tail Logical of length 1. When `TRUE`, the default, `p` is
#'   \eqn{P(Y \le q)}; when `FALSE` it is \eqn{P(Y > q)}.
#' @param log.p Logical of length 1. When `TRUE`, `p` is read as a logarithm.
#'   Defaults to `FALSE`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of counts in \eqn{\{0, \dots, n\}}, of length
#'   `max(length(p), length(mu), length(sigma))`.
#'
#' @seealso [distrib_cdf.BetaBinom1Distrib()] for the function inverted here,
#'   [distrib_rng.BetaBinom1Distrib()] for draws, and [distrib_quantile()] for
#'   the generic.
#'
#' @examples
#' d <- betabinom1_distrib(size = 10)
#' th <- list(mu = 0.3, sigma = 0.5)
#'
#' # The deciles are counts.
#' distrib_quantile(d, c(0.1, 0.5, 0.9), th)
#'
#' # The round trip from a count through F and back is the identity.
#' all.equal(distrib_quantile(d, distrib_cdf(d, 0:10, th), th), as.numeric(0:10))
#'
#' # The other direction only reaches at least p, because F is a step function.
#' rbind(p = c(0.2, 0.45, 0.8),
#'       F_at_Q = distrib_cdf(d, distrib_quantile(d, c(0.2, 0.45, 0.8), th), th))
S7::method(distrib_quantile, BetaBinom1Distrib) <- function(distrib, p, theta,
                                                            lower.tail = TRUE,
                                                            log.p = FALSE, ...) {
  n <- distrib@size
  fuzz <- 64 * .Machine$double.eps
  # the lower tail is inverted against F(k), the upper tail against
  # P(Y > k) summed from the top, so that neither probability is formed as
  # one minus the other
  cum <- betabinom1_cum(n, theta, length(p), upper = !lower.tail)
  if (is.matrix(cum)) p <- rep_len(p, nrow(cum))
  vapply(seq_along(p), function(i) {
    pp <- p[i]
    if (is.na(pp)) return(pp)
    pr <- if (log.p) exp(pp) else pp
    if (pr < 0 || pr > 1) return(NaN)
    ci <- if (is.matrix(cum)) cum[i, ] else cum
    if (lower.tail) {
      if (pr == 0) return(0)
      hit <- if (log.p) log(ci) >= pp - fuzz else ci >= pr * (1 - fuzz)
    } else {
      if (pr == 1) return(0)
      hit <- if (log.p) log(ci) <= pp + fuzz else ci <= pr * (1 + fuzz)
    }
    (0:n)[which(hit)[1L]]
  }, numeric(1))
}

#' The Beta-Binomial's Cumulative Mass Over Its Support
#'
#' @description
#' Returns the cumulative mass over \eqn{\{0, \dots, n\}}: a vector where both
#' parameters are scalar, and otherwise a matrix with one row per observation,
#' in which each observation's parameters are crossed with the support. With
#' `upper = TRUE` the entry at \eqn{k} is \eqn{P(Y > k)}, the mass above
#' \eqn{k} accumulated from the top of the support, in place of
#' \eqn{P(Y \le k)}.
#'
#' @details
#' The kernel reads its parameters element by element, so the support and a
#' vector of parameters are not passed to it together: each observation's
#' parameters are crossed with the support, and the sum is taken along each
#' row.
#'
#' @param n The size.
#' @param theta A named list with the two parameters.
#' @param len The length of the quantity at which the caller evaluates the
#'   cumulative mass.
#' @param upper Logical of length 1. When `TRUE` the upper tail is returned.
#'   Defaults to `FALSE`.
#'
#' @return A numeric vector of length `n + 1`, or a matrix with `n + 1`
#'   columns.
#'
#' @keywords internal
betabinom1_cum <- function(n, theta, len, upper = FALSE) {
  a <- theta[[1]]
  b <- theta[[2]]
  if (length(a) == 1L && length(b) == 1L) {
    f <- exp(betabinom_logpmf_cpp(0:n, a, b, n))
    if (!upper) return(cumsum(f))
    return(c(rev(cumsum(rev(f)))[-1L], 0))
  }
  nn <- max(len, length(a), length(b))
  m <- n + 1L
  mass <- exp(betabinom_logpmf_cpp(rep(0:n, each = nn),
                                   rep(rep_len(a, nn), times = m),
                                   rep(rep_len(b, nn), times = m), n))
  mass <- matrix(mass, nn, m)
  if (!upper) {
    if (m > 1L) for (j in 2:m) mass[, j] <- mass[, j] + mass[, j - 1L]
    return(mass)
  }
  # P(Y > k): the mass above k, accumulated from the top of the support
  out <- matrix(0, nn, m)
  if (m > 1L) for (j in (m - 1L):1L) out[, j] <- out[, j + 1L] + mass[, j + 1L]
  out
}

#' @title Beta-Binomial Random Generation
#' @name distrib_rng.BetaBinom1Distrib
#' @description
#' Draws `n` independent beta-binomial counts by the two-stage hierarchy that
#' defines the family: a success probability from
#' \eqn{\mathrm{Beta}(\alpha, \beta)} with \eqn{\alpha = \mu/\sigma} and
#' \eqn{\beta = (1-\mu)/\sigma}, then a count from
#' \eqn{\mathrm{Binomial}(n_{\mathrm{trials}}, p)} at that probability, one
#' fresh probability per draw. The draws are made by one call to
#' [stats::rbeta()] followed by one call to [stats::rbinom()], so they depend
#' on `.Random.seed` in the usual way.
#'
#' @param distrib A `BetaBinom1Distrib` object, from [betabinom1_distrib()].
#' @param n A single positive integer, the number of draws. The number of
#'   trials is the object's `size` property, not this argument.
#' @param theta A named list with components `mu` and `sigma`, each a numeric
#'   vector of length 1 or of length `n`. A component of length 1 is recycled.
#'   `mu` must lie in \eqn{(0, 1)} and `sigma` be strictly positive.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of `n` counts in \eqn{\{0, \dots, size\}}.
#'
#' @seealso [distrib_quantile.BetaBinom1Distrib()] for the inverse-transform
#'   route, [fit_distrib()] to estimate the parameters back from a sample, and
#'   [distrib_rng()] for the generic.
#'
#' @examples
#' d <- betabinom1_distrib(size = 10)
#' th <- list(mu = 0.3, sigma = 0.5)
#'
#' # The sample moments recover the mean and the overdispersed variance.
#' set.seed(1)
#' z <- distrib_rng(d, 2e5, th)
#' rbind(sample = c(mean = mean(z), var = var(z)),
#'       theoretical = c(mean(d, th), variance(d, th)))
#'
#' # The counts are bounded by the trial count, which is the object's size and
#' # not the argument n.
#' range(z)
S7::method(distrib_rng, BetaBinom1Distrib) <- function(distrib, n, theta, ...) {
  a <- theta[[1]] / theta[[2]]
  b <- (1 - theta[[1]]) / theta[[2]]
  stats::rbinom(n, size = distrib@size, prob = stats::rbeta(n, a, b))
}

#' @title Beta-Binomial Score
#' @name distrib_gradient.BetaBinom1Distrib
#' @description
#' Computes the first derivatives of the beta-binomial log-mass with respect to
#' the mean proportion \eqn{\mu} and the dispersion \eqn{\sigma}, one value per
#' observation, in closed form. The parameters enter only through the two beta
#' shapes, where each derivative is a difference of digammas:
#' \deqn{\dfrac{\partial \ell}{\partial \alpha}
#'       = \psi(y+\alpha) - \psi(\alpha) - \psi(n+S) + \psi(S), \qquad
#'       S = \alpha + \beta,}
#' and likewise in \eqn{\beta} with \eqn{n - y} in place of \eqn{y}. The
#' reported components follow by the chain rule of
#' \eqn{(\alpha, \beta) = (\mu/\sigma, (1-\mu)/\sigma)}.
#'
#' The kernel uses these shape forms for \eqn{\sigma \ge 0.1}, where they are
#' accurate. As \eqn{\sigma} goes to zero (the binomial limit) the derivatives
#' of the map grow like powers of \eqn{1/\sigma} and the terms of the chain
#' rule cancel, so for \eqn{\sigma < 0.1} the kernel uses the direct sums of
#' [betabinom1_distrib()] instead: with \eqn{u_i = \mu + i\sigma},
#' \eqn{v_i = 1 - \mu + i\sigma} and \eqn{w_i = 1 + i\sigma},
#' \deqn{\dfrac{\partial \ell}{\partial \mu}
#'       = \sum_{i<y} \dfrac{1}{u_i} - \sum_{i<n-y} \dfrac{1}{v_i}, \qquad
#'       \dfrac{\partial \ell}{\partial \sigma}
#'       = \sum_{i<y} \dfrac{i}{u_i} + \sum_{i<n-y} \dfrac{i}{v_i}
#'         - \sum_{i<n} \dfrac{i}{w_i}.}
#' Each sum has terms of one sign, so both components keep their accuracy at
#' any concentration; the cost is \eqn{O(n)} per observation.
#'
#' With `scale = "link"` the generic applies the chain rule for the links that
#' the family carries before returning. This method always returns the parameter
#' scale. The arithmetic runs in a compiled kernel decomposed over the elements
#' of the output, so the result does not depend on the thread count.
#'
#' @param distrib A `BetaBinom1Distrib` object, from [betabinom1_distrib()].
#' @param y A numeric vector of counts in \eqn{\{0, \dots, n\}}.
#' @param theta A named list with components `mu` and `sigma`, each a numeric
#'   vector of length 1 or of the length of `y`. A component of length 1 is
#'   recycled. `mu` must lie in \eqn{(0, 1)} and `sigma` be strictly positive.
#' @param scale One of `"parameter"` (the default) or `"link"`, matched by
#'   [base::match.arg()]. The generic reads it, and this method does not.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#' @param threads A single positive integer, how many threads the kernel may
#'   use. For fewer observations than an internal threshold the kernel runs
#'   sequentially whatever the count. Defaults to `1L`.
#'
#' @return A named list of two numeric vectors, `mu` and `sigma`, each of
#'   length `max(length(y), length(mu), length(sigma))`.
#'
#' @section Notation:
#' \eqn{\ell} is the log-mass of one observation, \eqn{\mu \in (0,1)} the mean
#' proportion, \eqn{\sigma > 0} the dispersion, \eqn{n} the trial count and
#' \eqn{\psi} the digamma function, `digamma()` in R.
#'
#' @seealso [distrib_hessian.BetaBinom1Distrib()] for the second derivatives,
#'   [distrib_expected_hessian.BetaBinom1Distrib()] for their expectation,
#'   [distrib_gradient.BetaBinom2Distrib()] for the same quantity in the
#'   shapes, and [distrib_gradient()] for the generic.
#'
#' @examples
#' d <- betabinom1_distrib(size = 10)
#' th <- list(mu = 0.3, sigma = 0.5)
#' g <- distrib_gradient(d, 0:10, th)
#'
#' # It is the derivative of the log-mass, so a central difference of it
#' # reproduces the mu component.
#' eps <- 1e-6
#' all.equal((distrib_pdf(d, 0:10, list(mu = 0.3 + eps, sigma = 0.5), log = TRUE) -
#'            distrib_pdf(d, 0:10, list(mu = 0.3 - eps, sigma = 0.5), log = TRUE)) /
#'             (2 * eps), g$mu, tolerance = 1e-6)
#'
#' # The score has mean zero over the support: the first Bartlett identity.
#' w <- distrib_pdf(d, 0:10, th)
#' vapply(g, function(v) sum(w * v), numeric(1))
#'
#' # At a tiny dispersion the family is a binomial, and the mu score is the
#' # binomial one, y/mu - (n-y)/(1-mu), evaluated here at y = 3.
#' c(kernel = distrib_gradient(d, 3, list(mu = 0.4, sigma = 1e-12))$mu,
#'   binomial = 3 / 0.4 - 7 / 0.6)
S7::method(distrib_gradient, BetaBinom1Distrib) <- function(distrib, y, theta,
                                                            scale = c("parameter", "link"), ...,
                                       threads = 1L) {
  betabinom_gradient_cpp(y, theta[[1]], theta[[2]], distrib@size, threads)
}

#' @title Beta-Binomial Observed Hessian
#' @name distrib_hessian.BetaBinom1Distrib
#' @description
#' Computes the three distinct second derivatives of the beta-binomial log-mass
#' with respect to \eqn{\mu} and \eqn{\sigma}, one value per observation, in
#' closed form. In the shapes each second derivative is a difference of
#' trigammas,
#' \deqn{\dfrac{\partial^2 \ell}{\partial \alpha^2}
#'       = \psi_1(y+\alpha) - \psi_1(\alpha) - \psi_1(n+S) + \psi_1(S), \qquad
#'       S = \alpha + \beta,}
#' and the mixed shape component carries only the \eqn{S} part,
#' \eqn{-\psi_1(n+S) + \psi_1(S)}, because the two shapes otherwise enter the
#' mass separately. The reported components follow by the two-variable chain rule of
#' \eqn{(\alpha, \beta) = (\mu/\sigma, (1-\mu)/\sigma)}, whose own second
#' derivatives contribute the rest. This route is taken for
#' \eqn{\sigma \ge 0.1}; for \eqn{\sigma < 0.1} the components are the direct
#' sums of order two given on [betabinom1_distrib()], for the reason given on
#' [distrib_gradient.BetaBinom1Distrib()].
#'
#' The arithmetic runs in a compiled kernel decomposed over the elements of the
#' output, so the result does not depend on the thread count.
#'
#' @param distrib A `BetaBinom1Distrib` object, from [betabinom1_distrib()].
#' @param y A numeric vector of counts in \eqn{\{0, \dots, n\}}.
#' @param theta A named list with components `mu` and `sigma`, each a numeric
#'   vector of length 1 or of the length of `y`. A component of length 1 is
#'   recycled. `mu` must lie in \eqn{(0, 1)} and `sigma` be strictly positive.
#' @param scale One of `"parameter"` (the default) or `"link"`, matched by
#'   [base::match.arg()]. The generic reads it, and this method does not.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#' @param threads A single positive integer, how many threads the kernel may
#'   use. For fewer observations than an internal threshold the kernel runs
#'   sequentially whatever the count. Defaults to `1L`.
#'
#' @return A named list of three numeric vectors, `mu_mu`, `mu_sigma` and
#'   `sigma_sigma`, each of length
#'   `max(length(y), length(mu), length(sigma))`. The three name the distinct
#'   entries of a symmetric \eqn{2 \times 2} matrix per observation.
#'
#' @section Notation:
#' \eqn{\ell} is the log-mass of one observation, \eqn{\mu \in (0,1)} the mean
#' proportion, \eqn{\sigma > 0} the dispersion, \eqn{n} the trial count and
#' \eqn{\psi_1} the trigamma function, `trigamma()` in R.
#'
#' @seealso [distrib_gradient.BetaBinom1Distrib()] for the score,
#'   [distrib_expected_hessian.BetaBinom1Distrib()] for the expectation of this
#'   quantity, [distrib_deriv3.BetaBinom1Distrib()] for the order above, and
#'   [distrib_hessian()] for the generic.
#'
#' @examples
#' d <- betabinom1_distrib(size = 10)
#' th <- list(mu = 0.3, sigma = 0.5)
#' h <- distrib_hessian(d, 0:10, th)
#' names(h)
#'
#' # A central difference of the score reproduces the pure-mu component.
#' eps <- 1e-5
#' up <- distrib_gradient(d, 0:10, list(mu = 0.3 + eps, sigma = 0.5))$mu
#' dn <- distrib_gradient(d, 0:10, list(mu = 0.3 - eps, sigma = 0.5))$mu
#' all.equal((up - dn) / (2 * eps), h$mu_mu, tolerance = 1e-6)
#'
#' # The mass-weighted sum over the support is the expected Hessian.
#' w <- distrib_pdf(d, 0:10, th)
#' rbind(summed = vapply(h, function(v) sum(w * v), numeric(1)),
#'       expected = vapply(distrib_expected_hessian(d, 0:10, th),
#'                         function(v) v[1], numeric(1)))
S7::method(distrib_hessian, BetaBinom1Distrib) <- function(distrib, y, theta,
                                                           scale = c("parameter", "link"), ...,
                                       threads = 1L) {
  betabinom_hessian_cpp(y, theta[[1]], theta[[2]], distrib@size, threads)
}

#' @title Beta-Binomial Expected Hessian
#' @name distrib_expected_hessian.BetaBinom1Distrib
#' @description
#' Returns the expectation of the observed Hessian under the model,
#' \eqn{\sum_{k=0}^{n} P(Y=k)\, \partial^2\ell/\partial\theta_i\partial\theta_j}
#' evaluated at \eqn{y = k}. Because the support is finite, the expectation is
#' an exact finite sum of \eqn{n+1} terms rather than a quadrature or a
#' sample. The summand is the observed Hessian of
#' [distrib_hessian.BetaBinom1Distrib()], which comes from the direct sums for
#' \eqn{\sigma < 0.1} and from the shapes above.
#'
#' `approx` is therefore ignored; every approximation returns the same
#' three numbers. The compiled kernel is sequential: it evaluates the sum once
#' for each run of consecutive observations that share their parameters, and
#' `threads` is ignored.
#'
#' @param distrib A `BetaBinom1Distrib` object, from [betabinom1_distrib()].
#' @param y A numeric vector of counts. Only its length is read, because the
#'   expectation does not depend on the data; the values are ignored.
#' @param theta A named list with components `mu` and `sigma`, each a numeric
#'   vector of length 1 or of the length of `y`. A component of length 1 is
#'   recycled. `mu` must lie in \eqn{(0, 1)} and `sigma` be strictly positive.
#' @param scale One of `"parameter"` (the default) or `"link"`, matched by
#'   [base::match.arg()]. The generic reads it, and this method does not.
#' @param approx An object inheriting from [ExpectApprox]. Ignored here, because
#'   the expectation is an exact sum. It is accepted so that the signature
#'   matches the generic's, where it selects the approximation of the
#'   expectation. See [expected_derivative_methods].
#' @param ... Unused, and accepted so that the signature matches the generic's.
#' @param threads A single positive integer. It is accepted so that the
#'   signature matches the other derivative methods; this kernel is
#'   sequential and ignores it. Defaults to `1L`.
#'
#' @return A named list of three numeric vectors, `mu_mu`, `mu_sigma` and
#'   `sigma_sigma`, each of length `length(y)` and each constant along it.
#'
#' @seealso [distrib_hessian.BetaBinom1Distrib()] for the quantity this is the
#'   expectation of, [distrib_gradient.BetaBinom1Distrib()] for the score, and
#'   [distrib_expected_hessian()] for the generic.
#'
#' @examples
#' d <- betabinom1_distrib(size = 10)
#' th <- list(mu = 0.3, sigma = 0.5)
#' eh <- distrib_expected_hessian(d, 0:10, th)
#' vapply(eh, function(v) v[1], numeric(1))
#'
#' # It is the mass-weighted sum of the observed Hessian over the support,
#' # written out here by hand and agreeing exactly.
#' w <- distrib_pdf(d, 0:10, th)
#' vapply(distrib_hessian(d, 0:10, th), function(v) sum(w * v), numeric(1))
#'
#' # The mean proportion and the dispersion are not orthogonal: the mixed
#' # entry does not vanish, so their estimates are asymptotically correlated.
#' eh$mu_sigma[1]
#'
#' # The approx argument has no effect, because the expectation is an exact
#' # sum.
#' identical(eh, distrib_expected_hessian(d, 0:10, th, approx = expect_mc(nsim = 50)))
S7::method(distrib_expected_hessian, BetaBinom1Distrib) <- function(distrib, y, theta,
                                                                    scale = c("parameter", "link"),
                                                                    approx = expect_opg(), ...,
                                       threads = 1L) {
  betabinom_expected_hessian_cpp(y, theta[[1]], theta[[2]], distrib@size, threads)
}

# --- CONSTRUCTOR WRAPPER ---

#' Beta-Binomial Distribution, Mean Proportion and Dispersion
#'
#' @description
#' Builds the distribution object for the beta-binomial family parametrized by
#' a mean proportion \eqn{\mu \in (0, 1)} and a dispersion \eqn{\sigma > 0}, on
#' the finite support \eqn{\{0, 1, \dots, n\}}. The returned object carries
#' closed-form derivatives of the log-mass to fourth order and expectations
#' that are exact finite sums over the support, so the object does not use
#' quadrature. The fifth derivative comes from the numerical fallback of the base class.
#'
#' The family is the binomial with its success probability drawn from a beta,
#' which makes it the standard model for a proportion whose trials are not
#' independent. For \eqn{n > 1} it is overdispersed relative to a binomial
#' with the same mean at every \eqn{\sigma > 0}, and it approaches the
#' binomial as \eqn{\sigma \to 0}. For \eqn{n = 1} it is the Bernoulli
#' distribution at every \eqn{\sigma}.
#'
#' @param size The number of trials \eqn{n}, a single positive integer. It is a
#'   constant of the distribution and not a parameter, as for
#'   [binomial_distrib()], so an object cannot be reused across data sets whose
#'   group sizes differ. Anything else signals an error naming the argument.
#' @param link_mu A `link` object from `linkfunctions7` for the mean proportion
#'   \eqn{\mu}. Defaults to [linkfunctions7::logit_link()], which maps
#'   \eqn{(0, 1)} onto the line.
#' @param link_sigma A `link` object from `linkfunctions7` for the dispersion
#'   \eqn{\sigma}. Defaults to [linkfunctions7::log_link()], which maps
#'   \eqn{(0, \infty)} onto the line and so keeps every fitted value positive.
#'
#' @details
#' # The parametrization
#'
#' Write \eqn{\alpha = \mu/\sigma} and \eqn{\beta = (1-\mu)/\sigma} for the two
#' beta shapes. The mass on \eqn{y \in \{0, \dots, n\}} is
#' \deqn{P(Y = y) = \binom{n}{y}
#'       \dfrac{B(y+\alpha,\; n-y+\beta)}{B(\alpha, \beta)},}
#' with \eqn{B} the beta function. [betabinom2_distrib()] is the same law
#' written in \eqn{(\alpha, \beta)} directly; the two agree mass for mass at
#' \eqn{\mu = \alpha/(\alpha+\beta)} and \eqn{\sigma = 1/(\alpha+\beta)}.
#'
#' The family is a binomial whose success probability is drawn from a beta
#' distribution. It is not a binomial with a parameter held fixed and not one
#' of the wrappers over [binomial_distrib()], because the mixing acts on the
#' success probability and not on the outcome.
#'
#' # Overdispersion
#'
#' The mean is \eqn{n\mu} and the variance
#' \deqn{\operatorname{Var}(Y) = n\mu(1-\mu)
#'       \left(1 + (n-1)\dfrac{\sigma}{1+\sigma}\right),}
#' so the inflation over a binomial of the same mean is
#' \eqn{1 + (n-1)\sigma/(1+\sigma)}, above 1 for \eqn{n > 1} and rising with
#' both the trial count and the dispersion. The factor \eqn{\sigma/(1+\sigma)} is the
#' intraclass correlation between two trials of the same group. At
#' \eqn{n = 10} and \eqn{\sigma = 0.5} that correlation is \eqn{1/3} and the
#' variance is four times a binomial's (8.4 against 2.1 at \eqn{\mu = 0.3}).
#'
#' As \eqn{\sigma \to 0} the mass converges to the binomial's at rate
#' \eqn{O(\sigma)}, and the compiled kernel stays accurate all the way there;
#' see the cancellation note below.
#'
#' # Derivatives
#'
#' The parameters enter only through the two shapes, where every derivative of
#' order \eqn{k} is a difference of \eqn{\psi^{(k-1)}}:
#' \deqn{\dfrac{\partial \ell}{\partial \alpha}
#'       = \psi(y+\alpha) - \psi(\alpha) - \psi(n+S) + \psi(S), \qquad
#'       S = \alpha + \beta,}
#' and likewise in \eqn{\beta}. The reported components follow by the chain
#' rule of \eqn{(\alpha, \beta) = (\mu/\sigma, (1-\mu)/\sigma)}, which is
#' linear in \eqn{\mu} at fixed \eqn{\sigma}, so every partial of the map
#' carrying two or more \eqn{\mu} is exactly zero and the third and fourth
#' orders are a short partition sum. Orders three and four are therefore
#' closed form as well, in [distrib_deriv3.BetaBinom1Distrib()], and not
#' the numerical fallback.
#'
#' The shifts of the beta functions are integers, so with
#' \eqn{u_i = \mu + i\sigma}, \eqn{v_i = 1 - \mu + i\sigma} and
#' \eqn{w_i = 1 + i\sigma} the log-mass is also
#' \deqn{\ell = \log\binom{n}{y} + \sum_{i<y}\log u_i + \sum_{i<n-y}\log v_i
#'   - \sum_{i<n}\log w_i,}
#' and a derivative with \eqn{a} differentiations in \eqn{\mu} and \eqn{b} in
#' \eqn{\sigma}, of order \eqn{K = a + b}, is
#' \deqn{c_K\Big[\sum_{i<y}\frac{i^b}{u_i^K} + (-1)^a\sum_{i<n-y}\frac{i^b}{v_i^K}
#'   - [a = 0]\sum_{i<n}\frac{i^b}{w_i^K}\Big], \qquad c_K = (-1)^{K-1}(K-1)!.}
#' Each sum has terms of one sign. These direct sums are used for
#' \eqn{\sigma < 0.1} at every order, observed and expected, and the route
#' through the shapes above it (see the next section).
#'
#' Every expectation, including those of the third and fourth derivatives, is
#' an exact finite sum over \eqn{\{0, \dots, n\}}, because the support is
#' bounded.
#'
#' # The cancellation at a small dispersion
#'
#' The two beta functions of the mass are of magnitude
#' \eqn{(\alpha+\beta)\log(\alpha+\beta)} while their difference is of order
#' one, so a mass computed as that difference loses one digit per factor of
#' ten in the concentration \eqn{1/\sigma}. The kernel computes the log-mass
#' as the sum of logarithms above instead. The route through the shapes has
#' the same defect in its derivatives: the polygamma differences are of order
#' one while the derivatives of the map \eqn{(\alpha, \beta)} in \eqn{\sigma}
#' are of order \eqn{\sigma^{-(K+1)}}, and their products cancel among
#' themselves as \eqn{\sigma} goes to zero. Below \eqn{\sigma = 0.1} the
#' derivatives and the expectations built from them therefore come from the
#' direct sums, at a cost of \eqn{O(n)} per observation; above it the shape
#' route is accurate, and the direct sums would lose about
#' \eqn{\log_{10} n} digits there.
#'
#' # Estimation
#'
#' [fit_distrib()] maximizes the log-likelihood on the link scale. Both
#' estimates lack a closed form. For a sample with no extra-binomial
#' variation, \eqn{\hat\sigma} tends to zero, the boundary of the parameter
#' space.
#'
#' @section Notation:
#' \eqn{\ell} is the log-mass of one observation, \eqn{\mu \in (0,1)} the mean
#' proportion, \eqn{\sigma > 0} the dispersion, \eqn{n} the trial count,
#' \eqn{B} the beta function and \eqn{\psi} the digamma function. \eqn{\eta} is
#' a parameter on the unconstrained scale of its link, with
#' \eqn{\theta = g^{-1}(\eta)}.
#'
#' @return An S7 object of class `BetaBinom1Distrib`, inheriting from
#'   `discrete_distrib`, with `size` the trial count, `distrib_name`
#'   `"beta-binomial [size=n]"`, `dimension` `"univariate"`, `bounds`
#'   `c(0, size)`, `params` `c("mu", "sigma")`, `n_params` `2`,
#'   `params_bounds` the domains \eqn{(0, 1)} and \eqn{(0, \infty)}, and
#'   `link_params` the two links given here.
#'
#' @references
#' Skellam, J. G. (1948). A probability distribution derived from the binomial
#' distribution by regarding the probability of success as variable between the
#' sets of trials. *Journal of the Royal Statistical Society, Series B*,
#' **10**(2), 257-261.
#'
#' Johnson, N. L., Kemp, A. W. and Kotz, S. (2005).
#' *Univariate Discrete Distributions*, 3rd edition, Section 6.9.
#' Wiley, Hoboken.
#'
#' @importFrom linkfunctions7 logit_link log_link
#' @importFrom stats rbeta rbinom
#'
#' @examples
#' d <- betabinom1_distrib(size = 10)
#' d
#'
#' # The mass over the support sums to one.
#' th <- list(mu = 0.3, sigma = 0.5)
#' sum(distrib_pdf(d, 0:10, th))
#'
#' # The variance is four times that of a binomial with the same mean: the
#' # intraclass correlation is sigma / (1 + sigma) = 1/3, and 1 + 9/3 = 4.
#' c(mean = mean(d, th), var = variance(d, th),
#'   binomial_var = variance(binomial_distrib(size = 10), list(mu = 0.3)),
#'   icc = 0.5 / 1.5)
#'
#' # The dispersion going to zero is the binomial, at rate O(sigma).
#' vapply(c(1e-2, 1e-4, 1e-6), function(s)
#'   max(abs(distrib_pdf(d, 0:10, list(mu = 0.3, sigma = s)) -
#'           dbinom(0:10, 10, 0.3))), numeric(1))
#'
#' # The mass equals that of betabinom2 at the implied shapes.
#' all.equal(distrib_pdf(d, 0:10, list(mu = 0.4, sigma = 0.2)),
#'           distrib_pdf(betabinom2_distrib(size = 10), 0:10,
#'                       list(alpha = 2, beta = 3)))
#'
#' # Fitting recovers both parameters.
#' set.seed(3)
#' z <- distrib_rng(d, 2000, list(mu = 0.35, sigma = 0.4))
#' coef(fit_distrib(d, z))
#'
#' @seealso
#' [betabinom2_distrib()] for the same law in its two beta shapes;
#' [binomial_distrib()] for the limit at \eqn{\sigma \to 0} and
#' [beta1_distrib()] for the mixing law; [negbin2_distrib()] for the
#' unbounded-count analogue, a Poisson mixed over a gamma; [fit_distrib()] to
#' estimate the parameters; [check_distrib()] to validate a user-defined family
#' against the same battery that this one passes; [BetaBinom1Distrib] for the
#' class.
#' @export
betabinom1_distrib <- function(size, link_mu = logit_link(),
                              link_sigma = log_link()) {
  if (!is.numeric(size) || length(size) != 1L || !is.finite(size) ||
      size < 1 || size != round(size)) {
    stop("'size' must be a single positive integer.", call. = FALSE)
  }
  BetaBinom1Distrib(
    size = as.numeric(size),
    distrib_name = sprintf("beta-binomial [size=%g]", size),
    dimension = "univariate", bounds = c(0, size),
    params = c("mu", "sigma"),
    params_interpretation = c(mu = "mean proportion", sigma = "dispersion"),
    n_params = 2, params_bounds = list(mu = c(0, 1), sigma = c(0, Inf)),
    link_params = list(mu = link_mu, sigma = link_sigma)
  )
}
