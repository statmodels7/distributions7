#' @include distrib.R generics.R
NULL

#' @title Beta-Binomial Distribution Class, Mean Proportion and Dispersion
#' @name BetaBinom1Distrib
#'
#' @description
#' The S7 class of the beta-binomial family parametrized by a mean proportion
#' \eqn{\mu \in (0, 1)} and a dispersion \eqn{\sigma > 0}, on the finite
#' support \eqn{\{0, 1, \dots, n\}}. It inherits from `discrete_distrib`, so it
#' answers every generic of the `distrib` contract; the seven methods listed
#' below are registered on it in this file, two more in `betabinom1_higher.R`,
#' and everything else comes from the parent.
#'
#' The class carries an extra property beyond the parent's, `size`: the number
#' of trials \eqn{n}, fixed at construction as it is for
#' [BinomialDistrib()]. Build one with [betabinom1_distrib()], which
#' validates `size`, supplies the two link functions and fills the properties
#' in. This page documents the raw S7 constructor, which validates none of the
#' relationships between them.
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
#' Registered on this class in this file:
#'   [`distrib_cdf()`][distrib_cdf.BetaBinom1Distrib],
#'   [`distrib_expected_hessian()`][distrib_expected_hessian.BetaBinom1Distrib],
#'   [`distrib_gradient()`][distrib_gradient.BetaBinom1Distrib],
#'   [`distrib_hessian()`][distrib_hessian.BetaBinom1Distrib],
#'   [`distrib_pdf()`][distrib_pdf.BetaBinom1Distrib],
#'   [`distrib_quantile()`][distrib_quantile.BetaBinom1Distrib],
#'   [`distrib_rng()`][distrib_rng.BetaBinom1Distrib].
#'
#' Six more are registered on the class from other files: the closed-form
#'   [`distrib_deriv3()`][distrib_deriv3.BetaBinom1Distrib] and
#'   [`distrib_deriv4()`][distrib_deriv4.BetaBinom1Distrib] in
#'   `betabinom1_higher.R`, and the four moments
#'   [`mean()`][mean.BetaBinom1Distrib], [`variance()`][variance.BetaBinom1Distrib],
#'   [`skewness()`][skewness.BetaBinom1Distrib] and
#'   [`kurtosis()`][kurtosis.BetaBinom1Distrib] in `moments.R`.
#'
#' The response derivatives are refused, as for every discrete family: a mass
#' function has no derivative in its argument.
#'
#' @examples
#' d <- betabinom1_distrib(size = 10)
#' S7::S7_inherits(d, discrete_distrib)
#'
#' # The trial count is a property of the object, not an entry of theta.
#' d@size
#' d@bounds
#'
#' # The properties a consumer reads to drive the family without knowing it.
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
#' while their difference is of order one, so writing the mass as that
#' difference loses every digit at a small \eqn{\sigma}. The compiled kernel
#' switches to a sum of logarithms instead, the shifts being integers, and
#' stays exact to the binomial limit.
#'
#' @param distrib A `BetaBinom1Distrib` object, from [betabinom1_distrib()].
#' @param y A numeric vector of counts. A value outside
#'   \eqn{\{0, \dots, n\}} or not an integer gives a mass of 0, or `-Inf` with
#'   `log = TRUE`.
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
#' \eqn{\{0, \dots, n\}}. The support being finite, that sum is exact rather
#' than an approximation, and it costs \eqn{n+1} evaluations of the mass
#' whatever the length of `q`. Below the support the probability is 0 and at or
#' above \eqn{n} it is 1.
#'
#' @param distrib A `BetaBinom1Distrib` object, from [betabinom1_distrib()].
#' @param q A numeric vector of quantiles. A non-integer value is floored, so
#'   `F(2.7)` is `F(2)`.
#' @param theta A named list with components `mu` and `sigma`, each a numeric
#'   vector of length 1. `mu` must lie in \eqn{(0, 1)} and `sigma` be strictly
#'   positive. A parameter varying by observation is not supported here: one
#'   cumulative table is built for the whole call.
#' @param lower.tail Logical of length 1. When `TRUE`, the default,
#'   probabilities are \eqn{P(Y \le q)}; when `FALSE` they are \eqn{P(Y > q)},
#'   formed as \eqn{1 - F}.
#' @param log.p Logical of length 1. When `TRUE` the logarithm of the
#'   probability is returned. Defaults to `FALSE`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of probabilities in \eqn{[0, 1]}, of length
#'   `length(q)`. With `log.p = TRUE` the values are logarithms and are
#'   non-positive.
#'
#' @seealso [distrib_quantile.BetaBinom1Distrib()] for the generalized
#'   inverse, [distrib_pdf.BetaBinom1Distrib()] for the mass, and
#'   [distrib_cdf()] for the generic.
#'
#' @examples
#' d <- betabinom1_distrib(size = 10)
#' th <- list(mu = 0.3, sigma = 0.5)
#'
#' # The cumulative sum of the mass, exactly.
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
  supp <- 0:n
  mass <- exp(betabinom_logpmf_cpp(supp, theta[[1]], theta[[2]], n))
  cum <- cumsum(mass)
  k <- floor(q)
  p <- ifelse(k < 0, 0, ifelse(k >= n, 1, cum[pmin(pmax(k, 0), n) + 1L]))
  p <- pmin(pmax(p, 0), 1)
  if (!lower.tail) p <- 1 - p
  if (log.p) log(p) else p
}

#' @title Beta-Binomial Quantile Function
#' @name distrib_quantile.BetaBinom1Distrib
#' @description
#' Computes the generalized inverse \eqn{Q(p) = \min\{y : F(y) \ge p\}} by
#' walking the exact cumulative sum over the support. The result is an integer
#' count, the distribution function being a step function, so `Q(F(y))` returns
#' `y` while `F(Q(p))` is at least `p` and generally exceeds it. A tolerance of
#' `1e-12` is allowed on the comparison so that a `p` obtained from
#' [distrib_cdf.BetaBinom1Distrib()] maps back to the count it came from.
#'
#' @param distrib A `BetaBinom1Distrib` object, from [betabinom1_distrib()].
#' @param p A numeric vector of probabilities in \eqn{[0, 1]}, or of their
#'   logarithms when `log.p = TRUE`. `NA` is returned for `NA`; a value at or
#'   below 0 gives 0 and one at or above 1 gives `n`.
#' @param theta A named list with components `mu` and `sigma`, each a numeric
#'   vector of length 1. `mu` must lie in \eqn{(0, 1)} and `sigma` be strictly
#'   positive. One cumulative table is built for the whole call, so a parameter
#'   varying by observation is not supported here.
#' @param lower.tail Logical of length 1. When `TRUE`, the default, `p` is
#'   \eqn{P(Y \le q)}; when `FALSE` it is \eqn{P(Y > q)}.
#' @param log.p Logical of length 1. When `TRUE`, `p` is read as a logarithm.
#'   Defaults to `FALSE`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of counts in \eqn{\{0, \dots, n\}}, of length
#'   `length(p)`.
#'
#' @seealso [distrib_cdf.BetaBinom1Distrib()] for the function inverted here,
#'   [distrib_rng.BetaBinom1Distrib()] for draws, and [distrib_quantile()] for
#'   the generic.
#'
#' @examples
#' d <- betabinom1_distrib(size = 10)
#' th <- list(mu = 0.3, sigma = 0.5)
#'
#' # The deciles, which are counts.
#' distrib_quantile(d, c(0.1, 0.5, 0.9), th)
#'
#' # The round trip from a count through F and back is the identity.
#' all.equal(distrib_quantile(d, distrib_cdf(d, 0:10, th), th), as.numeric(0:10))
#'
#' # The other direction only reaches at least p, F being a step function.
#' rbind(p = c(0.2, 0.45, 0.8),
#'       F_at_Q = distrib_cdf(d, distrib_quantile(d, c(0.2, 0.45, 0.8), th), th))
S7::method(distrib_quantile, BetaBinom1Distrib) <- function(distrib, p, theta,
                                                            lower.tail = TRUE,
                                                            log.p = FALSE, ...) {
  if (log.p) p <- exp(p)
  if (!lower.tail) p <- 1 - p
  n <- distrib@size
  cum <- cumsum(exp(betabinom_logpmf_cpp(0:n, theta[[1]], theta[[2]], n)))
  vapply(p, function(pp) {
    if (is.na(pp)) return(NA_real_)
    if (pp <= 0) return(0)
    if (pp >= 1) return(n)
    (0:n)[which(cum >= pp - 1e-12)[1L]]
  }, numeric(1))
}

#' @title Beta-Binomial Random Generation
#' @name distrib_rng.BetaBinom1Distrib
#' @description
#' Draws `n` independent beta-binomial counts by the two-stage hierarchy the
#' family is defined by: a success probability from
#' \eqn{\mathrm{Beta}(\alpha, \beta)} with \eqn{\alpha = \mu/\sigma} and
#' \eqn{\beta = (1-\mu)/\sigma}, then a count from
#' \eqn{\mathrm{Binomial}(n_{\mathrm{trials}}, p)} at that probability, one
#' fresh probability per draw. The draws depend on `.Random.seed` in the usual
#' way and consume two of R's streams per variate.
#'
#' @param distrib A `BetaBinom1Distrib` object, from [betabinom1_distrib()].
#' @param n A single positive integer, the number of draws. Note that the
#'   number of **trials** is the object's `size` property, not this argument.
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
#' Each digamma difference cancels to leading order as the shapes grow, so the
#' compiled kernel forms it as a sum of reciprocals, the shifts being integers.
#' The location component then holds its accuracy to a concentration of
#' \eqn{10^{15}}, where the direct difference has lost every digit.
#'
#' With `scale = "link"` the generic applies the chain rule for the links the
#' family carries before returning. This method always returns the parameter
#' scale. The arithmetic runs in a compiled kernel decomposed over the elements
#' of the output, so the result does not depend on the thread count.
#'
#' @param distrib A `BetaBinom1Distrib` object, from [betabinom1_distrib()].
#' @param y A numeric vector of counts in \eqn{\{0, \dots, n\}}.
#' @param theta A named list with components `mu` and `sigma`, each a numeric
#'   vector of length 1 or of the length of `y`. A component of length 1 is
#'   recycled. `mu` must lie in \eqn{(0, 1)} and `sigma` be strictly positive.
#' @param scale One of `"parameter"` (the default) or `"link"`, matched by
#'   [base::match.arg()]. Read by the generic, not by this method.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#' @param threads A single positive integer, how many threads the kernel may
#'   use. Below the measured internal threshold the kernel stays sequential
#'   whatever the count says. Defaults to `1L`.
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
#' and the **mixed** shape component carries only the \eqn{S} part,
#' \eqn{-\psi_1(n+S) + \psi_1(S)}, the two shapes entering the mass separately
#' otherwise. The reported components follow by the two-variable chain rule of
#' \eqn{(\alpha, \beta) = (\mu/\sigma, (1-\mu)/\sigma)}, whose own second
#' derivatives contribute the rest.
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
#'   [base::match.arg()]. Read by the generic, not by this method.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#' @param threads A single positive integer, how many threads the kernel may
#'   use. Below the measured internal threshold the kernel stays sequential
#'   whatever the count says. Defaults to `1L`.
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
#' evaluated at \eqn{y = k}. The support being finite, that is an **exact
#' finite sum** of \eqn{n+1} terms rather than a quadrature or a sample. The
#' answer is the expectation to machine precision, and a bounded support is
#' what makes that available.
#'
#' `approx` and `nsim` are therefore ignored; every strategy returns the same
#' three numbers. The arithmetic runs in a compiled kernel decomposed over the
#' elements of the output, so the result does not depend on the thread count.
#'
#' @param distrib A `BetaBinom1Distrib` object, from [betabinom1_distrib()].
#' @param y A numeric vector of counts. Only its length is read, the
#'   expectation not depending on the data; the values are ignored.
#' @param theta A named list with components `mu` and `sigma`, each a numeric
#'   vector of length 1 or of the length of `y`. A component of length 1 is
#'   recycled. `mu` must lie in \eqn{(0, 1)} and `sigma` be strictly positive.
#' @param scale One of `"parameter"` (the default) or `"link"`, matched by
#'   [base::match.arg()]. Read by the generic, not by this method.
#' @param approx Ignored here, the expectation being an exact sum. Accepted so
#'   that the signature matches the generic's, where it selects between
#'   `"bartlett"`, `"integrate"`, `"mc"` and `"opg"`.
#' @param nsim Ignored here, for the same reason. Defaults to `10000`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#' @param threads A single positive integer, how many threads the kernel may
#'   use. Below the measured internal threshold the kernel stays sequential
#'   whatever the count says. Defaults to `1L`.
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
#' # The strategy argument is inert, the expectation being an exact sum.
#' identical(eh, distrib_expected_hessian(d, 0:10, th, approx = "mc",
#'                                        nsim = 50))
S7::method(distrib_expected_hessian, BetaBinom1Distrib) <- function(distrib, y, theta,
                                                                    scale = c("parameter", "link"),
                                                                    approx = c("opg", "bartlett", "integrate", "mc"),
                                                                    nsim = 10000, ...,
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
#' that are exact finite sums over the support, so every generic of the toolkit
#' answers without a quadrature.
#'
#' The family is the binomial with its success probability drawn from a beta,
#' which makes it the standard model for a proportion whose trials are not
#' independent. It is overdispersed relative to a binomial at every
#' \eqn{\sigma > 0} and approaches one as \eqn{\sigma \to 0}.
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
#' The family is **not** reachable from anything else in the package. It is
#' neither a binomial with a parameter held fixed nor a wrapper over one, the
#' mixing being over the success probability and not over the outcome.
#'
#' # Overdispersion
#'
#' The mean is \eqn{n\mu} and the variance
#' \deqn{\operatorname{Var}(Y) = n\mu(1-\mu)
#'       \left(1 + (n-1)\dfrac{\sigma}{1+\sigma}\right),}
#' so the inflation over a binomial of the same mean is
#' \eqn{1 + (n-1)\sigma/(1+\sigma)}, always above 1 and rising with both the
#' trial count and the dispersion. The factor \eqn{\sigma/(1+\sigma)} is the
#' **intraclass correlation** between two trials of the same group. At
#' \eqn{n = 10} and \eqn{\sigma = 0.5} that correlation is \eqn{1/3} and the
#' variance is four times a binomial's, 8.4 against 2.1.
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
#' **closed form** as well, in [distrib_deriv3.BetaBinom1Distrib()], and not
#' the numerical fallback.
#'
#' Every expectation is an **exact finite sum** over \eqn{\{0, \dots, n\}}:
#' the support is bounded, so there is nothing to integrate over.
#'
#' # The cancellation at a small dispersion
#'
#' The two beta functions of the mass are of magnitude
#' \eqn{(\alpha+\beta)\log(\alpha+\beta)} while their difference is of order
#' one, so writing the mass as that difference loses one digit per factor of
#' ten in the concentration \eqn{1/\sigma}. At \eqn{\sigma = 10^{-14}} the
#' direct route is wrong in the third decimal of the log-mass. The kernel
#' switches to a sum of logarithms instead, the shifts being integers, and the
#' same rewrite runs through the score: measured at \eqn{n = 10},
#' \eqn{\mu = 0.4}, \eqn{y = 3}, the mean component holds the binomial value
#' \eqn{-25/6} to nine figures at a concentration of \eqn{10^{15}}.
#'
#' # Estimation
#'
#' [fit_distrib()] maximizes the log-likelihood on the link scale. Neither
#' estimate is closed form. A sample with no extra-binomial variation drives
#' \eqn{\hat\sigma} towards zero, which is the correct answer and the boundary
#' of the parameter space.
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
#' # Four times the variance of a binomial with the same mean: the intraclass
#' # correlation is sigma / (1 + sigma) = 1/3, and 1 + 9/3 = 4.
#' c(mean = mean(d, th), var = variance(d, th),
#'   binomial_var = variance(binomial_distrib(size = 10), list(mu = 0.3)),
#'   icc = 0.5 / 1.5)
#'
#' # The dispersion going to zero is the binomial, at rate O(sigma).
#' vapply(c(1e-2, 1e-4, 1e-6), function(s)
#'   max(abs(distrib_pdf(d, 0:10, list(mu = 0.3, sigma = s)) -
#'           dbinom(0:10, 10, 0.3))), numeric(1))
#'
#' # The same law as betabinom2 at the implied shapes.
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
#' unbounded-count analogue, a Poisson mixed over a gamma;
#' [fit_distrib()] to estimate the parameters; [check_distrib()] to validate a
#' family of your own against the same battery this one passes;
#' [BetaBinom1Distrib] for the class.
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
