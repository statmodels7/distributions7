#' @include distrib.R generics.R
NULL

#' @title Negative Binomial Distribution Class, NB1
#' @name NegBin1Distrib
#'
#' @description
#' The S7 class of the negative binomial family on the non-negative integers
#' whose variance is **linear** in the mean: with a mean \eqn{\mu > 0} and a
#' dispersion \eqn{\theta > 0}, \eqn{\operatorname{Var}(Y) = \mu(1+\theta)}, so
#' the variance-to-mean ratio is \eqn{1+\theta} at every mean. It inherits from
#' `discrete_distrib`, so it answers every generic of the `distrib` contract;
#' the seven methods listed below are registered on it in this file and
#' everything else comes from the parent.
#'
#' Build one with [negbin1_distrib()], which supplies the two link functions
#' and fills the properties in. This page documents the raw S7 constructor,
#' which takes the parent's properties and validates none of the relationships
#' between them.
#'
#' @inheritParams distrib
#'
#' @return An S7 object of class `NegBin1Distrib`, inheriting from
#'   `discrete_distrib` and from `distrib`. Its properties are the parent's:
#'   `distrib_name`, `dimension`, `bounds`, `params`, `params_interpretation`,
#'   `n_params`, `params_bounds`, `link_params` and `params_smooth`. For an
#'   object built by [negbin1_distrib()] they hold `"negbin1"`,
#'   `"univariate"`, `c(0, Inf)`, `c("mu", "theta")`, the interpretations
#'   `c(mu = "mean", theta = "dispersion")`, `2`, the domain
#'   \eqn{(0, \infty)} for both parameters, and the two links.
#'
#' @seealso [negbin1_distrib()] to build one;
#'   [negbin2_distrib()] for the quadratic-variance family, which is a
#'   different family and not a reparametrization of this one;
#'   [poisson_distrib()] for the limit as \eqn{\theta} goes to zero;
#'   [distrib_pdf.NegBin1Distrib()] and [distrib_gradient.NegBin1Distrib()] for
#'   the closed forms this class supplies.
#'
#' @section Methods:
#' Registered on this class in this file:
#'   [`distrib_cdf()`][distrib_cdf.NegBin1Distrib],
#'   [`distrib_expected_hessian()`][distrib_expected_hessian.NegBin1Distrib],
#'   [`distrib_gradient()`][distrib_gradient.NegBin1Distrib],
#'   [`distrib_hessian()`][distrib_hessian.NegBin1Distrib],
#'   [`distrib_pdf()`][distrib_pdf.NegBin1Distrib],
#'   [`distrib_quantile()`][distrib_quantile.NegBin1Distrib],
#'   [`distrib_rng()`][distrib_rng.NegBin1Distrib]
#'
#' The third and fourth orders,
#'   [`distrib_deriv3()`][distrib_deriv3.NegBin1Distrib] and
#'   [`distrib_deriv4()`][distrib_deriv4.NegBin1Distrib], are registered
#' elsewhere in the package, on top of [negbin1_components()]. A discrete
#' family has no derivatives in the response. Everything else is inherited from
#' [discrete_distrib()].
#'
#' @examples
#' d <- negbin1_distrib()
#' S7::S7_inherits(d, discrete_distrib)
#'
#' # The properties a consumer reads to drive the family without knowing it.
#' d@params
#' d@params_interpretation
#' d@bounds
#'
#' # The variance-to-mean ratio is 1 + theta at every mean, which is what
#' # separates this family from the quadratic one.
#' vapply(c(1, 10, 100),
#'        function(m) variance(d, list(mu = m, theta = 4)) / m, numeric(1))
NegBin1Distrib <- S7::new_class("NegBin1Distrib", parent = discrete_distrib)

#' The Size Behind an NB1 Mean
#'
#' @description
#' Returns \eqn{r = \mu/\theta}, the number of successes the base R negative
#' binomial functions take as `size`. Requiring the variance to be
#' \eqn{\mu(1+\theta)} fixes the success probability at \eqn{1/(1+\theta)}, and
#' the mean then determines the size. It is this that puts \eqn{\mu} inside the
#' gamma functions of the mass and makes NB1 a different family from the
#' quadratic-variance one.
#'
#' Nothing is validated; a non-positive `theta` gives an infinite or negative
#' size and the caller sees the failure at the base R function.
#'
#' @param mu The mean, a positive numeric vector.
#' @param theta The dispersion, a positive numeric vector. The two are used
#'   elementwise, so vectors of different lengths recycle in the usual way.
#'
#' @return A numeric vector of sizes, of the length the division produces.
#'
#' @seealso [negbin1_distrib()] for the family, [nb1_prob()] for the success
#'   probability that goes with it, and [distrib_pdf.NegBin1Distrib()] for the
#'   mass they feed.
#'
#' @keywords internal
nb1_size <- function(mu, theta) mu / theta

#' The Success Probability Behind an NB1 Dispersion
#'
#' @description
#' Returns \eqn{1/(1+\theta)}, the success probability the base R negative
#' binomial functions take as `prob`. It is the value that makes the
#' variance-to-mean ratio \eqn{1+\theta} at every mean, so the dispersion
#' relative to a Poisson does not change with the level of the counts.
#'
#' Nothing is validated; a `theta` at or below \eqn{-1} gives a value outside
#' \eqn{(0, 1)} and the caller sees the failure at the base R function.
#'
#' @param theta The dispersion, a positive numeric vector.
#'
#' @return A numeric vector in \eqn{(0, 1)}, of the length of `theta`, falling
#'   towards 0 as the dispersion grows and to 1 as it goes to zero.
#'
#' @seealso [negbin1_distrib()] for the family, [nb1_size()] for the size that
#'   goes with it, and [distrib_pdf.NegBin1Distrib()] for the mass they feed.
#'
#' @keywords internal
nb1_prob <- function(theta) 1 / (1 + theta)

# --- S7 METHODS IMPLEMENTATION ---

#' @title NB1 Probability Mass Function
#' @name distrib_pdf.NegBin1Distrib
#' @description
#' Computes the negative binomial mass at size \eqn{r = \mu/\theta} and success
#' probability \eqn{1/(1+\theta)},
#' \deqn{P(Y = y; \mu, \theta) = \dfrac{\Gamma(y+r)}{\Gamma(r)\,y!}
#'       \left(\dfrac{1}{1+\theta}\right)^{r}
#'       \left(\dfrac{\theta}{1+\theta}\right)^{y}, \qquad
#'       y = 0, 1, 2, \ldots,}
#' which is the pairing that makes the variance \eqn{\mu(1+\theta)}. The
#' log-mass is formed in a compiled kernel and exponentiated when `log` is
#' `FALSE`, so `log = TRUE` is the accurate route far into the tail.
#'
#' The mean sits **inside** the gamma functions here, through the size. In
#' [negbin2_distrib()] the size is \eqn{\theta} and the mean stays outside
#' them, which is the concrete difference between the two families.
#'
#' @param distrib A `NegBin1Distrib` object, from [negbin1_distrib()].
#' @param y A numeric vector of counts. A negative value gives a mass of 0.
#' @param theta A named list with components `mu` and `theta`, each a numeric
#'   vector of length 1 or of the length of `y`. A component of length 1 is
#'   recycled. Both must be strictly positive.
#' @param log Logical of length 1. When `TRUE` the log-mass is returned.
#'   Defaults to `FALSE`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of probabilities, of length
#'   `max(length(y), length(mu), length(theta))`, one value per observation.
#'
#' @seealso [distrib_cdf.NegBin1Distrib()] for the distribution function,
#'   [distrib_gradient.NegBin1Distrib()] for the derivatives of the log-mass,
#'   [nb1_size()] and [nb1_prob()] for the pairing this uses, and
#'   [distrib_pdf()] for the generic.
#'
#' @examples
#' d <- negbin1_distrib()
#' y <- c(0, 2, 6)
#' th <- list(mu = 4, theta = 4)
#'
#' # The mass is stats::dnbinom at size mu/theta and prob 1/(1 + theta).
#' all.equal(distrib_pdf(d, y, th), dnbinom(y, size = 1, prob = 1 / 5))
#'
#' # At mu = theta the size is 1 and the family is the geometric.
#' all.equal(distrib_pdf(d, 0:4, th), dgeom(0:4, prob = 1 / 5))
#'
#' # As theta goes to zero the mass tends to the Poisson's.
#' rbind(nb1 = distrib_pdf(d, 0:4, list(mu = 3, theta = 1e-6)),
#'       poisson = dpois(0:4, 3))
#'
#' # The mass sums to one over the support.
#' sum(distrib_pdf(d, 0:600, th))
S7::method(distrib_pdf, NegBin1Distrib) <- function(distrib, y, theta, log = FALSE, ...) {
  out <- negbin1_logpmf_cpp(y, theta[[1]], theta[[2]])
  if (log) out else exp(out)
}

#' @title NB1 Cumulative Distribution Function
#' @name distrib_cdf.NegBin1Distrib
#' @description
#' Computes the negative binomial distribution function, the partial sum of the
#' mass, by calling [stats::pnbinom()] at size \eqn{r = \mu/\theta} and success
#' probability \eqn{1/(1+\theta)}. That function evaluates it through the
#' incomplete beta function, so nothing is summed. The result is a step
#' function, constant between consecutive integers. Both tails are available
#' exactly, and `log.p = TRUE` returns a logarithm that stays finite where the
#' probability itself underflows.
#'
#' @param distrib A `NegBin1Distrib` object, from [negbin1_distrib()].
#' @param q A numeric vector of quantiles. A non-integer value is floored, and
#'   a value below zero gives a lower-tail probability of 0.
#' @param theta A named list with components `mu` and `theta`, each a numeric
#'   vector of length 1 or of the length of `q`. A component of length 1 is
#'   recycled. Both must be strictly positive.
#' @param lower.tail Logical of length 1. When `TRUE`, the default,
#'   probabilities are \eqn{P(Y \le q)}; when `FALSE` they are \eqn{P(Y > q)}.
#' @param log.p Logical of length 1. When `TRUE` the logarithm of the
#'   probability is returned. Defaults to `FALSE`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of probabilities in \eqn{[0, 1]}, of length
#'   `max(length(q), length(mu), length(theta))`. With `log.p = TRUE` the
#'   values are logarithms and are non-positive.
#'
#' @seealso [distrib_quantile.NegBin1Distrib()] for the generalized inverse,
#'   [distrib_pdf.NegBin1Distrib()] for the mass, [nb1_size()] and
#'   [nb1_prob()] for the pairing this uses, and [distrib_cdf()] for the
#'   generic.
#'
#' @examples
#' d <- negbin1_distrib()
#' th <- list(mu = 4, theta = 4)
#'
#' # The method is stats::pnbinom at size mu/theta and prob 1/(1 + theta).
#' all.equal(distrib_cdf(d, c(0, 2, 6), th),
#'           pnbinom(c(0, 2, 6), size = 1, prob = 1 / 5))
#'
#' # A step function: it does not move between two consecutive integers.
#' distrib_cdf(d, c(2, 2.5, 2.999), th)
#'
#' # It is the partial sum of the mass.
#' all.equal(distrib_cdf(d, 6, th), sum(distrib_pdf(d, 0:6, th)))
#'
#' # Far in the upper tail the probability underflows and its log does not.
#' distrib_cdf(d, 2000, th, lower.tail = FALSE)
#' distrib_cdf(d, 2000, th, lower.tail = FALSE, log.p = TRUE)
S7::method(distrib_cdf, NegBin1Distrib) <- function(distrib, q, theta,
                                                     lower.tail = TRUE,
                                                     log.p = FALSE, ...) {
  stats::pnbinom(q, size = nb1_size(theta[[1]], theta[[2]]),
                 prob = nb1_prob(theta[[2]]),
                 lower.tail = lower.tail, log.p = log.p)
}

#' @title NB1 Quantile Function
#' @name distrib_quantile.NegBin1Distrib
#' @description
#' Computes the generalized inverse of the distribution function,
#' \deqn{Q(p; \mu, \theta) =
#'       \min\left\{y \in \mathbb{N}_0 : F(y; \mu, \theta) \ge p\right\},}
#' by calling [stats::qnbinom()] at size \eqn{r = \mu/\theta} and success
#' probability \eqn{1/(1+\theta)}. The support is a lattice, so the round trip
#' through [distrib_cdf.NegBin1Distrib()] does **not** return `p`: it returns
#' the mass up to the smallest integer whose cumulative probability reaches
#' `p`, which overshoots.
#'
#' @param distrib A `NegBin1Distrib` object, from [negbin1_distrib()].
#' @param p A numeric vector of probabilities in \eqn{[0, 1]}, or of their
#'   logarithms when `log.p = TRUE`. A value outside the range gives `NaN` with
#'   a warning; `p = 1` gives `Inf`.
#' @param theta A named list with components `mu` and `theta`, each a numeric
#'   vector of length 1 or of the length of `p`. A component of length 1 is
#'   recycled. Both must be strictly positive.
#' @param lower.tail Logical of length 1. When `TRUE`, the default, `p` is
#'   \eqn{P(Y \le q)}; when `FALSE` it is \eqn{P(Y > q)}.
#' @param log.p Logical of length 1. When `TRUE` the values in `p` are read as
#'   logarithms of probabilities. Defaults to `FALSE`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of non-negative integers, of length
#'   `max(length(p), length(mu), length(theta))`.
#'
#' @seealso [distrib_cdf.NegBin1Distrib()], which this inverts;
#'   [distrib_rng.NegBin1Distrib()], which does not use it;
#'   [distrib_quantile()] for the generic.
#'
#' @examples
#' d <- negbin1_distrib()
#' th <- list(mu = 4, theta = 4)
#'
#' # A nominal central 95 percent interval on the counts.
#' p <- c(0.025, 0.5, 0.975)
#' q <- distrib_quantile(d, p, th)
#' q
#'
#' # On a lattice the round trip overshoots: the cumulative probability at the
#' # returned integer is at or above the one asked for, not equal to it.
#' rbind(asked = p, reached = distrib_cdf(d, q, th))
#'
#' # The overdispersion shows in the width: the same mean under a Poisson
#' # reaches far less far.
#' c(nb1 = distrib_quantile(d, 0.975, th),
#'   poisson = distrib_quantile(poisson_distrib(), 0.975, list(mu = 4)))
S7::method(distrib_quantile, NegBin1Distrib) <- function(distrib, p, theta,
                                                          lower.tail = TRUE,
                                                          log.p = FALSE, ...) {
  stats::qnbinom(p, size = nb1_size(theta[[1]], theta[[2]]),
                 prob = nb1_prob(theta[[2]]),
                 lower.tail = lower.tail, log.p = log.p)
}

#' @title NB1 Random Number Generator
#' @name distrib_rng.NegBin1Distrib
#' @description
#' Draws `n` independent counts by calling [stats::rnbinom()] at size
#' \eqn{r = \mu/\theta} and success probability \eqn{1/(1+\theta)}, so the
#' draws come from R's own generator and depend on `.Random.seed` in the usual
#' way. The cumulative-table fallback the base class supplies for a discrete
#' family is bypassed.
#'
#' @param distrib A `NegBin1Distrib` object, from [negbin1_distrib()].
#' @param n A single positive integer, the number of draws.
#' @param theta A named list with components `mu` and `theta`, each a numeric
#'   vector of length 1 or of length `n`. A component of length 1 is recycled,
#'   so a vector of length `n` draws one count per parameter setting. Both must
#'   be strictly positive.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of `n` non-negative integers.
#'
#' @seealso [distrib_quantile.NegBin1Distrib()] for the inverse-transform
#'   route, [fit_distrib()] to estimate the parameters back from a sample, and
#'   [distrib_rng()] for the generic.
#'
#' @examples
#' d <- negbin1_distrib()
#'
#' # Same generator as stats::rnbinom at size mu/theta and prob 1/(1 + theta).
#' set.seed(2)
#' a <- distrib_rng(d, 5, list(mu = 4, theta = 4))
#' set.seed(2)
#' identical(a, rnbinom(5, size = 1, prob = 1 / 5))
#'
#' # The sample moments recover the parameters: the mean directly, and the
#' # dispersion as var/mean - 1.
#' set.seed(5)
#' z <- distrib_rng(d, 2e4, list(mu = 4, theta = 4))
#' c(mu = mean(z), theta = var(z) / mean(z) - 1)
S7::method(distrib_rng, NegBin1Distrib) <- function(distrib, n, theta, ...) {
  stats::rnbinom(n, size = nb1_size(theta[[1]], theta[[2]]),
                 prob = nb1_prob(theta[[2]]))
}

#' @title NB1 Score
#' @name distrib_gradient.NegBin1Distrib
#' @description
#' Computes the first derivatives of the NB1 log-mass with respect to
#' \eqn{\mu} and \eqn{\theta}, one value per observation, in closed form. Both
#' come from the chain rule through the size \eqn{r = \mu/\theta}. Writing
#' \eqn{P = \psi(y+r) - \psi(r) - \log(1+\theta)} with \eqn{\psi} the digamma
#' function,
#' \deqn{\dfrac{\partial \ell}{\partial \mu} = \dfrac{P}{\theta}, \qquad
#'       \dfrac{\partial \ell}{\partial \theta} = -\dfrac{\mu}{\theta^2}P
#'         - \dfrac{r}{1+\theta} + \dfrac{y}{\theta} - \dfrac{y}{1+\theta}.}
#' The arithmetic runs in a compiled kernel decomposed over the elements of the
#' output, so the result does not depend on the thread count.
#'
#' # The dispersion at small theta
#'
#' As \eqn{\theta} goes to zero the family tends to the Poisson, so the
#' dispersion component tends to a finite limit while its individual terms run
#' away: \eqn{r = \mu/\theta} grows without bound and the chain rule divides by
#' \eqn{\theta^2}. The kernel computes the digamma difference in a form that
#' performs the cancellation symbolically, and the value converges onto
#' \deqn{\lim_{\theta \to 0} \dfrac{\partial \ell}{\partial \theta}
#'       = \dfrac{(y-\mu)^2 - y}{2\mu},}
#' which the example below checks at three settings. What survives is about
#' five significant figures at \eqn{\theta = 10^{-8}} and no more, the powers
#' of \eqn{r} in the assembly carrying a cancellation of their own that is not
#' removed.
#'
#' With `scale = "link"` the generic applies the chain rule for the links the
#' family carries before returning. This method always returns the parameter
#' scale.
#'
#' @param distrib A `NegBin1Distrib` object, from [negbin1_distrib()].
#' @param y A numeric vector of counts.
#' @param theta A named list with components `mu` and `theta`, each a numeric
#'   vector of length 1 or of the length of `y`. A component of length 1 is
#'   recycled. Both must be strictly positive.
#' @param scale One of `"parameter"` (the default) or `"link"`, matched by
#'   [base::match.arg()]. Read by the generic, not by this method.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#' @param threads A single positive integer, how many threads the kernel may
#'   use. Below the measured internal threshold the kernel stays sequential
#'   whatever the count says. Defaults to `1L`.
#'
#' @return A named list of two numeric vectors, `mu` and `theta`, each of
#'   length `max(length(y), length(mu), length(theta))`.
#'
#' @section Notation:
#' \eqn{\ell} is the log-mass of one observation, \eqn{\mu > 0} the mean and
#' \eqn{\theta > 0} the dispersion, with
#' \eqn{\operatorname{Var}(Y) = \mu(1+\theta)}. \eqn{r = \mu/\theta} is the
#' size and \eqn{\psi} the digamma function.
#'
#' @seealso [distrib_hessian.NegBin1Distrib()] for the second derivatives,
#'   [distrib_expected_hessian.NegBin1Distrib()] for their expectation and for
#'   what it does at small \eqn{\theta}, [distrib_gradient.NegBin2Distrib()]
#'   for the quadratic-variance family, and [distrib_gradient()] for the
#'   generic.
#'
#' @examples
#' d <- negbin1_distrib()
#' y <- c(0, 2, 6)
#' th <- list(mu = 4, theta = 4)
#' g <- distrib_gradient(d, y, th)
#'
#' # The two closed forms, written out.
#' r <- 4 / 4
#' P <- digamma(y + r) - digamma(r) - log(1 + 4)
#' all.equal(g$mu, P / 4)
#' all.equal(g$theta, -(4 / 4^2) * P - r / (1 + 4) + y / 4 - y / (1 + 4))
#'
#' # As theta goes to zero the dispersion component converges onto
#' # {(y - mu)^2 - y}/(2 mu).
#' lim <- function(y, mu) ((y - mu)^2 - y) / (2 * mu)
#' cmp <- function(y, mu) {
#'   c(vapply(c(1e-4, 1e-6, 1e-8),
#'            function(t) distrib_gradient(d, y, list(mu = mu, theta = t))$theta,
#'            numeric(1)),
#'     limit = lim(y, mu))
#' }
#' rbind(`y=3,mu=4` = cmp(3, 4), `y=0,mu=3` = cmp(0, 3), `y=7,mu=2` = cmp(7, 2))
#'
#' # Summed over a fitted sample the score is at the optimizer's tolerance.
#' set.seed(5)
#' z <- distrib_rng(d, 2000, th)
#' fit <- fit_distrib(d, z)
#' vapply(distrib_gradient(d, z, as.list(coef(fit))), sum, numeric(1))
S7::method(distrib_gradient, NegBin1Distrib) <- function(distrib, y, theta,
                                                          scale = c("parameter", "link"), ...,
                                       threads = 1L) {
  negbin1_gradient_cpp(y, theta[[1]], theta[[2]], threads)
}

#' @title NB1 Observed Hessian
#' @name distrib_hessian.NegBin1Distrib
#' @description
#' Computes the three distinct second derivatives of the NB1 log-mass with
#' respect to \eqn{\mu} and \eqn{\theta}, one value per observation, in closed
#' form. They are the two-variable chain rule through the size
#' \eqn{r = \mu/\theta} applied once more: the second derivative of the
#' log-mass in \eqn{r} is \eqn{\psi'(y+r) - \psi'(r)}, with \eqn{\psi'} the
#' trigamma function, and the term in which \eqn{\theta} appears outside
#' \eqn{r} contributes \eqn{-1/(1+\theta)} to the mixed entry.
#'
#' Because \eqn{\partial r/\partial\mu = 1/\theta} and
#' \eqn{\partial r/\partial\theta = -r/\theta}, every component divides by a
#' power of \eqn{\theta}, and the components in \eqn{\theta} lose their digits
#' at small \eqn{\theta} for the reason given on
#' [distrib_gradient.NegBin1Distrib()].
#'
#' @param distrib A `NegBin1Distrib` object, from [negbin1_distrib()].
#' @param y A numeric vector of counts.
#' @param theta A named list with components `mu` and `theta`, each a numeric
#'   vector of length 1 or of the length of `y`. A component of length 1 is
#'   recycled. Both must be strictly positive.
#' @param scale One of `"parameter"` (the default) or `"link"`, matched by
#'   [base::match.arg()]. Read by the generic, not by this method.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#' @param threads A single positive integer, how many threads the kernel may
#'   use. Below the measured internal threshold the kernel stays sequential
#'   whatever the count says. Defaults to `1L`.
#'
#' @return A named list of three numeric vectors, `mu_mu`, `mu_theta` and
#'   `theta_theta`, in that order, each of length
#'   `max(length(y), length(mu), length(theta))`. The three name the distinct
#'   entries of a symmetric \eqn{2 \times 2} matrix per observation.
#'
#' @section Notation:
#' \eqn{\ell^{(ij)}} is the second derivative of the log-mass in parameters
#' \eqn{i} and \eqn{j}; parenthesized superscripts name derivatives.
#' \eqn{r = \mu/\theta} is the size, \eqn{\psi} the digamma function and
#' \eqn{\psi'} the trigamma.
#'
#' @seealso [distrib_gradient.NegBin1Distrib()] for the score,
#'   [distrib_expected_hessian.NegBin1Distrib()] for the expectation of this
#'   quantity, [distrib_deriv3.NegBin1Distrib()] for the order above, and
#'   [distrib_hessian()] for the generic.
#'
#' @examples
#' d <- negbin1_distrib()
#' y <- c(0, 2, 6)
#' th <- list(mu = 4, theta = 4)
#' h <- distrib_hessian(d, y, th)
#' h
#'
#' # The curvature in mu is exactly zero at y = 0, the size being 1 there, so
#' # the observed information is singular at that count.
#' h$mu_mu
#'
#' # It is the second derivative of the log-mass, so a central difference of
#' # the score reproduces it.
#' eps <- 1e-6
#' up <- distrib_gradient(d, y, list(mu = 4 + eps, theta = 4))$mu
#' dn <- distrib_gradient(d, y, list(mu = 4 - eps, theta = 4))$mu
#' all.equal((up - dn) / (2 * eps), h$mu_mu, tolerance = 1e-5)
#'
#' # And the mixed entry against a difference of the other component.
#' up <- distrib_gradient(d, y, list(mu = 4, theta = 4 + eps))$mu
#' dn <- distrib_gradient(d, y, list(mu = 4, theta = 4 - eps))$mu
#' all.equal((up - dn) / (2 * eps), h$mu_theta, tolerance = 1e-5)
S7::method(distrib_hessian, NegBin1Distrib) <- function(distrib, y, theta,
                                                         scale = c("parameter", "link"), ...,
                                       threads = 1L) {
  negbin1_hessian_cpp(y, theta[[1]], theta[[2]], threads)
}

#' @title NB1 Expected Hessian
#' @name distrib_expected_hessian.NegBin1Distrib
#' @description
#' Returns the expectation of the observed Hessian under the model. Every term
#' carrying \eqn{P = \psi(y+r) - \psi(r) - \log(1+\theta)} drops out, its
#' expectation vanishing by the first Bartlett identity, and what remains needs
#' only \eqn{\mathbb{E}[\psi'(Y+r)]}. That has no closed form: it is summed
#' against the exact mass out to a far-tail quantile, so this is a truncated
#' exact sum and not a quadrature or a simulation. `approx` and `nsim` are
#' ignored, and `y` is read only for its length.
#'
#' **The mixed entry does not vanish**, so the mean and the dispersion are not
#' orthogonal in this family. Measured at four settings it is 0.0172, 0.0364,
#' 0.0108 and 0.0157; in [negbin2_distrib()] the same entry is exactly zero.
#' That is a difference between the two negative binomials rather than a
#' difference of parametrization.
#'
#' # A caveat at small theta
#'
#' The dispersion entry inherits the cancellation of
#' [distrib_gradient.NegBin1Distrib()] and **is not rewritten to remove it**.
#' The chain rule through \eqn{r = \mu/\theta} divides by \eqn{\theta^2} and
#' \eqn{\theta^4}, so what is left of the digits runs out early. Measured at
#' \eqn{\mu = 4}: the entry reads \eqn{-0.489} at \eqn{\theta = 10^{-2}} and
#' \eqn{-0.500} at \eqn{10^{-4}}, then \eqn{+2.1\times 10^{2}} at \eqn{10^{-6}}
#' and \eqn{+2.9\times 10^{8}} at \eqn{10^{-8}}. The sign is impossible for an
#' expected second derivative, and the matrix is indefinite there, its
#' determinant turning negative. A Fisher scoring step taken in that regime is
#' not reliable, and a nearly equidispersed sample drives a fit into it.
#'
#' @param distrib A `NegBin1Distrib` object, from [negbin1_distrib()].
#' @param y A numeric vector of counts. Only its length is used.
#' @param theta A named list with components `mu` and `theta`, each a numeric
#'   vector of length 1 or of the length of `y`. A component of length 1 is
#'   recycled. Both must be strictly positive.
#' @param scale One of `"parameter"` (the default) or `"link"`, matched by
#'   [base::match.arg()]. Read by the generic, not by this method.
#' @param approx Ignored. The surviving expectation is a truncated exact sum.
#'   Accepted so that the signature matches the generic's.
#' @param nsim Ignored, for the same reason. Defaults to `10000`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#' @param threads A single positive integer, how many threads the kernel may
#'   use. Below the measured internal threshold the kernel stays sequential
#'   whatever the count says. Defaults to `1L`.
#'
#' @return A named list of three numeric vectors, `mu_mu`, `mu_theta` and
#'   `theta_theta`, in that order, each of length
#'   `max(length(y), length(mu), length(theta))` and constant within itself
#'   when the parameters are.
#'
#' @section Notation:
#' The **expected information** is
#' \eqn{\mathbb{E}[-\partial^2\ell/\partial\theta\,\partial\theta^\top]}, the
#' expectation of the **observed information** under the model. The first
#' Bartlett identity is \eqn{\mathbb{E}[\partial\ell/\partial\theta] = 0}.
#'
#' @seealso [distrib_hessian.NegBin1Distrib()] for the observed quantity this
#'   is the expectation of, [distrib_expected_hessian.NegBin2Distrib()] for the
#'   quadratic-variance family, whose mixed entry is exactly zero,
#'   [fisher_scoring()], which inverts it at each step, and
#'   [distrib_expected_hessian()] for the generic.
#'
#' @examples
#' d <- negbin1_distrib()
#' th <- list(mu = 4, theta = 4)
#' e <- distrib_expected_hessian(d, c(0, 2, 6), th)
#' lapply(e, unique)
#'
#' # Negative definite at a moderate dispersion.
#' M <- matrix(c(e$mu_mu[1], e$mu_theta[1], e$mu_theta[1], e$theta_theta[1]), 2)
#' eigen(M, only.values = TRUE)$values
#'
#' # The mixed entry is not zero, so the mean and the dispersion are not
#' # orthogonal here; in the quadratic-variance family it is exactly zero.
#' c(nb1 = e$mu_theta[1],
#'   nb2 = distrib_expected_hessian(negbin2_distrib(), 0, th)$mu_theta)
#'
#' # The dispersion entry loses its digits as theta goes to zero, and turns
#' # positive, which an expected second derivative cannot be.
#' vapply(c(1e-2, 1e-4, 1e-6, 1e-8),
#'        function(t) distrib_expected_hessian(d, 0,
#'                      list(mu = 4, theta = t))$theta_theta,
#'        numeric(1))
S7::method(distrib_expected_hessian, NegBin1Distrib) <- function(distrib, y, theta,
                                                                  scale = c("parameter", "link"),
                                                                  approx = c("bartlett", "integrate", "mc", "opg"),
                                                                  nsim = 10000, ...,
                                       threads = 1L) {
  negbin1_expected_hessian_cpp(y, theta[[1]], theta[[2]], threads)
}

# --- CONSTRUCTOR WRAPPER ---

#' Negative Binomial Distribution, NB1
#'
#' @description
#' Builds the distribution object for the negative binomial family on the
#' non-negative integers whose variance is **linear** in the mean: with a mean
#' \eqn{\mu > 0} and a dispersion \eqn{\theta > 0},
#' \eqn{\operatorname{Var}(Y) = \mu(1+\theta)}. The returned object carries
#' closed-form derivatives of the log-mass to fourth order and closed-form
#' moments.
#'
#' The two arguments choose the links that carry each parameter to the
#' unconstrained scale an optimizer works on. Both default to the logarithm,
#' both parameters being positive.
#'
#' @param link_mu A `link` object from `linkfunctions7` for the mean
#'   \eqn{\mu}. Defaults to [linkfunctions7::log_link()], which maps
#'   \eqn{(0, \infty)} onto the line and so keeps every fitted mean positive.
#' @param link_theta A `link` object from `linkfunctions7` for the dispersion
#'   \eqn{\theta}. Defaults to [linkfunctions7::log_link()], for the same
#'   reason.
#'
#' @details
#' # Two negative binomials, and they are two families
#'
#' Two negative binomials are in common use and they are **different
#' families**, not two parametrizations of one. Here the variance is
#' \eqn{\mu(1+\theta)}, growing in proportion to the mean, so the dispersion
#' relative to a Poisson is the same at every mean; [negbin2_distrib()] has
#' \eqn{\mu + \mu^2/\theta}, growing quadratically. Fitting one is not fitting
#' the other, and a likelihood ratio between them is not a test of nested
#' models.
#'
#' The difference is visible in where the mean sits. The size is
#' \eqn{r = \mu/\theta} and the success probability \eqn{1/(1+\theta)}, so
#' \eqn{\mu} appears **inside** the gamma functions of the mass; in the
#' quadratic form the size is \eqn{\theta} and the mean stays outside them.
#'
#' # The parametrization
#'
#' The mass on \eqn{y = 0, 1, 2, \ldots} is
#' \deqn{P(Y = y) = \dfrac{\Gamma(y + \mu/\theta)}{\Gamma(\mu/\theta)\,y!}
#'       \left(\dfrac{1}{1+\theta}\right)^{\mu/\theta}
#'       \left(\dfrac{\theta}{1+\theta}\right)^{y},}
#' the distribution function is the partial sum and the quantile function its
#' generalized inverse. The mean is \eqn{\mu} and the variance
#' \eqn{\mu(1+\theta)}. Two settings are worth recognizing: at \eqn{\mu =
#' \theta} the size is 1 and the law is the geometric, and as \eqn{\theta} goes
#' to zero it tends to the Poisson, where the quadratic form needs
#' \eqn{\theta \to \infty} instead.
#'
#' # Derivatives
#'
#' Everything follows from the chain rule through the size. With
#' \eqn{P = \psi(y+r) - \psi(r) - \log(1+\theta)} and \eqn{\psi} the digamma
#' function, the score is
#' \deqn{\dfrac{\partial \ell}{\partial \mu} = \dfrac{P}{\theta}, \qquad
#'       \dfrac{\partial \ell}{\partial \theta} = -\dfrac{\mu}{\theta^2}P
#'         - \dfrac{r}{1+\theta} + \dfrac{y}{\theta} - \dfrac{y}{1+\theta},}
#' and the Hessian is the same chain rule at second order. In the expected
#' information every term carrying \eqn{P} drops out, its expectation vanishing
#' by the first Bartlett identity, and only \eqn{\mathbb{E}[\psi'(Y+r)]}
#' remains, which is summed against the exact mass out to a far-tail quantile.
#'
#' **The mean and the dispersion are not orthogonal here.** The mixed entry of
#' the expected information is small but non-zero at every setting measured,
#' where [negbin2_distrib()] has exactly zero there.
#'
#' # What happens as theta goes to zero
#'
#' The family tends to the Poisson, so every derivative in \eqn{\theta} tends
#' to a finite limit while the pieces it is assembled from run away: the size
#' \eqn{r = \mu/\theta} grows without bound and the chain rule divides by
#' powers of \eqn{\theta}. The score's digamma difference is computed in a form
#' that performs its own cancellation symbolically, and the value converges
#' onto \eqn{\{(y-\mu)^2 - y\}/(2\mu)}, holding to about five significant
#' figures at \eqn{\theta = 10^{-8}}. Two quantities are **not** rewritten, and
#' both pages say so:
#'
#' - the expected information in \eqn{\theta}, which turns positive from about
#'   \eqn{\theta = 10^{-6}} and leaves the matrix indefinite
#'   ([distrib_expected_hessian.NegBin1Distrib()]);
#' - the third and fourth derivatives in \eqn{\theta}, whose own cancellation
#'   in the powers of \eqn{r} is untouched
#'   ([distrib_deriv3.NegBin1Distrib()], [distrib_deriv4.NegBin1Distrib()]).
#'
#' A fit reaches that regime routinely: on 2,000 Poisson counts with mean 4 the
#' dispersion is estimated at about \eqn{1.7\times 10^{-8}} and the run reports
#' convergence.
#'
#' # Estimation
#'
#' [fit_distrib()] maximizes the log-likelihood on the link scale. Neither
#' estimate is closed form. The method of moments supplies the starting values
#' \eqn{\hat\mu = \bar y} and \eqn{\hat\theta = s^2/\bar y - 1}, with \eqn{s^2}
#' the sample variance.
#'
#' @section Notation:
#' \eqn{\ell} is the log-mass of one observation, \eqn{\mu > 0} the mean and
#' \eqn{\theta > 0} the dispersion. \eqn{r = \mu/\theta} is the size, \eqn{\psi}
#' the digamma function and \eqn{\psi'} the trigamma. \eqn{\eta} is a parameter
#' on the unconstrained scale of its link, with \eqn{\theta_j = g^{-1}(\eta_j)}.
#'
#' @return An S7 object of class `NegBin1Distrib`, inheriting from
#'   `discrete_distrib`, with `distrib_name` `"negbin1"`, `dimension`
#'   `"univariate"`, `bounds` `c(0, Inf)`, `params` `c("mu", "theta")`,
#'   `n_params` `2`, `params_bounds` the domain \eqn{(0, \infty)} for both, and
#'   `link_params` the two links given here.
#'
#' @seealso
#' [negbin2_distrib()] for the quadratic-variance family; [poisson_distrib()]
#' for the limit as \eqn{\theta} goes to zero and [geometric_distrib()] for the
#' case \eqn{\mu = \theta}; [pig1_distrib()] for a Poisson mixed over an
#' inverse Gaussian; [zero_inflated()] and [zero_adjusted()] for counts with
#' excess zeros; [fit_distrib()] to estimate the parameters; [check_distrib()]
#' to validate a family of your own against the same battery this one passes;
#' [NegBin1Distrib] for the class.
#'
#' @references
#' Cameron, A. C. and Trivedi, P. K. (1986). Econometric models based
#' on count data: comparisons and applications of some estimators and
#' tests. *Journal of Applied Econometrics* **1**, 29-53.
#'
#' @importFrom linkfunctions7 log_link
#' @importFrom stats pnbinom qnbinom rnbinom
#'
#' @examples
#' d <- negbin1_distrib()
#' d
#'
#' # The two negative binomials are different families: at the same (mu, theta)
#' # this one has variance mu(1 + theta) = 20 and the other mu + mu^2/theta = 8.
#' th <- list(mu = 4, theta = 4)
#' c(nb1 = variance(d, th), nb2 = variance(negbin2_distrib(), th))
#'
#' # The variance-to-mean ratio is 1 + theta at every mean.
#' vapply(c(1, 10, 100),
#'        function(m) variance(d, list(mu = m, theta = 4)) / m, numeric(1))
#'
#' # At mu = theta the size is 1 and the law is the geometric.
#' all.equal(distrib_pdf(d, 0:4, th), dgeom(0:4, prob = 1 / 5))
#'
#' # As theta goes to zero it is the Poisson.
#' rbind(nb1 = distrib_pdf(d, 0:4, list(mu = 3, theta = 1e-6)),
#'       poisson = dpois(0:4, 3))
#'
#' # Fitting recovers the parameters; the moment estimates start it off.
#' set.seed(5)
#' z <- distrib_rng(d, 2000, th)
#' fit <- fit_distrib(d, z)
#' rbind(fitted  = coef(fit),
#'       moments = c(mu = mean(z), theta = var(z) / mean(z) - 1))
#'
#' # On an equidispersed sample the dispersion runs to its boundary, which is
#' # the regime the expected information page warns about.
#' set.seed(4)
#' coef(fit_distrib(d, rpois(2000, 4)))
#'
#' @export
negbin1_distrib <- function(link_mu = log_link(), link_theta = log_link()) {
  NegBin1Distrib(
    distrib_name = "negbin1", dimension = "univariate", bounds = c(0, Inf),
    params = c("mu", "theta"),
    params_interpretation = c(mu = "mean", theta = "dispersion"),
    n_params = 2, params_bounds = list(mu = c(0, Inf), theta = c(0, Inf)),
    link_params = list(mu = link_mu, theta = link_theta)
  )
}
