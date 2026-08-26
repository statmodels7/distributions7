#' @include distrib.R generics.R
NULL

#' @title Negative Binomial Distribution Class, NB2
#' @name NegBin2Distrib
#'
#' @description
#' The S7 class of the negative binomial family on the non-negative integers in
#' the NB2 parametrization: the mean \eqn{\mu > 0} and a dispersion
#' \eqn{\theta > 0}, so that \eqn{\operatorname{Var}(Y) = \mu + \mu^2/\theta}
#' and the variance is quadratic in the mean. It inherits from
#' `discrete_distrib`, so it answers every generic of the `distrib` contract;
#' the nine methods listed below are registered on it directly and everything
#' else comes from the parent.
#'
#' Build one with [negbin2_distrib()], which supplies the two link functions
#' and fills the properties in. This page documents the raw S7 constructor,
#' which takes the parent's properties and validates none of the relationships
#' between them.
#'
#' @inheritParams distrib
#'
#' @return An S7 object of class `NegBin2Distrib`, inheriting from
#'   `discrete_distrib` and from `distrib`. Its properties are the parent's:
#'   `distrib_name`, `dimension`, `bounds`, `params`, `params_interpretation`,
#'   `n_params`, `params_bounds`, `link_params` and `params_smooth`. For an
#'   object built by [negbin2_distrib()] they hold `"negbin2"`,
#'   `"univariate"`, `c(0, Inf)`, `c("mu", "theta")`, the interpretations
#'   `c(mu = "mean", theta = "dispersion")`, `2`, the domain
#'   \eqn{(0, \infty)} for both parameters, and the two links.
#'
#' @seealso [negbin2_distrib()] to build one;
#'   [negbin1_distrib()] for the NB1 parametrization, whose variance is linear
#'   in the mean; [poisson_distrib()] for the limit as \eqn{\theta} grows;
#'   [distrib_pdf.NegBin2Distrib()] and [distrib_gradient.NegBin2Distrib()] for
#'   the closed forms this class supplies.
#'
#' @section Methods:
#' Registered on this class:
#'   [`distrib_cdf()`][distrib_cdf.NegBin2Distrib],
#'   [`distrib_deriv3()`][distrib_deriv3.NegBin2Distrib],
#'   [`distrib_deriv4()`][distrib_deriv4.NegBin2Distrib],
#'   [`distrib_expected_hessian()`][distrib_expected_hessian.NegBin2Distrib],
#'   [`distrib_gradient()`][distrib_gradient.NegBin2Distrib],
#'   [`distrib_hessian()`][distrib_hessian.NegBin2Distrib],
#'   [`distrib_pdf()`][distrib_pdf.NegBin2Distrib],
#'   [`distrib_quantile()`][distrib_quantile.NegBin2Distrib],
#'   [`distrib_rng()`][distrib_rng.NegBin2Distrib]
#'
#' A discrete family has no derivatives in the response, so there are nine
#' here where a continuous family has eleven. Everything else is inherited from
#' [discrete_distrib()].
#'
#' @examples
#' d <- negbin2_distrib()
#' S7::S7_inherits(d, discrete_distrib)
#'
#' # The properties a consumer reads to drive the family without knowing it.
#' d@params
#' d@params_interpretation
#' d@bounds
#'
#' # theta is a dispersion read the other way round from a variance: the
#' # smaller it is, the more overdispersed the counts.
#' vapply(c(0.5, 2, 1e6), function(t) variance(d, list(mu = 4, theta = t)),
#'        numeric(1))
NegBin2Distrib <- S7::new_class("NegBin2Distrib", parent = discrete_distrib)

# --- S7 METHODS IMPLEMENTATION ---

#' @title Negative Binomial Probability Mass Function, NB2
#' @name distrib_pdf.NegBin2Distrib
#' @description
#' Computes the negative binomial mass
#' \deqn{P(Y = y; \mu, \theta) = \dfrac{\Gamma(y+\theta)}{y!\,\Gamma(\theta)}
#'       \left(\dfrac{\theta}{\theta+\mu}\right)^{\theta}
#'       \left(\dfrac{\mu}{\theta+\mu}\right)^{y}, \qquad
#'       y = 0, 1, 2, \ldots,}
#' by calling [stats::dnbinom()] at `size = theta` and `mu = mu`, so the
#' accuracy is R's own. With `log = TRUE` the logarithm is formed inside
#' `dnbinom()` and stays finite far into the tail.
#'
#' The mass is that of a Poisson whose rate is itself gamma distributed with
#' mean \eqn{\mu} and shape \eqn{\theta}, which is where the overdispersion
#' comes from. As \eqn{\theta} grows the gamma concentrates and the mass tends
#' to the Poisson's; at \eqn{\theta = 1} it is the geometric.
#'
#' @param distrib A `NegBin2Distrib` object, from [negbin2_distrib()].
#' @param y A numeric vector of counts. A non-integer value gives 0 with a
#'   warning from [stats::dnbinom()], and a negative value gives 0.
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
#' @seealso [distrib_cdf.NegBin2Distrib()] for the distribution function,
#'   [distrib_gradient.NegBin2Distrib()] for the derivatives of the log-mass,
#'   [distrib_pdf.NegBin1Distrib()] for the NB1 parametrization, and
#'   [distrib_pdf()] for the generic.
#'
#' @examples
#' d <- negbin2_distrib()
#' y <- c(0, 2, 6)
#' th <- list(mu = 4, theta = 2)
#'
#' # The method is stats::dnbinom at size = theta, mu = mu.
#' all.equal(distrib_pdf(d, y, th), dnbinom(y, size = 2, mu = 4))
#'
#' # At theta = 1 the family is the geometric.
#' all.equal(distrib_pdf(d, y, list(mu = 4, theta = 1)),
#'           distrib_pdf(geometric_distrib(), y, list(mu = 4)))
#'
#' # At a large theta it is the Poisson, to six figures.
#' rbind(negbin = distrib_pdf(d, 0:4, list(mu = 3, theta = 1e6)),
#'       poisson = dpois(0:4, 3))
#'
#' # The mass sums to one over the support.
#' sum(distrib_pdf(d, 0:400, th))
S7::method(distrib_pdf, NegBin2Distrib) <- function(distrib, y, theta, log = FALSE, ...) {
  stats::dnbinom(
    x = y,
    mu = theta[[1]],
    size = theta[[2]],
    log = log
  )
}

#' @title Negative Binomial Cumulative Distribution Function, NB2
#' @name distrib_cdf.NegBin2Distrib
#' @description
#' Computes the negative binomial distribution function, the partial sum of the
#' mass
#' \deqn{F(q; \mu, \theta) = \sum_{k=0}^{\lfloor q \rfloor}
#'       P(Y = k; \mu, \theta),}
#' by calling [stats::pnbinom()] at `size = theta` and `mu = mu`, which
#' evaluates it through the incomplete beta function rather than by summing.
#' The function is a step function, constant between consecutive integers.
#' Both tails are available exactly, and `log.p = TRUE` returns a logarithm
#' that stays finite where the probability itself underflows.
#'
#' @param distrib A `NegBin2Distrib` object, from [negbin2_distrib()].
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
#' @seealso [distrib_quantile.NegBin2Distrib()] for the generalized inverse,
#'   [distrib_pdf.NegBin2Distrib()] for the mass, [distrib_grad_cdf()] for the
#'   derivatives of this function in the parameters, which are an exact finite
#'   sum for a discrete family, and [distrib_cdf()] for the generic.
#'
#' @examples
#' d <- negbin2_distrib()
#' th <- list(mu = 4, theta = 2)
#'
#' # The method is stats::pnbinom at size = theta, mu = mu.
#' all.equal(distrib_cdf(d, c(0, 2, 6), th), pnbinom(c(0, 2, 6), 2, mu = 4))
#'
#' # A step function: it does not move between two consecutive integers.
#' distrib_cdf(d, c(2, 2.5, 2.999), th)
#'
#' # It is the partial sum of the mass.
#' all.equal(distrib_cdf(d, 6, th), sum(distrib_pdf(d, 0:6, th)))
#'
#' # The two tails sum to one.
#' distrib_cdf(d, 6, th) + distrib_cdf(d, 6, th, lower.tail = FALSE)
#'
#' # Far in the upper tail the probability underflows and its log does not.
#' distrib_cdf(d, 2000, th, lower.tail = FALSE)
#' distrib_cdf(d, 2000, th, lower.tail = FALSE, log.p = TRUE)
S7::method(distrib_cdf, NegBin2Distrib) <- function(distrib, q, theta, lower.tail = TRUE, log.p = FALSE, ...) {
  stats::pnbinom(
    q = q,
    mu = theta[[1]],
    size = theta[[2]],
    lower.tail = lower.tail,
    log.p = log.p
  )
}

#' @title Negative Binomial Quantile Function, NB2
#' @name distrib_quantile.NegBin2Distrib
#' @description
#' Computes the generalized inverse of the distribution function,
#' \deqn{Q(p; \mu, \theta) =
#'       \min\left\{y \in \mathbb{N}_0 : F(y; \mu, \theta) \ge p\right\},}
#' by calling [stats::qnbinom()] at `size = theta` and `mu = mu`. The support
#' is a lattice, so the round trip through
#' [distrib_cdf.NegBin2Distrib()] does **not** return `p`: it returns the mass
#' up to the smallest integer whose cumulative probability reaches `p`, which
#' overshoots.
#'
#' @param distrib A `NegBin2Distrib` object, from [negbin2_distrib()].
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
#' @seealso [distrib_cdf.NegBin2Distrib()], which this inverts;
#'   [distrib_rng.NegBin2Distrib()], which does not use it;
#'   [distrib_quantile()] for the generic.
#'
#' @examples
#' d <- negbin2_distrib()
#' th <- list(mu = 4, theta = 2)
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
#' # The generalized inverse is a step function of p.
#' distrib_quantile(d, c(0.3, 0.4, 0.5), th)
S7::method(distrib_quantile, NegBin2Distrib) <- function(distrib, p, theta, lower.tail = TRUE, log.p = FALSE, ...) {
  stats::qnbinom(
    p = p,
    mu = theta[[1]],
    size = theta[[2]],
    lower.tail = lower.tail,
    log.p = log.p
  )
}

#' @title Negative Binomial Random Number Generator, NB2
#' @name distrib_rng.NegBin2Distrib
#' @description
#' Draws `n` independent negative binomial counts by calling
#' [stats::rnbinom()] at `size = theta` and `mu = mu`, so the draws come from
#' R's own generator and depend on `.Random.seed` in the usual way. The
#' cumulative-table fallback the base class supplies for a discrete family is
#' bypassed.
#'
#' @param distrib A `NegBin2Distrib` object, from [negbin2_distrib()].
#' @param n A single positive integer, the number of draws.
#' @param theta A named list with components `mu` and `theta`, each a numeric
#'   vector of length 1 or of length `n`. A component of length 1 is recycled,
#'   so a vector of length `n` draws one count per parameter setting. Both must
#'   be strictly positive.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of `n` non-negative integers.
#'
#' @seealso [distrib_quantile.NegBin2Distrib()] for the inverse-transform
#'   route, [fit_distrib()] to estimate the parameters back from a sample, and
#'   [distrib_rng()] for the generic.
#'
#' @examples
#' d <- negbin2_distrib()
#'
#' # Same generator as stats::rnbinom at size = theta, mu = mu.
#' set.seed(2)
#' a <- distrib_rng(d, 5, list(mu = 4, theta = 2))
#' set.seed(2)
#' identical(a, rnbinom(5, size = 2, mu = 4))
#'
#' # The sample moments recover the parameters: the mean directly, and the
#' # dispersion as mu^2/(var - mu).
#' set.seed(7)
#' z <- distrib_rng(d, 2e4, list(mu = 4, theta = 2))
#' c(mu = mean(z), theta = mean(z)^2 / (var(z) - mean(z)))
#'
#' # Overdispersion is the point: the variance is three times the mean here.
#' c(mean = mean(z), var = var(z))
S7::method(distrib_rng, NegBin2Distrib) <- function(distrib, n, theta, ...) {
  stats::rnbinom(
    n = n,
    mu = theta[[1]],
    size = theta[[2]]
  )
}

#' @title Negative Binomial Score, NB2
#' @name distrib_gradient.NegBin2Distrib
#' @description
#' Computes the first derivatives of the negative binomial log-mass with
#' respect to \eqn{\mu} and \eqn{\theta}, one value per observation. Writing
#' \eqn{s = \theta + \mu} and \eqn{\psi} for the digamma function,
#' \deqn{\dfrac{\partial \ell}{\partial \mu} =
#'         \dfrac{\theta}{s}\left(\dfrac{y}{\mu} - 1\right), \qquad
#'       \dfrac{\partial \ell}{\partial \theta} = \psi(y+\theta) - \psi(\theta)
#'         + \log\dfrac{\theta}{s} + \dfrac{\mu - y}{s}.}
#' The mean component is the score of a generalized linear model with a
#' quadratic variance function.
#'
#' # The dispersion at large theta
#'
#' The dispersion component is **not** computed as written above. As
#' \eqn{\theta} grows the family tends to the Poisson and this derivative
#' vanishes, its three terms canceling to leading order:
#' \eqn{\psi(y+\theta) - \psi(\theta)} is \eqn{y/\theta},
#' \eqn{\log\{\theta/s\}} is \eqn{-\mu/\theta} and \eqn{(\mu-y)/s} is
#' \eqn{(\mu-y)/\theta}, and the three sum to zero, so the value is
#' \eqn{O(\theta^{-2})} computed from terms of size \eqn{\theta^{-1}}. Written
#' directly it is wrong by about one part in a thousand at \eqn{\theta = 10^6}
#' and **changes sign** at \eqn{10^8}.
#'
#' The kernel performs each cancellation symbolically instead. With
#' \eqn{a = \theta}, \eqn{b = \theta + y}, \eqn{c = \theta + \mu} and
#' \eqn{w = (y-\mu)/c},
#' \deqn{\dfrac{\partial \ell}{\partial \theta} =
#'       \left\{\psi(b) - \psi(a) - \log(1 + y/a)\right\}
#'       + \left\{\log(1 + w) - w\right\},}
#' each bracket evaluated by its own series where the arguments are large. To
#' leading order the result is \eqn{\{y - (y-\mu)^2\}/(2\theta^2)}, which the
#' example below reproduces. A fit reaches that regime routinely: 2,000 counts
#' drawn at a true \eqn{\theta} of 100 can report an estimate of order
#' \eqn{10^7}, and where such a run stops is decided by this arithmetic.
#'
#' With `scale = "link"` the generic applies the chain rule for the links the
#' family carries before returning. This method always returns the parameter
#' scale.
#'
#' @param distrib A `NegBin2Distrib` object, from [negbin2_distrib()].
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
#' \eqn{\operatorname{Var}(Y) = \mu + \mu^2/\theta}. \eqn{\psi} is the digamma
#' function, \eqn{\psi = (\log\Gamma)'}.
#'
#' @seealso [distrib_hessian.NegBin2Distrib()] for the second derivatives,
#'   [distrib_expected_hessian.NegBin2Distrib()] for their expectation and for
#'   what it does at large \eqn{\theta}, and [distrib_gradient()] for the
#'   generic.
#'
#' @examples
#' d <- negbin2_distrib()
#' y <- c(0, 2, 6)
#' th <- list(mu = 4, theta = 2)
#' g <- distrib_gradient(d, y, th)
#'
#' # The two closed forms, written out at a dispersion where the direct
#' # expression still has its digits.
#' s <- 2 + 4
#' all.equal(g$mu, (2 / s) * (y / 4 - 1))
#' all.equal(g$theta, digamma(y + 2) - digamma(2) + log(2 / s) + (4 - y) / s)
#'
#' # At large theta the direct form loses them and the kernel does not. The
#' # value is {y - (y - mu)^2}/(2 theta^2) to leading order; at y = 3, mu = 4
#' # that is 1/theta^2.
#' cmp <- function(t) {
#'   c(kernel = distrib_gradient(d, 3, list(mu = 4, theta = t))$theta,
#'     direct = digamma(3 + t) - digamma(t) + log(t / (t + 4)) + 1 / (t + 4),
#'     leading = 1 / t^2)
#' }
#' rbind(`1e+04` = cmp(1e4), `1e+06` = cmp(1e6), `1e+08` = cmp(1e8))
#'
#' # Summed over a fitted sample the score is at the optimizer's tolerance.
#' set.seed(5)
#' z <- distrib_rng(d, 2000, th)
#' fit <- fit_distrib(d, z)
#' vapply(distrib_gradient(d, z, as.list(coef(fit))), sum, numeric(1))
S7::method(distrib_gradient, NegBin2Distrib) <- function(distrib, y, theta, scale = c("parameter", "link"), ...,
                                       threads = 1L) {
  negbin_gradient_cpp(y, theta[[1]], theta[[2]], threads)
}

#' @title Negative Binomial Observed Hessian, NB2
#' @name distrib_hessian.NegBin2Distrib
#' @description
#' Computes the three distinct second derivatives of the negative binomial
#' log-mass with respect to \eqn{\mu} and \eqn{\theta}, one value per
#' observation. Writing \eqn{s = \theta + \mu} and \eqn{\psi_1} for the
#' trigamma function,
#' \deqn{\ell^{(\mu\mu)} = \dfrac{y+\theta}{s^2} - \dfrac{y}{\mu^2}, \qquad
#'       \ell^{(\mu\theta)} = \dfrac{y-\mu}{s^2}, \qquad
#'       \ell^{(\theta\theta)} = \psi_1(y+\theta) - \psi_1(\theta)
#'         + \dfrac{\mu}{\theta s} + \dfrac{y-\mu}{s^2}.}
#'
#' The dispersion entry cancels at large \eqn{\theta} for the same reason the
#' score does, and the kernel handles it the same way. Here the collapse is
#' better than a series: with \eqn{a = \theta}, \eqn{b = \theta + y} and
#' \eqn{c = \theta + \mu} the three leading terms combine **exactly**,
#' \deqn{-\dfrac{y}{ab} + \dfrac{\mu}{ac} + \dfrac{y-\mu}{c^2}
#'       = \dfrac{(y-\mu)^2}{b\,c^2},}
#' so what is computed is that quotient plus the remainder
#' \eqn{\psi_1(b) - \psi_1(a) + y/(ab)}, and only the remainder needs a series.
#' To leading order the result is \eqn{\{(y-\mu)^2 - y\}/\theta^3}, which is
#' the derivative of the score's leading term and is how the two derivations
#' check each other.
#'
#' @param distrib A `NegBin2Distrib` object, from [negbin2_distrib()].
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
#' @return A named list of three numeric vectors, `mu_mu`, `theta_theta` and
#'   `mu_theta`, in that order, each of length
#'   `max(length(y), length(mu), length(theta))`. The three name the distinct
#'   entries of a symmetric \eqn{2 \times 2} matrix per observation.
#'
#' @section Notation:
#' \eqn{\ell^{(ij)}} is the second derivative of the log-mass in parameters
#' \eqn{i} and \eqn{j}; parenthesized superscripts name derivatives. \eqn{\psi}
#' and \eqn{\psi_1} are the digamma and trigamma functions.
#'
#' @seealso [distrib_gradient.NegBin2Distrib()] for the score,
#'   [distrib_expected_hessian.NegBin2Distrib()] for the expectation of this
#'   quantity, [distrib_deriv3.NegBin2Distrib()] for the order above, and
#'   [distrib_hessian()] for the generic.
#'
#' @examples
#' d <- negbin2_distrib()
#' y <- c(0, 2, 6)
#' th <- list(mu = 4, theta = 2)
#' h <- distrib_hessian(d, y, th)
#' h
#'
#' # The three closed forms, written out at a dispersion where the direct
#' # expression still has its digits.
#' s <- 2 + 4
#' all.equal(h$mu_mu, (y + 2) / s^2 - y / 4^2)
#' all.equal(h$mu_theta, (y - 4) / s^2)
#' all.equal(h$theta_theta,
#'           trigamma(y + 2) - trigamma(2) + 4 / (2 * s) + (y - 4) / s^2)
#'
#' # The curvature in mu is positive at y = 0, so the observed information is
#' # not positive definite at every count.
#' h$mu_mu
#'
#' # It is the second derivative of the log-mass, so a central difference of
#' # the score reproduces it.
#' eps <- 1e-6
#' up <- distrib_gradient(d, y, list(mu = 4 + eps, theta = 2))$mu
#' dn <- distrib_gradient(d, y, list(mu = 4 - eps, theta = 2))$mu
#' all.equal((up - dn) / (2 * eps), h$mu_mu, tolerance = 1e-5)
S7::method(distrib_hessian, NegBin2Distrib) <- function(distrib, y, theta, scale = c("parameter", "link"), ...,
                                       threads = 1L) {
  negbin_hessian_cpp(y, theta[[1]], theta[[2]], threads)
}

#' @title Negative Binomial Expected Hessian, NB2
#' @name distrib_expected_hessian.NegBin2Distrib
#' @description
#' Returns the expectation of the observed Hessian under the model. With
#' \eqn{s = \theta + \mu} and \eqn{\psi_1} the trigamma function,
#' \deqn{\mathbb{E}\left[\ell^{(\mu\mu)}\right] = -\dfrac{\theta}{\mu s},
#'       \qquad
#'       \mathbb{E}\left[\ell^{(\mu\theta)}\right] = 0, \qquad
#'       \mathbb{E}\left[\ell^{(\theta\theta)}\right] =
#'         \mathbb{E}[\psi_1(Y+\theta)] - \psi_1(\theta) + \dfrac{\mu}{\theta s}.}
#' The mean entry and the mixed one are closed form, the latter vanishing
#' because the observed mixed entry is \eqn{(y-\mu)/s^2} and
#' \eqn{\mathbb{E}[Y] = \mu}. So the mean and the dispersion are orthogonal.
#'
#' \eqn{\mathbb{E}[\psi_1(Y+\theta)]} has no closed form. It is summed over the
#' support through the mass recurrence, carried in log scale until the terms
#' are representable, and stopped when the accumulated mass reaches
#' \eqn{1 - 10^{-12}}. The sum is exact to that mass, so this entry is a
#' truncated exact sum and not a quadrature or a simulation; `approx` and
#' `nsim` are ignored, and `y` is read only for its length.
#'
#' # A caveat at large theta
#'
#' The dispersion entry is a difference of two quantities that agree to
#' leading order as \eqn{\theta} grows, and **it is not rewritten to remove
#' that cancellation**, unlike the score and the observed Hessian. Its leading
#' order needs one term more of the observed Hessian than those two do, the
#' \eqn{\theta^{-3}} term vanishing under expectation, so it is a derivation of
#' its own. Measured at \eqn{\mu = 4}: it reads \eqn{-7.3\times 10^{-8}} at
#' \eqn{\theta = 10^2}, \eqn{-3.6\times 10^{-16}} at \eqn{10^4} and
#' \eqn{+1.7\times 10^{-16}} at \eqn{10^6}. The last is **positive**, which an
#' expected second derivative cannot be, so a Fisher scoring step taken there
#' is not reliable. The example below shows the sign turning.
#'
#' @param distrib A `NegBin2Distrib` object, from [negbin2_distrib()].
#' @param y A numeric vector of counts. Only its length is used.
#' @param theta A named list with components `mu` and `theta`, each a numeric
#'   vector of length 1 or of the length of `y`. A component of length 1 is
#'   recycled. Both must be strictly positive.
#' @param scale One of `"parameter"` (the default) or `"link"`, matched by
#'   [base::match.arg()]. Read by the generic, not by this method.
#' @param approx Ignored. The mean and mixed entries are closed form and the
#'   dispersion entry is a truncated exact sum. Accepted so that the signature
#'   matches the generic's.
#' @param nsim Ignored, for the same reason. Defaults to `10000`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#' @param threads A single positive integer, how many threads the kernel may
#'   use. Below the measured internal threshold the kernel stays sequential
#'   whatever the count says. Defaults to `1L`.
#'
#' @return A named list of three numeric vectors, `mu_mu`, `theta_theta` and
#'   `mu_theta`, in that order, each of length
#'   `max(length(y), length(mu), length(theta))` and constant within itself
#'   when the parameters are.
#'
#' @section Notation:
#' The **expected information** is
#' \eqn{\mathbb{E}[-\partial^2\ell/\partial\theta\,\partial\theta^\top]}, the
#' expectation of the **observed information** under the model. The negative
#' binomial is a regular family, so the second Bartlett identity holds and this
#' equals the variance of the score.
#'
#' @seealso [distrib_hessian.NegBin2Distrib()] for the observed quantity this
#'   is the expectation of, [fisher_scoring()], which inverts it at each step,
#'   and [distrib_expected_hessian()] for the generic.
#'
#' @examples
#' d <- negbin2_distrib()
#' th <- list(mu = 4, theta = 2)
#'
#' # The three entries, one value per observation.
#' lapply(distrib_expected_hessian(d, c(0, 2, 6), th), unique)
#'
#' # The mean entry is closed form and the mixed one is exactly zero.
#' c(closed = -2 / (4 * (2 + 4)),
#'   mixed = distrib_expected_hessian(d, 0, th)$mu_theta)
#'
#' # The observed Hessian averages onto them over a large sample.
#' set.seed(11)
#' z <- distrib_rng(d, 2e4, th)
#' rbind(observed = vapply(distrib_hessian(d, z, th), mean, numeric(1)),
#'       expected = vapply(distrib_expected_hessian(d, z, th),
#'                         function(v) v[1], numeric(1)))
#'
#' # The dispersion entry cancels at large theta and turns positive, which an
#' # expected second derivative cannot be.
#' vapply(c(1e2, 1e4, 1e6),
#'        function(t) distrib_expected_hessian(d, 0,
#'                      list(mu = 4, theta = t))$theta_theta,
#'        numeric(1))
S7::method(distrib_expected_hessian, NegBin2Distrib) <- function(distrib, y, theta, scale = c("parameter", "link"), approx = c("bartlett", "integrate", "mc", "opg"), nsim = 10000, ...,
                                       threads = 1L) {
  negbin_expected_hessian_cpp(y, theta[[1]], theta[[2]], threads)
}

#' @title Negative Binomial Third-Order Derivatives, NB2
#' @name distrib_deriv3.NegBin2Distrib
#' @description
#' Computes the four distinct third derivatives of the negative binomial
#' log-mass with respect to \eqn{\mu} and \eqn{\theta}, in closed form. Writing
#' \eqn{s = \theta + \mu}, the three components involving the mean are rational
#' in \eqn{(\mu, \theta)} and linear in \eqn{y}, and the pure dispersion
#' component carries \eqn{\psi_2(y+\theta) - \psi_2(\theta)}, with
#' \eqn{\psi_2} the second derivative of the digamma function.
#'
#' With `expected = TRUE` the expectations are returned. The three components
#' involving the mean need only \eqn{\mathbb{E}[Y] = \mu}; the pure dispersion
#' one needs \eqn{\mathbb{E}[\psi_2(Y+\theta)]}, which has no closed form and
#' is summed over the support with a far-tail correction, exactly as the
#' expected Hessian sums \eqn{\mathbb{E}[\psi_1(Y+\theta)]}. `approx` and
#' `nsim` are ignored either way.
#'
#' The cancellation of [distrib_gradient.NegBin2Distrib()] applies here too and
#' **is not removed at this order**: the pure dispersion components lose their
#' digits at large \eqn{\theta}, where the score itself is reliable to about
#' five figures and no further.
#'
#' @param distrib A `NegBin2Distrib` object, from [negbin2_distrib()].
#' @param y A numeric vector of counts. With `expected = TRUE` only its length
#'   is used.
#' @param theta A named list with components `mu` and `theta`, each a numeric
#'   vector of length 1 or of the length of `y`. A component of length 1 is
#'   recycled. Both must be strictly positive.
#' @param expected Logical of length 1. When `TRUE` the expectations under the
#'   model are returned in place of the observed values. Defaults to `FALSE`.
#' @param scale One of `"parameter"` (the default) or `"link"`, matched by
#'   [base::match.arg()]. Read by the generic, not by this method.
#' @param approx Ignored, the expected values being closed form or exact sums.
#' @param nsim Ignored, for the same reason. Defaults to `10000`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#' @param threads A single positive integer, how many threads the kernel may
#'   use. Below the measured internal threshold the kernel stays sequential
#'   whatever the count says. Defaults to `1L`.
#'
#' @return A named list of four numeric vectors, `mu_mu_mu`, `mu_mu_theta`,
#'   `mu_theta_theta` and `theta_theta_theta`, each of length
#'   `max(length(y), length(mu), length(theta))`. The names enumerate the
#'   distinct multi-indices of order three in two parameters.
#'
#' @section Notation:
#' \eqn{\ell^{(ijk)}} is the third derivative of the log-mass in parameters
#' \eqn{i}, \eqn{j} and \eqn{k}. \eqn{\psi} is the digamma function and
#' \eqn{\psi_m} its \eqn{m}th derivative.
#'
#' @seealso [distrib_hessian.NegBin2Distrib()] for the order below and
#'   [distrib_deriv4.NegBin2Distrib()] for the order above;
#'   [distrib_deriv3()] for the generic.
#'
#' @examples
#' d <- negbin2_distrib()
#' y <- c(0, 2, 6)
#' th <- list(mu = 4, theta = 2)
#' d3 <- distrib_deriv3(d, y, th)
#' names(d3)
#' d3$mu_mu_mu
#'
#' # The expected values are constants at a fixed parameter setting.
#' lapply(distrib_deriv3(d, y, th, expected = TRUE), unique)
#'
#' # A central difference of the Hessian reproduces the same component.
#' eps <- 1e-6
#' up <- distrib_hessian(d, y, list(mu = 4 + eps, theta = 2))$mu_mu
#' dn <- distrib_hessian(d, y, list(mu = 4 - eps, theta = 2))$mu_mu
#' all.equal((up - dn) / (2 * eps), d3$mu_mu_mu, tolerance = 1e-4)
S7::method(distrib_deriv3, NegBin2Distrib) <- function(distrib, y, theta, expected = FALSE, scale = c("parameter", "link"), approx = c("integrate", "bartlett", "mc", "opg"), nsim = 10000, ...,
                                       threads = 1L) {
  if (expected) negbin_deriv3_expected_cpp(y, theta[[1]], theta[[2]], threads)
  else negbin_deriv3_cpp(y, theta[[1]], theta[[2]], threads)
}

#' @title Negative Binomial Fourth-Order Derivatives, NB2
#' @name distrib_deriv4.NegBin2Distrib
#' @description
#' Computes the five distinct fourth derivatives of the negative binomial
#' log-mass with respect to \eqn{\mu} and \eqn{\theta}, in closed form, by the
#' same route the third order takes: the components involving the mean are
#' rational in \eqn{(\mu, \theta)} and linear in \eqn{y}, and the pure
#' dispersion component carries \eqn{\psi_3(y+\theta) - \psi_3(\theta)}, with
#' \eqn{\psi_3} the third derivative of the digamma function.
#'
#' With `expected = TRUE` the expectations are returned. The pure dispersion
#' one needs \eqn{\mathbb{E}[\psi_3(Y+\theta)]}, summed over the support with a
#' far-tail correction; the rest need only \eqn{\mathbb{E}[Y] = \mu}. `approx`
#' and `nsim` are ignored either way, and the caveat about large \eqn{\theta}
#' on [distrib_deriv3.NegBin2Distrib()] applies here as well.
#'
#' @param distrib A `NegBin2Distrib` object, from [negbin2_distrib()].
#' @param y A numeric vector of counts. With `expected = TRUE` only its length
#'   is used.
#' @param theta A named list with components `mu` and `theta`, each a numeric
#'   vector of length 1 or of the length of `y`. A component of length 1 is
#'   recycled. Both must be strictly positive.
#' @param expected Logical of length 1. When `TRUE` the expectations under the
#'   model are returned in place of the observed values. Defaults to `FALSE`.
#' @param scale One of `"parameter"` (the default) or `"link"`, matched by
#'   [base::match.arg()]. Read by the generic, not by this method.
#' @param approx Ignored, the expected values being closed form or exact sums.
#' @param nsim Ignored, for the same reason. Defaults to `10000`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#' @param threads A single positive integer, how many threads the kernel may
#'   use. Below the measured internal threshold the kernel stays sequential
#'   whatever the count says. Defaults to `1L`.
#'
#' @return A named list of five numeric vectors, `mu_mu_mu_mu`,
#'   `mu_mu_mu_theta`, `mu_mu_theta_theta`, `mu_theta_theta_theta` and
#'   `theta_theta_theta_theta`, each of length
#'   `max(length(y), length(mu), length(theta))`.
#'
#' @section Notation:
#' \eqn{\ell^{(ijkl)}} is the fourth derivative of the log-mass in parameters
#' \eqn{i}, \eqn{j}, \eqn{k} and \eqn{l}. \eqn{\psi} is the digamma function
#' and \eqn{\psi_m} its \eqn{m}th derivative.
#'
#' @seealso [distrib_deriv3.NegBin2Distrib()] for the order below,
#'   [distrib_deriv4()] for the generic, and
#'   [distrib_expected_hessian.NegBin2Distrib()] for the second-order
#'   expectation these extend.
#'
#' @examples
#' d <- negbin2_distrib()
#' y <- c(0, 2, 6)
#' th <- list(mu = 4, theta = 2)
#' d4 <- distrib_deriv4(d, y, th)
#' names(d4)
#'
#' # The expected values are constants at a fixed parameter setting.
#' lapply(distrib_deriv4(d, y, th, expected = TRUE), unique)
#'
#' # A central difference of the third order reproduces it.
#' eps <- 1e-6
#' up <- distrib_deriv3(d, y, list(mu = 4 + eps, theta = 2))$mu_mu_mu
#' dn <- distrib_deriv3(d, y, list(mu = 4 - eps, theta = 2))$mu_mu_mu
#' all.equal((up - dn) / (2 * eps), d4$mu_mu_mu_mu, tolerance = 1e-3)
S7::method(distrib_deriv4, NegBin2Distrib) <- function(distrib, y, theta, expected = FALSE, scale = c("parameter", "link"), approx = c("integrate", "bartlett", "mc", "opg"), nsim = 10000, ...,
                                       threads = 1L) {
  if (expected) negbin_deriv4_expected_cpp(y, theta[[1]], theta[[2]], threads)
  else negbin_deriv4_cpp(y, theta[[1]], theta[[2]], threads)
}

# --- CONSTRUCTOR WRAPPER ---

#' Negative Binomial Distribution, NB2
#'
#' @description
#' Builds the distribution object for the negative binomial family on the
#' non-negative integers in the NB2 parametrization: the mean \eqn{\mu > 0} and
#' a dispersion \eqn{\theta > 0}, so that
#' \eqn{\operatorname{Var}(Y) = \mu + \mu^2/\theta}. This is the count model a
#' regression reaches for when a Poisson is overdispersed. The returned object
#' carries closed-form derivatives of the log-mass to fourth order and
#' closed-form moments.
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
#' # The parametrization
#'
#' The mass on \eqn{y = 0, 1, 2, \ldots} is
#' \deqn{P(Y = y; \mu, \theta) = \dfrac{\Gamma(y+\theta)}{y!\,\Gamma(\theta)}
#'       \left(\dfrac{\theta}{s}\right)^{\theta}
#'       \left(\dfrac{\mu}{s}\right)^{y}, \qquad s = \theta + \mu,}
#' the distribution function is the partial sum and the quantile function its
#' generalized inverse. The mean is \eqn{\mu}, the variance
#' \eqn{\mu + \mu^2/\theta}, the skewness
#' \eqn{(\theta + 2\mu)/\sqrt{\mu\theta(\theta+\mu)}} and the excess kurtosis
#' \eqn{6/\theta + \theta/\{\mu(\theta+\mu)\}}.
#'
#' The law is a Poisson whose rate is gamma distributed with mean \eqn{\mu} and
#' shape \eqn{\theta}. Two limits follow: \eqn{\theta = 1} gives the geometric,
#' [geometric_distrib()], and \eqn{\theta \to \infty} gives the Poisson,
#' [poisson_distrib()]. **The dispersion reads the opposite way round from a
#' variance**: a small \eqn{\theta} is heavy overdispersion and a large one is
#' nearly a Poisson.
#'
#' The numbering follows Cameron and Trivedi: NB2 has a variance quadratic in
#' the mean, and [negbin1_distrib()] has one linear in it.
#'
#' # Derivatives
#'
#' With \eqn{\psi} the digamma function, the score is
#' \deqn{\dfrac{\partial \ell}{\partial \mu} =
#'         \dfrac{\theta}{s}\left(\dfrac{y}{\mu} - 1\right), \qquad
#'       \dfrac{\partial \ell}{\partial \theta} = \psi(y+\theta) - \psi(\theta)
#'         + \log\dfrac{\theta}{s} + \dfrac{\mu - y}{s},}
#' and the expected Hessian is
#' \deqn{\mathbb{E}\left[\ell^{(\mu\mu)}\right] = -\dfrac{\theta}{\mu s},
#'       \qquad
#'       \mathbb{E}\left[\ell^{(\mu\theta)}\right] = 0, \qquad
#'       \mathbb{E}\left[\ell^{(\theta\theta)}\right] =
#'         \mathbb{E}[\psi_1(Y+\theta)] - \psi_1(\theta) + \dfrac{\mu}{\theta s}.}
#' The zero off-diagonal makes the mean and the dispersion orthogonal.
#'
#' # Every derivative in theta cancels as theta grows
#'
#' The family tends to the Poisson as \eqn{\theta} grows, so every derivative
#' in \eqn{\theta} vanishes there and is written as a sum of terms that cancel
#' to leading order. The score and the observed Hessian are computed in forms
#' that perform those cancellations symbolically, and are reliable throughout;
#' the direct expression for the score is wrong by one part in a thousand at
#' \eqn{\theta = 10^6} and changes sign at \eqn{10^8}. Two quantities are
#' **not** rewritten, and both pages say so:
#'
#' - the expected information in \eqn{\theta}, which is measured positive from
#'   about \eqn{\theta = 10^6} and so is unusable there
#'   ([distrib_expected_hessian.NegBin2Distrib()]);
#' - the third and fourth derivatives in \eqn{\theta}
#'   ([distrib_deriv3.NegBin2Distrib()], [distrib_deriv4.NegBin2Distrib()]).
#'
#' A fit reaches that regime routinely, a nearly equidispersed sample driving
#' \eqn{\theta} towards its boundary.
#'
#' # Estimation
#'
#' [fit_distrib()] maximizes the log-likelihood on the link scale. The mean has
#' the closed-form estimate \eqn{\hat\mu = \bar y}; the dispersion solves an
#' equation in digamma functions and is reached numerically. The method of
#' moments supplies the starting value \eqn{\bar y^2/(s^2 - \bar y)}, with
#' \eqn{s^2} the sample variance.
#'
#' @section Notation:
#' \eqn{\ell} is the log-mass of one observation, \eqn{\mu > 0} the mean and
#' \eqn{\theta > 0} the dispersion. \eqn{\psi} is the digamma function and
#' \eqn{\psi_m} its \eqn{m}th derivative. \eqn{\eta} is a parameter on the
#' unconstrained scale of its link, with \eqn{\theta_j = g^{-1}(\eta_j)}.
#'
#' @return An S7 object of class `NegBin2Distrib`, inheriting from
#'   `discrete_distrib`, with `distrib_name` `"negbin2"`, `dimension`
#'   `"univariate"`, `bounds` `c(0, Inf)`, `params` `c("mu", "theta")`,
#'   `n_params` `2`, `params_bounds` the domain \eqn{(0, \infty)} for both, and
#'   `link_params` the two links given here.
#'
#' @seealso
#' [negbin1_distrib()] for the NB1 parametrization, with a variance linear in
#' the mean; [poisson_distrib()] for the equidispersed limit and
#' [geometric_distrib()] for \eqn{\theta = 1}; [pig1_distrib()] for the
#' Poisson mixed over an inverse Gaussian instead of a gamma;
#' [zero_inflated()] and [zero_adjusted()] for counts with excess zeros;
#' [fit_distrib()] to estimate the parameters; [check_distrib()] to validate a
#' family of your own against the same battery this one passes;
#' [NegBin2Distrib] for the class.
#'
#' @references
#' Cameron, A. C. and Trivedi, P. K. (2013). *Regression Analysis of Count
#' Data*, 2nd edition, Chapter 3. Cambridge University Press, Cambridge.
#'
#' @importFrom linkfunctions7 log_link
#' @importFrom stats dnbinom pnbinom qnbinom rnbinom
#'
#' @examples
#' d <- negbin2_distrib()
#' d
#'
#' # The mass is stats::dnbinom at size = theta, mu = mu.
#' y <- c(0, 2, 6)
#' th <- list(mu = 4, theta = 2)
#' all.equal(distrib_pdf(d, y, th), dnbinom(y, size = 2, mu = 4))
#'
#' # Moments: the variance exceeds the mean, by mu^2/theta.
#' c(mean = mean(d, th), var = variance(d, th),
#'   skew = skewness(d, th), kurt = kurtosis(d, th))
#' c(4 + 4^2 / 2, 6 / 2 + 2 / (4 * 6))
#'
#' # theta = 1 is the geometric and a large theta is the Poisson.
#' all.equal(distrib_pdf(d, y, list(mu = 4, theta = 1)),
#'           distrib_pdf(geometric_distrib(), y, list(mu = 4)))
#' rbind(negbin = distrib_pdf(d, 0:4, list(mu = 3, theta = 1e6)),
#'       poisson = dpois(0:4, 3))
#'
#' # Fitting recovers the parameters; the moment estimates start it off.
#' set.seed(5)
#' z <- distrib_rng(d, 2000, th)
#' fit <- fit_distrib(d, z)
#' rbind(fitted  = coef(fit),
#'       moments = c(mu = mean(z),
#'                   theta = mean(z)^2 / (var(z) - mean(z))))
#'
#' @export
negbin2_distrib <- function(link_mu = log_link(), link_theta = log_link()) {

  NegBin2Distrib(
    distrib_name = "negbin2",
    dimension = "univariate",
    bounds = c(0, Inf),

    params = c("mu", "theta"),
    params_interpretation = c(mu = "mean", theta = "dispersion"),
    n_params = 2,

    params_bounds = list(
      mu = c(0, Inf),
      theta = c(0, Inf)
    ),

    link_params = list(
      mu = link_mu,
      theta = link_theta
    )
  )

}
