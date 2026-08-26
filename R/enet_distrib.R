#' @include distrib.R generics.R cross_derivatives.R
NULL

# The two quantities every derivative of this family is written in.
#
#   x = a / sqrt(c),   M(x) = Phi(-x) / phi(x),   G = d log M / dx
#
# M is the Mills ratio proper and G = x - 1/M, where 1/M is what
# numericals7::mills_ratio() computes stably at -x. The difference cancels
# for large x -- both terms are of size x while G is of size 1/x -- so
# above a threshold the asymptotic series is used instead; the crossover
# is where the two agree to about 1e-11.
.enet_G <- function(x) {
  out <- numeric(length(x))
  big <- abs(x) > 1e3
  if (any(!big)) {
    xb <- x[!big]
    out[!big] <- xb - numericals7::mills_ratio(-xb)$r
  }
  if (any(big)) {
    xb <- x[big]
    u <- 1 / xb
    out[big] <- -u + 2 * u^3 - 10 * u^5 + 74 * u^7
  }
  out
}

# log M(x). The direct form adds x^2/2 to a log-probability of the same
# size and opposite sign, so it loses a digit for every factor of ten in
# x; past 30 the asymptotic series of the Mills ratio is used, where the
# first neglected term is below 1e-12 relative.
.enet_logM <- function(x) {
  out <- numeric(length(x))
  big <- x > 30
  if (any(!big)) {
    xb <- x[!big]
    out[!big] <- stats::pnorm(-xb, log.p = TRUE) + xb^2 / 2 +
      0.5 * base::log(2 * pi)
  }
  if (any(big)) {
    u <- 1 / x[big]
    u2 <- u^2
    out[big] <- -base::log(x[big]) +
      base::log1p(u2 * (-1 + u2 * (3 + u2 * (-15 + u2 * 105))))
  }
  out
}

# log Phi(-x), which underflows past about 38 on the natural scale
.enet_logQ <- function(x) stats::pnorm(-x, log.p = TRUE)

# the pieces the derivatives share: the rates, x, G and G'
.enet_parts <- function(theta) {
  lam <- theta[[2]]
  al <- theta[[3]]
  a <- lam * al
  cc <- lam * (1 - al)
  x <- a / sqrt(cc)
  g <- .enet_G(x)
  list(mu = theta[[1]], lam = lam, al = al, a = a, c = cc, x = x,
       g = g, dg = 1 + x * g - g^2)
}

#' @title Elastic-Net Distribution Class
#' @name EnetDistrib
#'
#' @description
#' The S7 class of the density whose negative logarithm is the elastic-net
#' penalty: a Laplace and a Gaussian at the same location, multiplied together
#' and normalized,
#' \deqn{f(y) \propto \exp\{-\lambda\alpha|y-\mu|
#'       - \lambda(1-\alpha)(y-\mu)^2/2\}.}
#' At \eqn{\alpha \to 1} it is the Laplace of [Laplace2Distrib] and at
#' \eqn{\alpha \to 0} the Gaussian; \eqn{\alpha} is confined to the open
#' interval, as every bounded parameter in this package is, so neither end is
#' a member. Both remain families of their own.
#'
#' Like the Laplace, its log-likelihood is **not differentiable in
#' \eqn{\mu}**: the absolute value carries a kink at the location, and
#' `params_smooth` records `mu = FALSE`.
#'
#' Build one with [enet_distrib()], which supplies the three link functions.
#' This page documents the raw S7 constructor, which validates none of the
#' relationships between its properties.
#'
#' @inheritParams distrib
#'
#' @return An S7 object of class `EnetDistrib`, inheriting from
#'   `continuous_distrib` and from `distrib`. For an object built by
#'   [enet_distrib()] the properties hold `"enet"`, `"univariate"`,
#'   `c(-Inf, Inf)`, `c("mu", "lambda", "alpha")`, the interpretations
#'   `c(mu = "location", lambda = "rate", alpha = "mixing weight")`, `3`, the
#'   domains \eqn{(-\infty,\infty)}, \eqn{(0,\infty)} and \eqn{(0,1)}, and
#'   `params_smooth = c(mu = FALSE, lambda = TRUE, alpha = TRUE)`.
#'
#' @section Methods:
#' Registered in this file:
#'   [`distrib_pdf()`][distrib_pdf.EnetDistrib],
#'   [`distrib_cdf()`][distrib_cdf.EnetDistrib],
#'   [`distrib_quantile()`][distrib_quantile.EnetDistrib],
#'   [`distrib_rng()`][distrib_rng.EnetDistrib],
#'   [`distrib_gradient()`][distrib_gradient.EnetDistrib],
#'   [`distrib_hessian()`][distrib_hessian.EnetDistrib],
#'   [`distrib_expected_hessian()`][distrib_expected_hessian.EnetDistrib],
#'   [`distrib_grad_y()`][distrib_grad_y.EnetDistrib],
#'   [`distrib_hess_y()`][distrib_hess_y.EnetDistrib],
#'   [`distrib_cross_y()`][distrib_cross_y.EnetDistrib],
#'   [`mean()`][mean.EnetDistrib], [`variance()`][variance.EnetDistrib],
#'   [`skewness()`][skewness.EnetDistrib].
#'
#' Registered in `enet_higher.R`, both closed form:
#'   [`distrib_deriv3()`][distrib_deriv3.EnetDistrib],
#'   [`distrib_deriv4()`][distrib_deriv4.EnetDistrib].
#'
#' The **kurtosis** has no closed form here and comes from the base class by
#' quadrature. So does the fourth response derivative and beyond, which are
#' zero: the log-density is quadratic in \eqn{y} away from the location.
#'
#' @section The two quantities everything is written in:
#' With \eqn{a = \lambda\alpha}, \eqn{c = \lambda(1-\alpha)} and
#' \eqn{x = a/\sqrt c}, every derivative in the two rates is a polynomial in
#' \eqn{x} and \eqn{G = \mathrm{d}\log M/\mathrm{d}x}, where \eqn{M} is the
#' Mills ratio. The map \eqn{(\lambda, \alpha) \mapsto (a, c)} is bilinear, so
#' the chain onto the reported parameters adds one cross term and nothing else.
#'
#' @seealso [enet_distrib()] to build one;
#'   [laplace2_distrib()] and [gaussian1_distrib()] for the two ends;
#'   [penalties7::elasticnet_penalty()], the consumer this family exists for.
#'
#' @examples
#' d <- enet_distrib()
#' S7::S7_inherits(d, continuous_distrib)
#'
#' d@params
#' d@params_bounds
#'
#' # The location carries a kink, so it is declared non-smooth. That is what
#' # switches off the finite-difference guard in check_distrib().
#' d@params_smooth
#'
#' # The two ends, approached but not reached.
#' th <- list(mu = 0, lambda = 2, alpha = 0.5)
#' c(laplace = distrib_pdf(d, 0.7, list(mu = 0, lambda = 2, alpha = 1 - 1e-10)),
#'   laplace_exact = distrib_pdf(laplace2_distrib(), 0.7,
#'                               list(mu = 0, lambda = 2)),
#'   gaussian = distrib_pdf(d, 0.7, list(mu = 0, lambda = 2, alpha = 1e-10)),
#'   gaussian_exact = dnorm(0.7, 0, 1 / sqrt(2)))
EnetDistrib <- S7::new_class("EnetDistrib", parent = continuous_distrib)

# --- S7 METHODS IMPLEMENTATION ---

#' @title Elastic-Net Density
#' @name distrib_pdf.EnetDistrib
#'
#' @description
#' Computes the elastic-net density
#' \deqn{f(y; \mu, \lambda, \alpha) = \frac{1}{Z}
#'   \exp\left\{-a|y-\mu| - \tfrac{c}{2}(y-\mu)^2\right\},}
#' with \eqn{a = \lambda\alpha}, \eqn{c = \lambda(1-\alpha)} and
#' \eqn{Z = 2M(a/\sqrt{c})/\sqrt{c}}, where \eqn{M} is the Mills ratio. The
#' constant is finite at both ends of the mixing weight: \eqn{2/a} as
#' \eqn{\alpha \to 1} and \eqn{\sqrt{2\pi/c}} as \eqn{\alpha \to 0}.
#'
#' The constant goes through \eqn{\log M}. Written directly, \eqn{\log M(x)}
#' adds \eqn{x^2/2} to a log-probability of the same size and opposite sign,
#' so it loses a digit for every factor of ten in \eqn{x}; past \eqn{x = 30}
#' the asymptotic series is used instead. At \eqn{\alpha = 1 - 10^{-12}} the
#' argument reaches \eqn{10^6}, which an ordinary elastic net does.
#'
#' @param distrib An `EnetDistrib` object, from [enet_distrib()].
#' @param y A numeric vector of observations, anywhere on the real line.
#' @param theta A named list with components `mu`, `lambda` and `alpha`, each a
#'   numeric vector of length 1 or of the length of `y`. `lambda` must be
#'   strictly positive and `alpha` strictly inside \eqn{(0, 1)}; at either
#'   endpoint the constant's argument is not defined and the value is `NaN`.
#' @param log Logical of length 1. When `TRUE` the log-density is returned.
#'   Defaults to `FALSE`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of densities, of the length of the recycled inputs.
#'
#' @section Notation:
#' \eqn{\mu} is the location, \eqn{\lambda > 0} the overall rate,
#' \eqn{\alpha \in (0,1)} the mixing weight, \eqn{M} the Mills ratio
#' \eqn{\Phi(-x)/\phi(x)}, and \eqn{Z} the normalizing constant.
#'
#' @seealso [distrib_cdf.EnetDistrib()] for the distribution function,
#'   [distrib_gradient.EnetDistrib()] for the score,
#'   [laplace2_distrib()] and [gaussian1_distrib()] for the two ends, and
#'   [distrib_pdf()] for the generic.
#'
#' @examples
#' d <- enet_distrib()
#' y <- c(-1.5, -0.3, 0.4, 2.1)
#' th <- list(mu = 0, lambda = 2, alpha = 0.5)
#'
#' # It integrates to one.
#' integrate(function(u) distrib_pdf(d, u, th), -Inf, Inf)$value
#'
#' # The kernel written out, up to the constant.
#' r <- distrib_pdf(d, y, th) /
#'      exp(-1 * abs(y) - 1 * y^2 / 2)     # a = c = 1 at alpha = 0.5
#' all.equal(r, rep(r[1], length(r)))
#'
#' # The two ends, approached to ten figures.
#' rbind(alpha = c(1 - 1e-10, 1e-10),
#'       enet = c(distrib_pdf(d, 0.7, list(mu = 0, lambda = 2,
#'                                         alpha = 1 - 1e-10)),
#'                distrib_pdf(d, 0.7, list(mu = 0, lambda = 2,
#'                                         alpha = 1e-10))),
#'       limit = c(distrib_pdf(laplace2_distrib(), 0.7,
#'                             list(mu = 0, lambda = 2)),
#'                 dnorm(0.7, 0, 1 / sqrt(2))))
S7::method(distrib_pdf, EnetDistrib) <- function(distrib, y, theta,
                                                 log = FALSE, ...) {
  p <- .enet_parts(theta)
  z <- y - p$mu
  log_z <- base::log(2) - 0.5 * base::log(p$c) + .enet_logM(p$x)
  log_d <- -p$a * abs(z) - p$c * z^2 / 2 - log_z
  if (log) log_d else exp(log_d)
}

#' @title Elastic-Net Distribution Function
#' @name distrib_cdf.EnetDistrib
#'
#' @description
#' Computes the distribution function in closed form. Each half of the density
#' is a truncated Gaussian, so with \eqn{z = q-\mu} and \eqn{x = a/\sqrt c},
#' \deqn{F(q) = \dfrac{\Phi(\sqrt{c}\,z - x)}{2\Phi(-x)} \quad (z \le 0),
#'       \qquad
#'       F(q) = 1 - \dfrac{\Phi(-\sqrt{c}\,z - x)}{2\Phi(-x)} \quad (z > 0).}
#'
#' Both ratios are taken on the log scale and exponentiated. The denominator
#' \eqn{\Phi(-x)} is exactly zero past \eqn{x = 38} in double precision, and
#' \eqn{x} reaches that at ordinary values of \eqn{\alpha}: at
#' \eqn{\lambda = 20}, \eqn{\alpha = 0.995} it is 63.
#'
#' @param distrib An `EnetDistrib` object, from [enet_distrib()].
#' @param q A numeric vector of quantiles.
#' @param theta A named list with components `mu`, `lambda` and `alpha`, each a
#'   numeric vector of length 1 or of the length of `q`. `lambda` must be
#'   strictly positive and `alpha` strictly inside \eqn{(0, 1)}.
#' @param lower.tail Logical of length 1. When `TRUE`, the default, the value
#'   is \eqn{P(Y \le q)}; when `FALSE` it is one minus that.
#' @param log.p Logical of length 1. When `TRUE` the logarithm of the
#'   probability is returned. Defaults to `FALSE`.
#'
#' @return A numeric vector of probabilities in \eqn{[0, 1]}, or their
#'   logarithms with `log.p = TRUE`, of the length of the recycled inputs.
#'
#' @section Notation:
#' \eqn{a = \lambda\alpha}, \eqn{c = \lambda(1-\alpha)}, \eqn{x = a/\sqrt c},
#' and \eqn{\Phi} the standard Gaussian distribution function.
#'
#' @seealso [distrib_quantile.EnetDistrib()], which inverts this in closed
#'   form, [distrib_pdf.EnetDistrib()] for the density, and [distrib_cdf()]
#'   for the generic.
#'
#' @examples
#' d <- enet_distrib()
#' q <- c(-2, -0.5, 0, 0.5, 2)
#' th <- list(mu = 0, lambda = 2, alpha = 0.5)
#'
#' # Against a direct quadrature of the density.
#' rbind(closed = distrib_cdf(d, q, th),
#'       quadrature = vapply(q, function(u)
#'         integrate(function(v) distrib_pdf(d, v, th), -Inf, u)$value, 0))
#'
#' # The density is symmetric about the location, so the median is mu.
#' distrib_cdf(d, 0, th)
#'
#' # An ordinary elastic net puts the Mills argument past 38, where the
#' # denominator of both halves has underflowed to zero on the natural scale.
#' hard <- list(mu = 0, lambda = 20, alpha = 0.995)
#' c(x = 20 * 0.995 / sqrt(20 * 0.005),
#'   phi_minus_x = pnorm(-20 * 0.995 / sqrt(20 * 0.005)),
#'   cdf = distrib_cdf(d, 0.1, hard))
S7::method(distrib_cdf, EnetDistrib) <- function(distrib, q, theta,
                                                 lower.tail = TRUE,
                                                 log.p = FALSE) {
  p <- .enet_parts(theta)
  z <- q - p$mu
  s <- sqrt(p$c)
  # each half is a ratio of normal tails, taken on the log scale: the
  # denominator Phi(-x) underflows past x = 38 and the mixing weight
  # reaches that at ordinary values of alpha
  lq <- .enet_logQ(p$x)
  half_lo <- exp(stats::pnorm(s * z - p$x, log.p = TRUE) - lq) / 2
  half_hi <- exp(stats::pnorm(-s * z - p$x, log.p = TRUE) - lq) / 2
  lower <- ifelse(z <= 0, half_lo, 1 - half_hi)
  out <- if (lower.tail) lower else 1 - lower
  if (log.p) base::log(out) else out
}

#' @title Elastic-Net Quantile Function
#' @name distrib_quantile.EnetDistrib
#'
#' @description
#' Inverts [distrib_cdf.EnetDistrib()] in closed form, each half of the density
#' being a truncated Gaussian. Below the median the quantile is
#' \eqn{\mu + \{x + \Phi^{-1}(2p\,\Phi(-x))\}/\sqrt c} and above it the
#' reflection. Nothing is inverted by root finding, so
#' [distrib_rng.EnetDistrib()] can use the inverse transform.
#'
#' The `qnorm` call is made on the **log** scale, its argument being
#' \eqn{\log(2p) + \log\Phi(-x)}, because \eqn{\Phi(-x)} underflows to zero
#' past \eqn{x = 38}.
#'
#' @param distrib An `EnetDistrib` object, from [enet_distrib()].
#' @param p A numeric vector of probabilities in \eqn{[0, 1]}, or their
#'   logarithms when `log.p = TRUE`.
#' @param theta A named list with components `mu`, `lambda` and `alpha`, each a
#'   numeric vector of length 1 or of the length of `p`.
#' @param lower.tail Logical of length 1. When `TRUE`, the default, `p` is
#'   \eqn{P(Y \le q)}; when `FALSE` it is the survival probability.
#' @param log.p Logical of length 1. When `TRUE`, `p` is given as a logarithm
#'   and is exponentiated first. Defaults to `FALSE`.
#'
#' @return A numeric vector of quantiles, of the length of the recycled
#'   inputs.
#'
#' @seealso [distrib_cdf.EnetDistrib()], which it inverts,
#'   [distrib_rng.EnetDistrib()], which draws from it, and
#'   [distrib_quantile()] for the generic.
#'
#' @examples
#' d <- enet_distrib()
#' th <- list(mu = 0, lambda = 2, alpha = 0.5)
#'
#' # The round trip is exact, both functions being closed form.
#' q <- c(-2, -0.5, 0.5, 2)
#' all.equal(distrib_quantile(d, distrib_cdf(d, q, th), th), q)
#'
#' # Symmetric about the location.
#' p <- c(0.1, 0.25, 0.5, 0.75, 0.9)
#' round(distrib_quantile(d, p, th), 10)
#'
#' # It still answers where the Mills argument has passed 38.
#' distrib_quantile(d, c(0.25, 0.75),
#'                  list(mu = 0, lambda = 20, alpha = 0.995))
S7::method(distrib_quantile, EnetDistrib) <- function(distrib, p, theta,
                                                      lower.tail = TRUE,
                                                      log.p = FALSE) {
  pr <- if (log.p) exp(p) else p
  if (!lower.tail) pr <- 1 - pr
  pa <- .enet_parts(theta)
  s <- sqrt(pa$c)
  lq <- .enet_logQ(pa$x)
  lo <- (pa$x + stats::qnorm(base::log(2 * pr) + lq, log.p = TRUE)) / s
  hi <- -(pa$x + stats::qnorm(base::log(2 * (1 - pr)) + lq, log.p = TRUE)) / s
  pa$mu + ifelse(pr <= 0.5, lo, hi)
}

#' @title Elastic-Net Random Generation
#' @name distrib_rng.EnetDistrib
#'
#' @description
#' Draws by the inverse transform, `distrib_quantile(distrib, runif(n),
#' theta)`. The quantile function of this family is closed form, so the
#' transform is exact and costs one uniform per draw; the base class's
#' ratio-of-uniforms fallback is not needed.
#'
#' @param distrib An `EnetDistrib` object, from [enet_distrib()].
#' @param n A single positive integer, the number of draws.
#' @param theta A named list with components `mu`, `lambda` and `alpha`, each a
#'   numeric vector of length 1 or of length `n`; a component of length 1 is
#'   recycled.
#'
#' @return A numeric vector of `n` draws, symmetric about `mu`.
#'
#' @seealso [distrib_quantile.EnetDistrib()], which it inverts through,
#'   [variance.EnetDistrib()] for the moment the draws reproduce, and
#'   [distrib_rng()] for the generic.
#'
#' @examples
#' d <- enet_distrib()
#' th <- list(mu = 0, lambda = 2, alpha = 0.5)
#'
#' set.seed(71)
#' x <- distrib_rng(d, 2e5, th)
#'
#' # Two moments against their closed forms, and the symmetry.
#' rbind(sample = c(mean(x), var(x), mean((x - mean(x))^3) / sd(x)^3),
#'       theory = c(mean(d, th), variance(d, th), skewness(d, th)))
#'
#' # As alpha approaches one the draws become Laplace.
#' set.seed(72)
#' xl <- distrib_rng(d, 1e5, list(mu = 0, lambda = 2, alpha = 1 - 1e-10))
#' c(sample_var = var(xl), laplace_var = 2 / 2^2)
S7::method(distrib_rng, EnetDistrib) <- function(distrib, n, theta) {
  distrib_quantile(distrib, stats::runif(n), theta)
}

# the derivatives of log Z with respect to the two rates
.enet_logz_derivs <- function(p) {
  g <- p$g
  dg <- p$dg
  x <- p$x
  cc <- p$c
  za <- g / sqrt(cc)
  zc <- -(1 + x * g) / (2 * cc)
  zaa <- dg / cc
  zac <- -(x * dg + g) / (2 * cc^1.5)
  zcc <- (2 * (1 + x * g) + x * (g + x * dg)) / (4 * cc^2)
  list(za = za, zc = zc, zaa = zaa, zac = zac, zcc = zcc)
}

#' @title Elastic-Net Score
#' @name distrib_gradient.EnetDistrib
#'
#' @description
#' Computes the three first derivatives of the log-density in closed form. With
#' \eqn{z = y-\mu}, \eqn{a = \lambda\alpha} and \eqn{c = \lambda(1-\alpha)},
#' \deqn{\dfrac{\partial\ell}{\partial\mu} = a\,\mathrm{sgn}(z) + cz,}
#' and the two rate components are the data terms less the derivatives of
#' \eqn{\log Z}, which are \eqn{G/\sqrt c} in \eqn{a} and \eqn{-(1+xG)/(2c)} in
#' \eqn{c}. The chain to \eqn{(\lambda, \alpha)} is linear, the map being
#' bilinear.
#'
#' The location component is **undefined at \eqn{y = \mu}**, where
#' \eqn{\mathrm{sgn}} jumps; `sign(0)` is 0 in R, so the value returned there
#' is \eqn{cz = 0}, one point of the subdifferential. `params_smooth`
#' records `mu = FALSE` so that
#' [check_distrib()]'s finite-difference guard knows.
#'
#' @param distrib An `EnetDistrib` object, from [enet_distrib()].
#' @param y A numeric vector of observations.
#' @param theta A named list with components `mu`, `lambda` and `alpha`, each a
#'   numeric vector of length 1 or of the length of `y`.
#' @param scale Either `"parameter"`, the default, or `"link"`. The
#'   transformation to the link scale is applied in the generic's body, so this
#'   method always returns the parameter scale.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A named list of three numeric vectors, `mu`, `lambda` and `alpha`,
#'   each of the length of the recycled inputs.
#'
#' @section Notation:
#' \eqn{a = \lambda\alpha}, \eqn{c = \lambda(1-\alpha)}, \eqn{x = a/\sqrt c},
#' \eqn{G = \mathrm{d}\log M/\mathrm{d}x} with \eqn{M} the Mills ratio, and
#' \eqn{Z} the normalizing constant.
#'
#' @seealso [distrib_hessian.EnetDistrib()] for the second derivatives,
#'   [distrib_expected_hessian.EnetDistrib()] for the information, which does
#'   not agree with the observed curvature at the kink, and
#'   [distrib_gradient()] for the generic.
#'
#' @examples
#' d <- enet_distrib()
#' y <- c(-1.5, -0.3, 0.4, 2.1)
#' th <- list(mu = 0, lambda = 2, alpha = 0.5)
#' g <- distrib_gradient(d, y, th)
#'
#' # The location component written out.
#' all.equal(g$mu, 1 * sign(y) + 1 * y)
#'
#' # The two rate components against numerical differentiation.
#' f <- function(p) sum(distrib_pdf(d, y, list(mu = 0, lambda = p[1],
#'                                             alpha = p[2]), log = TRUE))
#' rbind(analytic = c(sum(g$lambda), sum(g$alpha)),
#'       numeric = numDeriv::grad(f, c(2, 0.5)))
#'
#' # At the location the derivative does not exist; the value returned is
#' # one point of the subdifferential.
#' c(at_kink = distrib_gradient(d, 0, th)$mu,
#'   just_below = distrib_gradient(d, -1e-9, th)$mu,
#'   just_above = distrib_gradient(d, 1e-9, th)$mu)
S7::method(distrib_gradient, EnetDistrib) <- function(distrib, y, theta,
                                                      scale = c("parameter",
                                                                "link"), ...) {
  p <- .enet_parts(theta)
  zz <- .enet_logz_derivs(p)
  z <- y - p$mu
  s <- sign(z)
  al <- p$al
  lam <- p$lam
  d_a <- -abs(z) - zz$za
  d_c <- -z^2 / 2 - zz$zc
  list(
    mu = p$a * s + p$c * z,
    lambda = d_a * al + d_c * (1 - al),
    alpha = lam * (d_a - d_c)
  )
}

#' @title Elastic-Net Observed Hessian
#' @name distrib_hessian.EnetDistrib
#'
#' @description
#' Computes the six second derivatives of the log-density in closed form. The
#' data term of the log-density is **linear in the two rates**, so every rate
#' block is a second derivative of \eqn{\log Z} carried through the bilinear
#' map \eqn{(\lambda,\alpha) \mapsto (a, c)}, whose own cross term contributes
#' the \eqn{\lambda\alpha} entry. In the location the curvature is simply
#' \eqn{-c}.
#'
#' That \eqn{-c} is the observed curvature and **not** the information. It
#' misses the point mass \eqn{\mathrm{d}\,\mathrm{sgn}(z)/\mathrm{d}z =
#' 2\delta(z)} the absolute value carries at the location, exactly as the
#' Laplace does; see [distrib_expected_hessian.EnetDistrib()] for what the
#' information is instead.
#'
#' @param distrib An `EnetDistrib` object, from [enet_distrib()].
#' @param y A numeric vector of observations.
#' @param theta A named list with components `mu`, `lambda` and `alpha`, each a
#'   numeric vector of length 1 or of the length of `y`.
#' @param scale Either `"parameter"`, the default, or `"link"`. The
#'   transformation is applied in the generic's body.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A named list of six numeric vectors in [hess_names()]'s order:
#'   `mu_mu`, `lambda_lambda`, `alpha_alpha`, `mu_lambda`, `mu_alpha`,
#'   `lambda_alpha`. Two of them, `lambda_lambda` and `alpha_alpha`, carry no
#'   data at all and are their own expectations.
#'
#' @section Notation:
#' \eqn{a = \lambda\alpha}, \eqn{c = \lambda(1-\alpha)}, \eqn{z = y - \mu},
#' \eqn{x = a/\sqrt c}, and \eqn{Z} the normalizing constant.
#'
#' @seealso [distrib_gradient.EnetDistrib()] for the order below,
#'   [distrib_deriv3.EnetDistrib()] for the order above,
#'   [distrib_expected_hessian.EnetDistrib()] for the information, and
#'   [distrib_hessian()] for the generic.
#'
#' @examples
#' d <- enet_distrib()
#' y <- c(-1.5, -0.3, 0.4, 2.1)
#' th <- list(mu = 0, lambda = 2, alpha = 0.5)
#' h <- distrib_hessian(d, y, th)
#' names(h)
#'
#' # The curvature in the location is -c at every observation, the absolute
#' # value contributing nothing away from the kink.
#' c(observed = unique(h$mu_mu), minus_c = -(2 * (1 - 0.5)))
#'
#' # Two entries carry no data, so they equal their own expectations.
#' e <- distrib_expected_hessian(d, y, th)
#' rbind(observed = c(h$lambda_lambda[1], h$alpha_alpha[1], h$mu_mu[1]),
#'       expected = c(e$lambda_lambda[1], e$alpha_alpha[1], e$mu_mu[1]))
#'
#' # Against a central difference of the score.
#' eps <- 1e-5
#' rbind(analytic = h$lambda_alpha,
#'       numeric = (distrib_gradient(d, y, list(mu = 0, lambda = 2,
#'                                              alpha = 0.5 + eps))$lambda -
#'                  distrib_gradient(d, y, list(mu = 0, lambda = 2,
#'                                              alpha = 0.5 - eps))$lambda) /
#'                 (2 * eps))
S7::method(distrib_hessian, EnetDistrib) <- function(distrib, y, theta,
                                                     scale = c("parameter",
                                                               "link"), ...) {
  p <- .enet_parts(theta)
  zz <- .enet_logz_derivs(p)
  z <- y - p$mu
  s <- sign(z)
  al <- p$al
  lam <- p$lam
  one <- rep(1, length(z))
  list(
    mu_mu = -p$c * one,
    lambda_lambda = -(zz$zaa * al^2 + 2 * zz$zac * al * (1 - al) +
                      zz$zcc * (1 - al)^2) * one,
    alpha_alpha = -lam^2 * (zz$zaa - 2 * zz$zac + zz$zcc) * one,
    mu_lambda = al * s + (1 - al) * z,
    mu_alpha = lam * (s - z),
    lambda_alpha = (-abs(z) + z^2 / 2) -
      (lam * (zz$zaa * al + zz$zac * (1 - 2 * al) - zz$zcc * (1 - al)) +
       zz$za - zz$zc)
  )
}

#' @title Elastic-Net Expected Information
#' @name distrib_expected_hessian.EnetDistrib
#'
#' @description
#' Returns the expected Hessian in closed form. The family already carries
#' every piece of it: in the two rates the density is an exponential family, so
#' the information is the Hessian of its own log normalizing constant, and
#' `.enet_logz_derivs()` computes that for the observed Hessian already.
#'
#' `expected`, `approx` and `nsim` are all ignored, the answer being exact.
#'
#' @details
#' # The rate block
#'
#' With sufficient statistics \eqn{-|z|} and \eqn{-z^2/2}, \eqn{\log Z} is the
#' cumulant generating function and
#' \deqn{I_{aa} = \operatorname{Var}(|z|), \quad
#'   I_{ac} = \operatorname{Cov}\!\left(-|z|, -\tfrac{z^2}{2}\right), \quad
#'   I_{cc} = \tfrac{1}{4}\operatorname{Var}(z^2).}
#' The map to \eqn{(\lambda, \alpha)} is bilinear, so the information
#' transforms by \eqn{J^\top I J} with no second-derivative term.
#'
#' Two of the three rate entries of the observed Hessian carry no data and are
#' therefore their own expectations, which is why `lambda_lambda` and
#' `alpha_alpha` repeat them exactly. The third does not: `lambda_alpha`
#' carries \eqn{-|z| + z^2/2}, whose expectation
#' \eqn{\partial_a\log Z - \partial_c\log Z} cancels the constant of the same
#' value sitting beside it and leaves the \eqn{\log Z} term alone.
#'
#' # The location, where the kink lives
#'
#' The observed second derivative in \eqn{\mu} is \eqn{-c}, which misses the
#' point mass \eqn{\mathrm{d}\,\mathrm{sgn}(z)/\mathrm{d}z = 2\delta(z)} the
#' density carries at its own location, exactly as the Laplace does. The
#' information is defined there as the variance of the score, and
#' \deqn{I_{\mu\mu} = E\left[(a\,\mathrm{sgn}(z) + cz)^2\right]
#'   = a^2 + 2ac\,E|z| + c^2E[z^2]
#'   = a^2 - 2ac\,\partial_a \log Z - 2c^2 \partial_c \log Z,}
#' since \eqn{E|z| = -\partial_a\log Z} and \eqn{E[z^2] = -2\partial_c\log Z}.
#'
#' Measured at \eqn{\lambda = 2}: it is 3.99999996 as \eqn{\alpha \to 1}, which
#' is the Laplace's \eqn{\lambda^2 = 1/\sigma^2}, and exactly 2 as
#' \eqn{\alpha \to 0}, which is the Gaussian's \eqn{c}.
#'
#' The location stays **orthogonal** to both rates, every cross-expectation
#' vanishing by symmetry, so `mu_lambda` and `mu_alpha` are exactly zero.
#'
#' @param distrib An `EnetDistrib` object, from [enet_distrib()].
#' @param y A numeric vector. Its values do not enter the result, which is an
#'   expectation; only its length does, through recycling.
#' @param theta A named list with components `mu`, `lambda` and `alpha`.
#' @param scale Either `"parameter"`, the default, or `"link"`. The
#'   transformation is applied in the generic's body.
#' @param approx Ignored: the answer is exact. Accepted so that the signature
#'   matches the generic's.
#' @param nsim Ignored, for the same reason.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A named list of six numeric vectors in [hess_names()]'s order:
#'   `mu_mu`, `lambda_lambda`, `alpha_alpha`, `mu_lambda`, `mu_alpha`,
#'   `lambda_alpha`. The two mixed location entries are exactly zero.
#'
#' @section Notation:
#' \eqn{a = \lambda\alpha}, \eqn{c = \lambda(1-\alpha)}, \eqn{z = y-\mu} and
#' \eqn{Z} the normalizing constant.
#'
#' @seealso [distrib_hessian.EnetDistrib()] for the observed curvature, which
#'   differs from this only in the location;
#'   [distrib_expected_hessian.Laplace2Distrib()] for the same phenomenon in
#'   the family this one contains; and [distrib_expected_hessian()] for the
#'   generic.
#'
#' @examples
#' d <- enet_distrib()
#' th <- list(mu = 0, lambda = 2, alpha = 0.5)
#' e <- distrib_expected_hessian(d, 0, th)
#' names(e)
#'
#' # The location is orthogonal to both rates.
#' c(mu_lambda = e$mu_lambda, mu_alpha = e$mu_alpha)
#'
#' # The observed and expected curvature differ in the location and nowhere
#' # else among the entries that carry no data.
#' h <- distrib_hessian(d, 0, th)
#' rbind(observed = c(h$mu_mu, h$lambda_lambda, h$alpha_alpha),
#'       expected = c(e$mu_mu, e$lambda_lambda, e$alpha_alpha))
#'
#' # The location entry reaches both ends: lambda^2 at alpha -> 1, c at 0.
#' inf_mu <- function(a)
#'   -distrib_expected_hessian(d, 0, list(mu = 0, lambda = 2, alpha = a))$mu_mu
#' rbind(alpha = c(1 - 1e-8, 0.5, 1e-8),
#'       information = vapply(c(1 - 1e-8, 0.5, 1e-8), inf_mu, 0),
#'       limit = c(2^2, NA, 2 * (1 - 1e-8)))
#'
#' # The strategy argument is ignored, the answer being exact.
#' identical(e, distrib_expected_hessian(d, 0, th, approx = "mc", nsim = 5))
S7::method(distrib_expected_hessian, EnetDistrib) <- function(
    distrib, y, theta, scale = c("parameter", "link"),
    approx = c("bartlett", "integrate", "mc", "opg"), nsim = 10000, ...) {
  p <- .enet_parts(theta)
  zz <- .enet_logz_derivs(p)
  al <- p$al
  lam <- p$lam
  one <- rep(1, length(y))
  list(
    mu_mu = -(p$a^2 - 2 * p$a * p$c * zz$za - 2 * p$c^2 * zz$zc) * one,
    lambda_lambda = -(zz$zaa * al^2 + 2 * zz$zac * al * (1 - al) +
                      zz$zcc * (1 - al)^2) * one,
    alpha_alpha = -lam^2 * (zz$zaa - 2 * zz$zac + zz$zcc) * one,
    mu_lambda = 0 * one,
    mu_alpha = 0 * one,
    lambda_alpha = -lam * (zz$zaa * al + zz$zac * (1 - 2 * al) -
                           zz$zcc * (1 - al)) * one
  )
}

#' @title Elastic-Net Response Derivative
#' @name distrib_grad_y.EnetDistrib
#'
#' @description
#' Computes \eqn{\partial\ell/\partial y = -a\,\mathrm{sgn}(y-\mu) - c(y-\mu)},
#' with \eqn{a = \lambda\alpha} and \eqn{c = \lambda(1-\alpha)}. It is minus
#' the derivative in \eqn{\mu}, the response and the location entering only
#' through their difference.
#'
#' It is **undefined at the location**, where the sign jumps. `sign(0)` is 0 in
#' R, so the value returned there is \eqn{-c(y-\mu) = 0}, one point of the
#' subdifferential.
#'
#' The score is bounded below by \eqn{-a - c(y-\mu)} on either side, so the
#' Laplace part contributes a jump of \eqn{2a} across the location and the
#' Gaussian part the slope.
#'
#' @param distrib An `EnetDistrib` object, from [enet_distrib()].
#' @param y A numeric vector of observations.
#' @param theta A named list with components `mu`, `lambda` and `alpha`, each a
#'   numeric vector of length 1 or of the length of `y`.
#'
#' @return A numeric vector of the length of the recycled inputs.
#'
#' @section Notation:
#' \eqn{\ell} is the log-density of one observation, \eqn{a = \lambda\alpha}
#' and \eqn{c = \lambda(1-\alpha)}.
#'
#' @seealso [distrib_hess_y.EnetDistrib()] for the second derivative,
#'   [distrib_cross_y.EnetDistrib()] for the mixed one,
#'   [distrib_gradient.EnetDistrib()] for the parameter derivatives, and
#'   [distrib_grad_y()] for the generic.
#'
#' @examples
#' d <- enet_distrib()
#' y <- c(-1.5, -0.3, 0.4, 2.1)
#' th <- list(mu = 0, lambda = 2, alpha = 0.5)
#'
#' # It is minus the score in the location, exactly.
#' all.equal(distrib_grad_y(d, y, th), -distrib_gradient(d, y, th)$mu)
#'
#' # Against a central difference of the log-density, away from the kink.
#' eps <- 1e-6
#' rbind(analytic = distrib_grad_y(d, y, th),
#'       numeric = (distrib_pdf(d, y + eps, th, log = TRUE) -
#'                  distrib_pdf(d, y - eps, th, log = TRUE)) / (2 * eps))
#'
#' # The jump across the location is 2a, which is the Laplace part.
#' c(just_below = distrib_grad_y(d, -1e-9, th),
#'   just_above = distrib_grad_y(d, 1e-9, th),
#'   two_a = 2 * 2 * 0.5)
S7::method(distrib_grad_y, EnetDistrib) <- function(distrib, y, theta) {
  p <- .enet_parts(theta)
  z <- y - p$mu
  -p$a * sign(z) - p$c * z
}

#' @title Elastic-Net Second Response Derivative
#' @name distrib_hess_y.EnetDistrib
#'
#' @description
#' Computes \eqn{\partial^2\ell/\partial y^2 = -c}, with
#' \eqn{c = \lambda(1-\alpha)}: the absolute value contributes nothing away
#' from the location, so the curvature in the response is the Gaussian part
#' alone. It is constant, free of the data, and strictly negative.
#'
#' At the location itself the second derivative carries a point mass
#' \eqn{-2a\delta(y-\mu)}, which this value omits, exactly as
#' [distrib_hessian.EnetDistrib()]'s `mu_mu` does.
#'
#' Every response derivative of **third order and beyond is zero**: the
#' log-density is quadratic in \eqn{y} away from the kink.
#'
#' @param distrib An `EnetDistrib` object, from [enet_distrib()].
#' @param y A numeric vector of observations. Its values do not enter the
#'   result; only its length does.
#' @param theta A named list with components `mu`, `lambda` and `alpha`, each a
#'   numeric vector of length 1 or of the length of `y`.
#'
#' @return A numeric vector of the length of `y`, every entry \eqn{-c}.
#'
#' @section Notation:
#' \eqn{c = \lambda(1-\alpha)} is the Gaussian rate and \eqn{a = \lambda\alpha}
#' the Laplace one.
#'
#' @seealso [distrib_grad_y.EnetDistrib()] for the first derivative,
#'   [distrib_hessian.EnetDistrib()] for the parameter curvature it shares an
#'   expression with, and [distrib_hess_y()] for the generic.
#'
#' @examples
#' d <- enet_distrib()
#' y <- c(-1.5, -0.3, 0.4, 2.1)
#' th <- list(mu = 0, lambda = 2, alpha = 0.5)
#'
#' # Constant, and equal to the curvature in the location.
#' c(hess_y = unique(distrib_hess_y(d, y, th)),
#'   mu_mu = unique(distrib_hessian(d, y, th)$mu_mu),
#'   minus_c = -(2 * (1 - 0.5)))
#'
#' # Against a central difference of the response derivative, away from the
#' # kink.
#' eps <- 1e-5
#' rbind(analytic = distrib_hess_y(d, y, th),
#'       numeric = (distrib_grad_y(d, y + eps, th) -
#'                  distrib_grad_y(d, y - eps, th)) / (2 * eps))
#'
#' # The Gaussian rate alone, so it goes to zero as alpha approaches one.
#' vapply(c(0.1, 0.5, 0.9, 0.999),
#'        function(a) distrib_hess_y(d, 1, list(mu = 0, lambda = 2,
#'                                              alpha = a)), 0)
S7::method(distrib_hess_y, EnetDistrib) <- function(distrib, y, theta) {
  p <- .enet_parts(theta)
  rep(-p$c, length(y))
}

#' @title Elastic-Net Mixed Response-Parameter Derivatives
#' @name distrib_cross_y.EnetDistrib
#'
#' @description
#' Computes \eqn{\partial^2\ell/\partial y\,\partial\theta_i}, one component
#' per parameter, by differentiating
#' \eqn{\partial\ell/\partial y = -a\,\mathrm{sgn}(z) - cz} in each. With
#' \eqn{z = y-\mu} the components are \eqn{c} in the location,
#' \eqn{-\alpha\,\mathrm{sgn}(z) - (1-\alpha)z} in the rate and
#' \eqn{\lambda(z - \mathrm{sgn}(z))} in the mixing weight.
#'
#' The normalizing constant does not appear: it carries no \eqn{y}, so
#' differentiating in the response removes it before any parameter is
#' differentiated.
#'
#' This is the off-diagonal block of the joint Hessian in \eqn{(y, \theta)},
#' whose diagonals are [distrib_hess_y.EnetDistrib()] and
#' [distrib_hessian.EnetDistrib()]. `penalties7` consumes it: a penalty is a
#' negative log-density at the coefficients, and estimating coefficients and
#' hyperparameters together needs this block.
#'
#' @param distrib An `EnetDistrib` object, from [enet_distrib()].
#' @param y A numeric vector of observations. At `y == mu` the sign is 0 in R,
#'   so the value returned is one point of the subdifferential.
#' @param theta A named list with components `mu`, `lambda` and `alpha`, each a
#'   numeric vector of length 1 or of the length of `y`.
#' @param scale Either `"parameter"`, the default, or `"link"`. On the link
#'   scale the chain rule is the gradient's own first-order diagonal one, the
#'   response derivative not interacting with a reparametrization of `theta`;
#'   the transformation is applied in the generic's body.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A named list of three numeric vectors, `mu`, `lambda` and `alpha`,
#'   each of the length of the recycled inputs.
#'
#' @section Notation:
#' \eqn{z = y - \mu}, \eqn{a = \lambda\alpha}, \eqn{c = \lambda(1-\alpha)}.
#'
#' @seealso [distrib_grad_y.EnetDistrib()], which it differentiates,
#'   [distrib_hessian.EnetDistrib()] for the other diagonal block, and
#'   [distrib_cross_y()] for the generic.
#'
#' @examples
#' d <- enet_distrib()
#' y <- c(-1.5, -0.3, 0.4, 2.1)
#' th <- list(mu = 0, lambda = 2, alpha = 0.5)
#' cy <- distrib_cross_y(d, y, th)
#' names(cy)
#'
#' # Against a central difference of the response derivative in a parameter.
#' eps <- 1e-6
#' rbind(analytic = cy$lambda,
#'       numeric = (distrib_grad_y(d, y, list(mu = 0, lambda = 2 + eps,
#'                                            alpha = 0.5)) -
#'                  distrib_grad_y(d, y, list(mu = 0, lambda = 2 - eps,
#'                                            alpha = 0.5))) / (2 * eps))
#'
#' # The location component is the constant c, and equals minus the second
#' # response derivative.
#' c(cross_mu = unique(cy$mu), minus_hess_y = -unique(distrib_hess_y(d, y, th)))
S7::method(distrib_cross_y, EnetDistrib) <- function(distrib, y, theta,
                                                     scale = c("parameter",
                                                               "link"), ...) {
  p <- .enet_parts(theta)
  z <- y - p$mu
  s <- sign(z)
  list(
    mu = rep(p$c, length(z)),
    lambda = -p$al * s - (1 - p$al) * z,
    alpha = p$lam * (z - s)
  )
}

#' @title Mean of the Elastic-Net Distribution
#' @name mean.EnetDistrib
#'
#' @description
#' Returns \eqn{\mu}, the location. The density is symmetric about it, both
#' factors of the product being symmetric about the same point, so the location
#' is the mean, the median and the mode at once.
#'
#' The value is recycled to the length the three parameters imply, so a `theta`
#' whose components vary by observation gets one mean per observation.
#'
#' @param x An `EnetDistrib` object, from [enet_distrib()].
#' @param theta A named list with components `mu`, `lambda` and `alpha`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of means, of the length the recycled
#'   parameters imply. The value is the location; the two rates do not
#'   enter it, so a setting that varies one of them repeats one number.
#'
#' @seealso [variance.EnetDistrib()] and [skewness.EnetDistrib()] for the
#'   other two closed-form moments, [kurtosis()] for the fourth, which is not
#'   closed form here, and [distrib_quantile.EnetDistrib()], which confirms
#'   the median.
#'
#' @examples
#' d <- enet_distrib()
#' th <- list(mu = 3, lambda = 2, alpha = 0.5)
#'
#' # Mean, median and mode coincide.
#' c(mean = mean(d, th), median = distrib_quantile(d, 0.5, th))
#'
#' # One value per observation when a parameter varies.
#' mean(d, list(mu = c(0, 3, 7), lambda = 2, alpha = 0.5))
#'
#' @keywords internal
S7::method(mean, EnetDistrib) <- function(x, theta, ...) {
  .enet_parts(theta)$mu + moment_const(theta, 3L, 0)
}

#' @title Variance of the Elastic-Net Distribution
#' @name variance.EnetDistrib
#'
#' @description
#' Returns \eqn{(1 + xG)/c} in closed form, with \eqn{c = \lambda(1-\alpha)},
#' \eqn{x = a/\sqrt c} and \eqn{G = \mathrm{d}\log M/\mathrm{d}x}. It is
#' \eqn{-2\,\partial\log Z/\partial c}, the first cumulant of the sufficient
#' statistic \eqn{-z^2/2}, so it needs nothing the score does not already
#' compute.
#'
#' At \eqn{\alpha \to 0} it tends to \eqn{1/c}, the Gaussian's; at
#' \eqn{\alpha \to 1} to \eqn{2/a^2}, the Laplace's.
#'
#' @param x An `EnetDistrib` object, from [enet_distrib()].
#' @param theta A named list with components `mu`, `lambda` and `alpha`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of variances, strictly positive, of the
#'   length the recycled parameters imply. The location does not enter
#'   the value, so a setting that varies it repeats one number.
#'
#' @section Notation:
#' \eqn{a = \lambda\alpha}, \eqn{c = \lambda(1-\alpha)}, \eqn{x = a/\sqrt c},
#' \eqn{M} the Mills ratio and \eqn{Z} the normalizing constant.
#'
#' @seealso [mean.EnetDistrib()] and [skewness.EnetDistrib()], and
#'   [distrib_expected_hessian.EnetDistrib()], which uses the same
#'   \eqn{\log Z} derivatives.
#'
#' @examples
#' d <- enet_distrib()
#' th <- list(mu = 0, lambda = 2, alpha = 0.5)
#'
#' # Against a quadrature of the second central moment.
#' c(closed = variance(d, th),
#'   quadrature = integrate(function(u) u^2 * distrib_pdf(d, u, th),
#'                          -Inf, Inf)$value)
#'
#' # The two ends: 1/c at alpha -> 0 and 2/a^2 at alpha -> 1.
#' rbind(alpha = c(1e-10, 1 - 1e-10),
#'       ours = c(variance(d, list(mu = 0, lambda = 2, alpha = 1e-10)),
#'                variance(d, list(mu = 0, lambda = 2, alpha = 1 - 1e-10))),
#'       limit = c(1 / 2, 2 / 2^2))
#'
#' @keywords internal
S7::method(variance, EnetDistrib) <- function(x, theta, ...) {
  p <- .enet_parts(theta)
  (1 + p$x * p$g) / p$c + moment_const(theta, 3L, 0)
}

#' @title Skewness of the Elastic-Net Distribution
#' @name skewness.EnetDistrib
#'
#' @description
#' Returns 0. The density is symmetric about \eqn{\mu} at every parameter
#' value, so every odd central moment vanishes and the standardized third one
#' with them. The result is multiplied by a length-carrying quantity, so a
#' `theta` whose components vary by observation gets one zero per observation.
#'
#' The **fourth** moment is not zero and is not closed form here: [kurtosis()]
#' falls to the base class and integrates. Measured at \eqn{\mu = 0},
#' \eqn{\lambda = 2}, \eqn{\alpha = 0.5} the excess kurtosis is 0.766, between
#' the Gaussian's 0 and the Laplace's 3.
#'
#' @param x An `EnetDistrib` object, from [enet_distrib()].
#' @param theta A named list with components `mu`, `lambda` and `alpha`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of zeros, of the length the recycled
#'   parameters imply. The family is symmetric about its location at
#'   every setting of the two rates, so no parameter enters the value.
#'
#' @seealso [mean.EnetDistrib()] and [variance.EnetDistrib()] for the other two
#'   closed-form moments, and [kurtosis()] for the fourth.
#'
#' @examples
#' d <- enet_distrib()
#' th <- list(mu = 0, lambda = 2, alpha = 0.5)
#'
#' # Zero at every parameter value, and confirmed by a quadrature.
#' c(closed = skewness(d, th),
#'   quadrature = integrate(function(u) u^3 * distrib_pdf(d, u, th),
#'                          -Inf, Inf)$value)
#'
#' # The fourth moment is not zero, and sits between the two ends.
#' c(gaussian = 0, enet = kurtosis(d, th), laplace = 3)
#'
#' @keywords internal
S7::method(skewness, EnetDistrib) <- function(x, theta, ...) {
  moment_const(theta, 3L, 0)
}

#' @title Elastic-Net Distribution Object
#'
#' @description
#' Builds the density whose negative logarithm is the elastic-net penalty: a
#' Laplace and a Gaussian at the same location, multiplied together and
#' normalized. It exists so that [penalties7::elasticnet_penalty()] is the same
#' construction as ridge and lasso, a [penalties7::distrib_penalty()] over a
#' [fixed()] family, instead of a branch of its own.
#'
#' The parametrization is the one a reader of the penalty expects: an overall
#' rate \eqn{\lambda} and a mixing weight \eqn{\alpha}, with
#' \eqn{a = \lambda\alpha} the Laplace rate and \eqn{c = \lambda(1-\alpha)} the
#' Gaussian one.
#'
#' @param link_mu A `linkfunctions7` link object for the location \eqn{\mu},
#'   which is unconstrained. Defaults to [linkfunctions7::identity_link()].
#' @param link_lambda A link object for the overall rate \eqn{\lambda}, which
#'   must be strictly positive. Defaults to [linkfunctions7::log_link()].
#' @param link_alpha A link object for the mixing weight \eqn{\alpha}, which
#'   must lie strictly inside \eqn{(0, 1)}. Defaults to
#'   [linkfunctions7::logit_link()].
#'
#' @details
#' # The normalizing constant
#'
#' \eqn{Z = 2M(a/\sqrt{c})/\sqrt{c}} with \eqn{M} the Mills ratio, finite at
#' both ends: \eqn{2/a} as \eqn{\alpha \to 1} and \eqn{\sqrt{2\pi/c}} as
#' \eqn{\alpha \to 0}. Every derivative in the two rates is a polynomial in
#' \eqn{x = a/\sqrt c} and \eqn{G = \mathrm{d}\log M/\mathrm{d}x}, with
#' \eqn{G' = 1 + xG - G^2} closing the recursion, and the chain to
#' \eqn{(\lambda, \alpha)} is bilinear.
#'
#' # What is closed form
#'
#' All four derivative orders in the parameters, the expected information, the
#' distribution function, the quantile function, the mean, the variance and the
#' skewness. The **kurtosis** is not, and comes from the base class by
#' quadrature.
#'
#' The expected information is closed form for a reason worth naming: in the
#' two rates the density is an exponential family with sufficient statistics
#' \eqn{-|z|} and \eqn{-z^2/2}, so the information there is the Hessian of
#' \eqn{\log Z}, which the observed Hessian already computes. See
#' [distrib_expected_hessian.EnetDistrib()].
#'
#' # The kink
#'
#' Like the Laplace, the log-likelihood is not differentiable in \eqn{\mu}.
#' `params_smooth` records `mu = FALSE`, the observed curvature there is
#' \eqn{-c} and misses the point mass at the location, and the information is
#' defined as the variance of the score instead. The two agree at every other
#' entry.
#'
#' # Two numerical cautions
#'
#' Both were found by sweeping the parameters, and both bite at ordinary
#' settings.
#'
#' - \eqn{\log M(x)} written directly adds \eqn{x^2/2} to a log-probability of
#'   the same size and opposite sign, losing a digit per factor of ten in
#'   \eqn{x}. At \eqn{\alpha = 1 - 10^{-12}} the argument reaches \eqn{10^6}
#'   and the density was wrong in the fourth digit; past \eqn{x = 30} the
#'   asymptotic series of the Mills ratio is used instead.
#' - \eqn{\Phi(-x)} is exactly zero past \eqn{x = 38}, and the distribution
#'   function and the quantile divide by it. At \eqn{\lambda = 20},
#'   \eqn{\alpha = 0.995} the argument is 63, so both work through `log.p` and
#'   `qnorm(..., log.p = TRUE)`.
#'
#' # Comparison with glmnet
#'
#' `glmnet` standardizes the response internally for the Gaussian family and
#' corrects `lambda` by a single factor. That correction is exact for the
#' \eqn{\ell_1} part, homogeneous of degree one, and cannot be for the
#' \eqn{\ell_2} part, homogeneous of degree two. So a nominal
#' \eqn{(\lambda, \alpha)} names two different objectives in the two packages
#' whenever \eqn{\alpha < 1}. At \eqn{\alpha = 1} the two agree to
#' \eqn{1.4\times10^{-4}} whatever the scale of the response.
#'
#' # Parameter domains
#'
#' - \eqn{\mu \in (-\infty, \infty)}
#' - \eqn{\lambda \in (0, \infty)}
#' - \eqn{\alpha \in (0, 1)}, open at both ends: at 1 it is
#'   [laplace2_distrib()] and at 0 [gaussian1_distrib()], both of which remain
#'   families of their own.
#'
#' @section The distribution:
#' \deqn{f(y) = \frac{1}{Z}\exp\!\left\{-a|y-\mu| - \frac{c}{2}(y-\mu)^{2}\right\}, \qquad a = \lambda\alpha, \quad c = \lambda(1-\alpha), \quad Z = \frac{2M(a/\sqrt{c})}{\sqrt{c}}}
#' on \eqn{y \in \mathbb{R}}.
#'
#' \deqn{\mathbb{E}[Y] = \mu, \qquad \operatorname{Var}(Y) = \frac{1 + xG}{c}, \qquad x = \frac{a}{\sqrt{c}}, \quad G = \frac{\mathrm{d}\log M}{\mathrm{d}x}}
#'
#' @return An S7 object of class [EnetDistrib], inheriting from
#'   `continuous_distrib`. Its `params` are `c("mu", "lambda", "alpha")`, its
#'   `bounds` `c(-Inf, Inf)`, its `params_smooth`
#'   `c(mu = FALSE, lambda = TRUE, alpha = TRUE)`, and its `link_params` the
#'   three links given here.
#'
#' @references
#' Zou, H. and Hastie, T. (2005). Regularization and variable selection
#' via the elastic net. *Journal of the Royal Statistical Society,
#' Series B* 67, 301-320.
#'
#' @seealso [laplace2_distrib()] and [gaussian1_distrib()] for the two ends,
#'   [penalties7::elasticnet_penalty()] for the consumer,
#'   [fixed()] for holding the location at zero, which is how a penalty uses
#'   this, and [EnetDistrib] for the class and its method list.
#'
#' @examples
#' d <- enet_distrib()
#' th <- list(mu = 0, lambda = 2, alpha = 0.5)
#'
#' distrib_pdf(d, c(-1, 0, 1), th)
#'
#' # It integrates to one, and the two ends are the two families.
#' integrate(function(u) distrib_pdf(d, u, th), -Inf, Inf)$value
#' rbind(enet = c(distrib_pdf(d, 0.7, list(mu = 0, lambda = 2,
#'                                         alpha = 1 - 1e-10)),
#'                distrib_pdf(d, 0.7, list(mu = 0, lambda = 2,
#'                                         alpha = 1e-10))),
#'       limit = c(distrib_pdf(laplace2_distrib(), 0.7,
#'                             list(mu = 0, lambda = 2)),
#'                 dnorm(0.7, 0, 1 / sqrt(2))))
#'
#' # The kurtosis sits between the two ends, and is the one moment that is
#' # not closed form.
#' c(gaussian = 0, enet = kurtosis(d, th), laplace = 3)
#'
#' # The location carries a kink, so its observed curvature and its
#' # information differ; every other entry agrees.
#' rbind(observed = unlist(distrib_hessian(d, 0, th)),
#'       expected = unlist(distrib_expected_hessian(d, 0, th)))
#'
#' # Held at zero, this is the prior a penalty is written from.
#' p <- fixed(d, mu = 0)
#' p@params
#'
#' @export
enet_distrib <- function(link_mu = identity_link(),
                         link_lambda = log_link(),
                         link_alpha = logit_link()) {
  EnetDistrib(
    distrib_name = "enet",
    dimension = "univariate",
    bounds = c(-Inf, Inf),

    params = c("mu", "lambda", "alpha"),
    params_interpretation = c(mu = "location", lambda = "rate",
                              alpha = "mixing weight"),
    n_params = 3,

    params_bounds = list(
      mu = c(-Inf, Inf),
      lambda = c(0, Inf),
      alpha = c(0, 1)
    ),

    link_params = list(
      mu = link_mu,
      lambda = link_lambda,
      alpha = link_alpha
    ),

    params_smooth = c(mu = FALSE, lambda = TRUE, alpha = TRUE)
  )
}
