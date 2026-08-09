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

#' @title S7 Class for the Elastic-Net Distribution
#' @name EnetDistrib
#'
#' @description A subclass of \code{continuous_distrib} for the density
#' whose negative logarithm is the elastic-net penalty,
#' \eqn{f(y) \propto \exp\{-\lambda\alpha|y-\mu| -
#' \lambda(1-\alpha)(y-\mu)^2/2\}}: the product of a Laplace and a
#' Gaussian, normalized. At \eqn{\alpha \to 1} it is the Laplace of
#' \code{\link{Laplace2Distrib}} and at \eqn{\alpha \to 0} the Gaussian,
#' both of which remain families of their own; \eqn{\alpha} is confined
#' to the open interval, as every bounded parameter in this package is.
#' Like the Laplace, its log-likelihood is not differentiable in
#' \eqn{\mu}.
#' @inheritParams distrib
#' @return An object of class \code{EnetDistrib}.
#' @seealso \code{\link{enet_distrib}}, \code{\link{laplace2_distrib}}
#'
#' @section Methods:
#' Methods implemented for this class:
#'   \code{\link[=distrib_cdf.EnetDistrib]{distrib_cdf()}},
#'   \code{\link[=distrib_cross_y.EnetDistrib]{distrib_cross_y()}},
#'   \code{\link[=distrib_grad_y.EnetDistrib]{distrib_grad_y()}},
#'   \code{\link[=distrib_gradient.EnetDistrib]{distrib_gradient()}},
#'   \code{\link[=distrib_hess_y.EnetDistrib]{distrib_hess_y()}},
#'   \code{\link[=distrib_hessian.EnetDistrib]{distrib_hessian()}},
#'   \code{\link[=distrib_pdf.EnetDistrib]{distrib_pdf()}},
#'   \code{\link[=distrib_quantile.EnetDistrib]{distrib_quantile()}},
#'   \code{\link[=distrib_rng.EnetDistrib]{distrib_rng()}},
#'   \code{\link[=mean.distrib]{mean()}},
#'   \code{\link[=skewness]{skewness()}},
#'   \code{\link[=variance]{variance()}}
#'
#' Everything else is inherited from \code{\link{continuous_distrib}}.
EnetDistrib <- S7::new_class("EnetDistrib", parent = continuous_distrib)

# --- S7 METHODS IMPLEMENTATION ---

#' @title Elastic-Net Density
#' @name distrib_pdf.EnetDistrib
#' @description
#' \deqn{f(y; \mu, \lambda, \alpha) = \frac{1}{Z}
#'   \exp\left\{-a|y-\mu| - \tfrac{c}{2}(y-\mu)^2\right\},}
#' with \eqn{a = \lambda\alpha}, \eqn{c = \lambda(1-\alpha)} and
#' \eqn{Z = 2M(a/\sqrt{c})/\sqrt{c}}, where \eqn{M} is the Mills ratio.
#' The constant is evaluated through the log Mills ratio, both factors of
#' which underflow separately.
#' @param distrib An \code{EnetDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing \code{mu}, \code{lambda} and \code{alpha}.
#' @param log Logical; if \code{TRUE}, returns the log-density.
#' @return A numeric vector of density values.
#' @seealso \code{\link{enet_distrib}}
S7::method(distrib_pdf, EnetDistrib) <- function(distrib, y, theta,
                                                 log = FALSE) {
  p <- .enet_parts(theta)
  z <- y - p$mu
  log_z <- base::log(2) - 0.5 * base::log(p$c) + .enet_logM(p$x)
  log_d <- -p$a * abs(z) - p$c * z^2 / 2 - log_z
  if (log) log_d else exp(log_d)
}

#' @title Elastic-Net Distribution Function
#' @name distrib_cdf.EnetDistrib
#' @description
#' Each half of the density is a truncated Gaussian, so
#' \eqn{F(q) = \Phi(\sqrt{c}\,z - x) / (2\Phi(-x))} for \eqn{z \le 0} and
#' its reflection above, with \eqn{z = q-\mu} and \eqn{x = a/\sqrt{c}}.
#' @param distrib An \code{EnetDistrib} object.
#' @param q A numeric vector of quantiles.
#' @param theta A list containing \code{mu}, \code{lambda} and \code{alpha}.
#' @param lower.tail Logical; if \code{TRUE} (default), \eqn{P(Y \le q)}.
#' @param log.p Logical; if \code{TRUE}, returns the log-probability.
#' @return A numeric vector of probabilities.
#' @seealso \code{\link{enet_distrib}}
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
#' @description The distribution function of
#' \code{\link[=distrib_cdf.EnetDistrib]{distrib_cdf()}} inverted in
#' closed form, each half being a truncated Gaussian.
#' @param distrib An \code{EnetDistrib} object.
#' @param p A numeric vector of probabilities.
#' @param theta A list containing \code{mu}, \code{lambda} and \code{alpha}.
#' @param lower.tail Logical; if \code{TRUE} (default), probabilities are
#'   \eqn{P(Y \le q)}.
#' @param log.p Logical; if \code{TRUE}, \code{p} is on the log scale.
#' @return A numeric vector of quantiles.
#' @seealso \code{\link{enet_distrib}}
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
#' @description Inverse transform on the closed-form quantile function.
#' @param distrib An \code{EnetDistrib} object.
#' @param n The number of draws.
#' @param theta A list containing \code{mu}, \code{lambda} and \code{alpha}.
#' @return A numeric vector of draws.
#' @seealso \code{\link{enet_distrib}}
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
#' @description
#' With \eqn{z = y-\mu}, \eqn{\ell^{(\mu)} = a\,\mathrm{sgn}(z) + cz}
#' and the two rate components are the data terms less the derivatives of
#' \eqn{\log Z}, which are \eqn{G/\sqrt{c}} and
#' \eqn{-(1+xG)/(2c)} for \eqn{a} and \eqn{c}; the chain to
#' \eqn{(\lambda, \alpha)} is linear, \eqn{a = \lambda\alpha} and
#' \eqn{c = \lambda(1-\alpha)}.
#' @param distrib An \code{EnetDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing \code{mu}, \code{lambda} and \code{alpha}.
#' @param scale Either \code{"parameter"} or \code{"link"}.
#' @param ... Ignored.
#' @return A named list of score components.
#' @seealso \code{\link{enet_distrib}}
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

#' @title Elastic-Net Observed Information
#' @name distrib_hessian.EnetDistrib
#' @description The second derivatives of the log-density. The data term
#' is linear in the two rates, so every rate block is a second derivative
#' of \eqn{\log Z} carried through the bilinear map
#' \eqn{(\lambda,\alpha) \mapsto (a, c)}, whose own cross term
#' contributes the \eqn{\lambda\alpha} entry.
#' @param distrib An \code{EnetDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing \code{mu}, \code{lambda} and \code{alpha}.
#' @param scale Either \code{"parameter"} or \code{"link"}.
#' @param ... Ignored.
#' @return A named list of Hessian components.
#' @seealso \code{\link{enet_distrib}}
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

#' @title Elastic-Net Response Derivatives
#' @name distrib_grad_y.EnetDistrib
#' @description \eqn{\ell^{(y)} = -a\,\mathrm{sgn}(y-\mu) - c(y-\mu)},
#' undefined at the location.
#' @param distrib An \code{EnetDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing \code{mu}, \code{lambda} and \code{alpha}.
#' @return A numeric vector.
#' @seealso \code{\link{enet_distrib}}
S7::method(distrib_grad_y, EnetDistrib) <- function(distrib, y, theta) {
  p <- .enet_parts(theta)
  z <- y - p$mu
  -p$a * sign(z) - p$c * z
}

#' @title Elastic-Net Second Response Derivative
#' @name distrib_hess_y.EnetDistrib
#' @description \eqn{\ell^{(yy)} = -c}: the absolute value contributes
#' nothing away from the location, so the curvature in the response is the
#' Gaussian part alone.
#' @param distrib An \code{EnetDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing \code{mu}, \code{lambda} and \code{alpha}.
#' @return A numeric vector.
#' @seealso \code{\link{enet_distrib}}
S7::method(distrib_hess_y, EnetDistrib) <- function(distrib, y, theta) {
  p <- .enet_parts(theta)
  rep(-p$c, length(y))
}

#' @title Elastic-Net Mixed Response-Parameter Derivatives
#' @name distrib_cross_y.EnetDistrib
#' @description \eqn{\partial^2\ell/\partial y\,\partial\theta_i},
#' obtained by differentiating \eqn{\ell^{(y)} = -a\,\mathrm{sgn}(z) - cz}
#' in each parameter.
#' @param distrib An \code{EnetDistrib} object.
#' @param y A numeric vector of observations.
#' @param theta A list containing \code{mu}, \code{lambda} and \code{alpha}.
#' @param scale Either \code{"parameter"} or \code{"link"}.
#' @param ... Ignored.
#' @return A named list, one component per parameter.
#' @seealso \code{\link{enet_distrib}}
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

#' @title Elastic-Net Moments
#' @name mean.EnetDistrib
#' @description The density is symmetric about \eqn{\mu}, so the mean is
#' the location and the skewness zero. The variance is
#' \eqn{(1+xG)/c}, which is \eqn{-2\,\partial \log Z/\partial c}.
#' @param x An \code{EnetDistrib} object.
#' @param theta A list containing \code{mu}, \code{lambda} and \code{alpha}.
#' @param ... Ignored.
#' @return A numeric vector.
#' @seealso \code{\link{enet_distrib}}
S7::method(mean, EnetDistrib) <- function(x, theta, ...) {
  .enet_parts(theta)$mu
}

#' @rdname mean.EnetDistrib
#' @name variance.EnetDistrib
S7::method(variance, EnetDistrib) <- function(x, theta, ...) {
  p <- .enet_parts(theta)
  (1 + p$x * p$g) / p$c
}

#' @rdname mean.EnetDistrib
#' @name skewness.EnetDistrib
S7::method(skewness, EnetDistrib) <- function(x, theta, ...) {
  0 * .enet_parts(theta)$c
}

#' @title Construct an Elastic-Net Distribution
#'
#' @description
#' The density whose negative logarithm is the elastic-net penalty: the
#' product of a Laplace and a Gaussian at the same location, normalized.
#' It exists so that \code{\link[penalties7]{elasticnet_penalty}} is the
#' same construction as ridge and lasso -- a
#' \code{\link[penalties7]{distrib_penalty}} over a
#' \code{\link{fixed}} family -- rather than a branch of its own.
#'
#' @details
#' Writing \eqn{a = \lambda\alpha} and \eqn{c = \lambda(1-\alpha)}, the
#' normalizing constant is \eqn{Z = 2M(a/\sqrt{c})/\sqrt{c}} with \eqn{M}
#' the Mills ratio, which is finite at both ends: \eqn{2/a} as
#' \eqn{\alpha \to 1} and \eqn{\sqrt{2\pi/c}} as \eqn{\alpha \to 0}. Every
#' derivative in the two rates is a polynomial in \eqn{x = a/\sqrt{c}} and
#' \eqn{G = \mathrm{d}\log M/\mathrm{d}x}, and the chain to
#' \eqn{(\lambda, \alpha)} is bilinear.
#'
#' Orders three and four in the parameters come from the numerical
#' fallback; the response derivatives are exact at every order, the third
#' and beyond being zero.
#'
#' @section The distribution:
#' \deqn{f(y) = \frac{1}{Z}\exp\!\left\{-a|y-\mu| - \frac{c}{2}(y-\mu)^{2}\right\}, \qquad a = \lambda\alpha, \quad c = \lambda(1-\alpha), \quad Z = \frac{2M(a/\sqrt{c})}{\sqrt{c}}}
#' on \eqn{y \in \mathbb{R}}.
#'
#' \deqn{\mathbb{E}[Y] = \mu, \qquad \operatorname{Var}(Y) = \frac{1 + xG}{c}, \qquad x = \frac{a}{\sqrt{c}}, \quad G = \frac{\mathrm{d}\log M}{\mathrm{d}x}}
#'
#' @param link_mu A \pkg{linkfunctions7} link for the location.
#' @param link_lambda A link for the overall rate, positive.
#' @param link_alpha A link for the mixing weight, in \eqn{(0,1)}.
#'
#' @return An object of class \code{\link{EnetDistrib}}.
#'
#' @examples
#' d <- enet_distrib()
#' theta <- list(mu = 0, lambda = 2, alpha = 0.5)
#' distrib_pdf(d, c(-1, 0, 1), theta)
#'
#' # the density integrates to one, and the two ends are the two families
#' stats::integrate(function(u) distrib_pdf(d, u, theta), -Inf, Inf)$value
#' distrib_pdf(d, 0.7, list(mu = 0, lambda = 2, alpha = 1 - 1e-10))
#' distrib_pdf(laplace2_distrib(), 0.7, list(mu = 0, lambda = 2))
#'
#' @seealso \code{\link{laplace2_distrib}}, \code{\link{gaussian1_distrib}}
#' @references
#' Zou, H. and Hastie, T. (2005). Regularization and variable selection
#' via the elastic net. \emph{Journal of the Royal Statistical Society,
#' Series B} 67, 301-320.
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
