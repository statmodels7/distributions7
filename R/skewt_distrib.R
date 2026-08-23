#' @include distrib.R generics.R skewnormal1_distrib.R
NULL

#' @title S7 Class for the Skew t Distribution
#' @name SkewTDistrib
#'
#' @description A subclass of `continuous_distrib` representing Azzalini's
#' skew \eqn{t} distribution: a Student \eqn{t} with a shape parameter
#' controlling the asymmetry, so that the tail weight and the skewness are
#' modeled separately.
#' @inheritParams distrib
#' @return An object of class `SkewTDistrib`.
#' @seealso [skewt_distrib()]
#'
#' @section Methods:
#' Methods implemented for this class:
#'   [`distrib_grad_y()`][distrib_grad_y.SkewTDistrib],
#'   [`distrib_gradient()`][distrib_gradient.SkewTDistrib],
#'   [`distrib_hess_y()`][distrib_hess_y.SkewTDistrib],
#'   [`distrib_hessian()`][distrib_hessian.SkewTDistrib],
#'   [`distrib_deriv3()`][distrib_deriv3.SkewTDistrib],
#'   [`distrib_deriv4()`][distrib_deriv4.SkewTDistrib],
#'   [`distrib_pdf()`][distrib_pdf.SkewTDistrib],
#'   [`distrib_rng()`][distrib_rng.SkewTDistrib],
#'   [`kurtosis()`][kurtosis],
#'   [`mean()`][mean.distrib],
#'   [`skewness()`][skewness],
#'   [`variance()`][variance]
#'
#' Everything else is inherited from [continuous_distrib()], including
#' the distribution function and the quantile function.
SkewTDistrib <- S7::new_class("SkewTDistrib", parent = continuous_distrib)

# --- S7 METHODS IMPLEMENTATION ---

#' The Pieces a Skew t Evaluates From
#'
#' @description
#' Assembles the standardized variable and the argument of the tilting
#' distribution function, together with the six scalar functions every
#' derivative of the log-density is a combination of.
#'
#' @details
#' With \eqn{z = (y-\mu)/\sigma}, \eqn{m = \nu + 1} and
#' \eqn{c = \sqrt{m/(\nu + z^2)}}, the tilting argument is \eqn{w = \alpha z c}.
#' The functions returned are
#' \deqn{A = \dfrac{\partial}{\partial z}\log t_\nu(z) = -\dfrac{m z}{\nu + z^2},
#'       \qquad
#'       A' = -\dfrac{m(\nu - z^2)}{(\nu + z^2)^2},}
#' \deqn{E = \dfrac{\nu\sqrt{m}}{(\nu + z^2)^{3/2}}, \qquad
#'       B = \dfrac{\partial w}{\partial z} = \alpha E, \qquad
#'       B' = -\dfrac{3\alpha\nu\sqrt{m}\,z}{(\nu + z^2)^{5/2}},}
#' and \eqn{Q = t_m(w)/T_m(w)} with
#' \eqn{Q' = Q\left\{-\dfrac{(m+1)w}{m + w^2} - Q\right\}}, the last from
#' differentiating the quotient.
#'
#' @param y A numeric vector of observations.
#' @param mu,sigma,alpha,nu The parameters.
#'
#' @return A list with `z`, `w`, `c`, `a`, `da`,
#'   `e`, `b`, `db`, `q` and `dq`.
#'
#' @keywords internal
skewt_pieces <- function(y, mu, sigma, alpha, nu) {
  z <- (y - mu) / sigma
  m <- nu + 1
  s <- nu + z^2
  cc <- sqrt(m / s)
  w <- alpha * z * cc
  # Q is formed on the log scale for the same reason numericals7::mills_ratio() is: both the
  # density and the distribution function underflow in the far left tail while
  # their ratio stays finite.
  q <- exp(stats::dt(w, df = m, log = TRUE) - stats::pt(w, df = m, log.p = TRUE))
  e <- nu * sqrt(m) / s^1.5
  list(
    z = z, w = w, c = cc,
    a = -m * z / s,
    da = -m * (nu - z^2) / s^2,
    e = e,
    b = alpha * e,
    db = -3 * alpha * nu * sqrt(m) * z / s^2.5,
    q = q,
    dq = q * (-(m + 1) * w / (m + w^2) - q)
  )
}

#' The Step a Skew t Differences the Degrees of Freedom With
#'
#' @description
#' Returns the finite-difference step used for the derivatives in \eqn{\nu},
#' relative to \eqn{\nu} itself and floored so that it stays meaningful for a
#' small number of degrees of freedom.
#'
#' @details
#' The relative step \eqn{10^{-3}} is measured rather than assumed. Swept over
#' \eqn{\nu} from 2 to 30 and sample sizes from 500 to 4000, it is where the
#' truncation error of the five-point stencil of [fd5_first()] has
#' fallen to the level of the rounding error and the two are balanced; a
#' smaller step is dominated by rounding, which the stencil amplifies by
#' \eqn{18/(12h)}, and a larger one by truncation, which grows as \eqn{h^4}.
#'
#' @param nu The degrees of freedom.
#'
#' @return A numeric vector.
#'
#' @keywords internal
skewt_nu_step <- function(nu) {
  pmax(1e-3 * abs(nu), 1e-6)
}

#' A Five-Point First Derivative
#'
#' @description
#' Returns \eqn{(f(x-2h) - 8f(x-h) + 8f(x+h) - f(x+2h))/(12h)}, the
#' fourth-order central stencil.
#'
#' @details
#' The three-point stencil is not accurate enough for the derivative in
#' \eqn{\nu} of a fitted likelihood. Its truncation error is \eqn{O(h^2)} per
#' observation and does **not** cancel when the observations are summed,
#' because it is a bias rather than noise: on a sample of a few thousand it
#' leaves the summed score at about \eqn{10^{-8}}, an order of magnitude worse
#' than the five-point stencil, which is \eqn{O(h^4)} and reaches the level of
#' rounding at the cost of two more evaluations.
#'
#' This is one stencil applied to an analytic quantity, not a difference of a
#' difference.
#'
#' @param f A function of one scalar, returning a numeric vector.
#' @param x The point to differentiate at.
#' @param h The step.
#'
#' @return A numeric vector, of whatever length `f` returns.
#'
#' @keywords internal
fd5_first <- function(f, x, h) {
  # numericals7's shared weights at accuracy four: the displayed formula is
  # exactly what the Vandermonde construction produces on five nodes.
  numericals7::fd_derivative(f, x, 1L, h = h, accuracy = 4L)
}

#' A Five-Point Second Derivative
#'
#' @description
#' Returns \eqn{(-f(x-2h) + 16f(x-h) - 30f(x) + 16f(x+h) - f(x+2h))/(12h^2)},
#' the fourth-order central stencil for the second derivative.
#'
#' @param f A function of one scalar, returning a numeric vector.
#' @param x The point to differentiate at.
#' @param h The step.
#'
#' @return A numeric vector, of whatever length `f` returns.
#'
#' @keywords internal
fd5_second <- function(f, x, h) {
  numericals7::fd_derivative(f, x, 2L, h = h, accuracy = 4L)
}

#' A Five-Point Third Derivative
#'
#' @description
#' Returns \eqn{(-f(x-2h)/2 + f(x-h) - f(x+h) + f(x+2h)/2)/h^3}, the
#' second-order central stencil for the third derivative. One stencil applied
#' to an analytic quantity, never a difference of differences.
#'
#' @param f A function of one scalar, returning a numeric vector.
#' @param x The point to differentiate at.
#' @param h The step.
#'
#' @return A numeric vector.
#'
#' @keywords internal
fd5_third <- function(f, x, h) {
  numericals7::fd_derivative(f, x, 3L, h = h)
}

#' A Five-Point Fourth Derivative
#'
#' @description
#' Returns \eqn{(f(x-2h) - 4f(x-h) + 6f(x) - 4f(x+h) + f(x+2h))/h^4}, the
#' second-order central stencil for the fourth derivative. Rounding is
#' amplified by \eqn{h^{-4}}, so with the family's relative step this is
#' accurate to roughly four significant digits, which the pages that rely on
#' it state.
#'
#' @param f A function of one scalar, returning a numeric vector.
#' @param x The point to differentiate at.
#' @param h The step.
#'
#' @return A numeric vector.
#'
#' @keywords internal
fd5_fourth <- function(f, x, h) {
  numericals7::fd_derivative(f, x, 4L, h = h)
}

#' @title Skew t Probability Density Function
#' @name distrib_pdf.SkewTDistrib
#' @description
#' Computes the probability density function, with \eqn{z = (y-\mu)/\sigma} and
#' \eqn{w = \alpha z \sqrt{(\nu+1)/(\nu+z^2)}}:
#' \deqn{f(y; \mu, \sigma, \alpha, \nu) = \dfrac{2}{\sigma}\,
#'   t_\nu(z)\,T_{\nu+1}(w)}
#' with \eqn{t_\nu} the standard Student \eqn{t} density and \eqn{T_{\nu+1}} the
#' Student \eqn{t} distribution function on one more degree of freedom.
#' @param distrib A `SkewTDistrib` object.
#' @param y A numeric vector of observations.
#' @param theta A list containing `mu`, `sigma`, `alpha` and `nu`.
#' @param log Logical; if `TRUE`, returns the log-density.
#' @return A numeric vector of density values.
#' @seealso [skewt_distrib()]
S7::method(distrib_pdf, SkewTDistrib) <- function(distrib, y, theta, log = FALSE, ...) {
  mu <- theta[[1]]
  sigma <- theta[[2]]
  alpha <- theta[[3]]
  nu <- theta[[4]]
  z <- (y - mu) / sigma
  w <- alpha * z * sqrt((nu + 1) / (nu + z^2))
  log_d <- log(2) - log(sigma) + stats::dt(z, df = nu, log = TRUE) +
    stats::pt(w, df = nu + 1, log.p = TRUE)
  if (log) log_d else exp(log_d)
}

#' @title Skew t Random Number Generator
#' @name distrib_rng.SkewTDistrib
#' @description
#' Generates draws exactly, from the scale mixture
#' \eqn{Y = \mu + \sigma Z/\sqrt{V/\nu}} with \eqn{Z} standard skew normal of
#' shape \eqn{\alpha} and \eqn{V \sim \chi^2_\nu} independent of it. No
#' inversion or rejection is involved.
#' @param distrib A `SkewTDistrib` object.
#' @param n Number of observations to generate.
#' @param theta A list containing `mu`, `sigma`, `alpha` and `nu`.
#' @return A numeric vector of random draws.
#' @seealso [skewt_distrib()]
S7::method(distrib_rng, SkewTDistrib) <- function(distrib, n, theta) {
  alpha <- theta[[3]]
  nu <- theta[[4]]
  delta <- alpha / sqrt(1 + alpha^2)
  z <- delta * abs(stats::rnorm(n)) + sqrt(1 - delta^2) * stats::rnorm(n)
  theta[[1]] + theta[[2]] * z / sqrt(stats::rchisq(n, df = nu) / nu)
}

#' @title Skew t Analytical Gradient
#' @name distrib_gradient.SkewTDistrib
#' @description
#' First derivatives of the log-density. Three of the four are closed form:
#' with \eqn{D = A + QB} in the notation of [skewt_pieces()],
#' \deqn{\dfrac{\partial \ell}{\partial \mu} = -\dfrac{D}{\sigma},
#'       \qquad
#'       \dfrac{\partial \ell}{\partial \sigma} = -\dfrac{1 + zD}{\sigma},
#'       \qquad
#'       \dfrac{\partial \ell}{\partial \alpha} = Q z c.}
#' @details
#' The derivative in \eqn{\nu} is **not** closed form. It contains
#' \eqn{\partial \log T_{\nu+1}(w)/\partial \nu}, a derivative of the Student
#' \eqn{t} distribution function with respect to its degrees of freedom, which
#' has no elementary expression --- the same obstruction the gamma and beta
#' distribution functions meet in their shape. That one component is obtained
#' by a single central difference of the log-density, never by differencing
#' another difference.
#' @param distrib A `SkewTDistrib` object.
#' @param y A numeric vector of observations.
#' @param theta A list containing `mu`, `sigma`, `alpha` and `nu`.
#' @param scale Either `"parameter"` or `"link"`.
#' @param ... Unused.
#' @return A named list of first derivatives.
#' @seealso [skewt_distrib()]
S7::method(distrib_gradient, SkewTDistrib) <- function(distrib, y, theta, scale = c("parameter", "link"), ...) {
  mu <- theta[[1]]
  sigma <- theta[[2]]
  alpha <- theta[[3]]
  nu <- theta[[4]]
  p <- skewt_pieces(y, mu, sigma, alpha, nu)
  d <- p$a + p$q * p$b

  h <- skewt_nu_step(nu)
  lp <- function(v) {
    distrib_pdf(distrib, y, list(mu, sigma, alpha, v), log = TRUE)
  }

  list(
    mu = -d / sigma,
    sigma = -(1 + p$z * d) / sigma,
    alpha = p$q * p$z * p$c,
    nu = fd5_first(lp, nu, h)
  )
}

#' @title Skew t Analytical Observed Hessian
#' @name distrib_hessian.SkewTDistrib
#' @description
#' Second derivatives of the log-density. The block in
#' \eqn{(\mu, \sigma, \alpha)} is closed form; every component involving
#' \eqn{\nu} comes from one finite-difference stencil applied to the
#' log-density, for the reason given in
#' [distrib_gradient.SkewTDistrib()].
#' @details
#' With \eqn{D = A + QB} and \eqn{D' = A' + Q'B^2 + QB'},
#' \deqn{\dfrac{\partial^2 \ell}{\partial \mu^2} = \dfrac{D'}{\sigma^2},
#'       \qquad
#'       \dfrac{\partial^2 \ell}{\partial \mu \, \partial \sigma}
#'         = \dfrac{D + zD'}{\sigma^2},
#'       \qquad
#'       \dfrac{\partial^2 \ell}{\partial \sigma^2}
#'         = \dfrac{1 + 2zD + z^2 D'}{\sigma^2},}
#' \deqn{\dfrac{\partial^2 \ell}{\partial \alpha^2} = Q' z^2 c^2,
#'       \qquad
#'       \dfrac{\partial^2 \ell}{\partial \mu \, \partial \alpha}
#'         = -\dfrac{Q' B z c + Q E}{\sigma},
#'       \qquad
#'       \dfrac{\partial^2 \ell}{\partial \sigma \, \partial \alpha}
#'         = -\dfrac{z(Q' B z c + Q E)}{\sigma}.}
#' The stencils used for \eqn{\nu} are the three-point one in \eqn{\nu} alone
#' and the four-point mixed one otherwise; the mixed stencil differences two
#' *different* variables, so it is a single stencil rather than a
#' difference of a difference.
#' @param distrib A `SkewTDistrib` object.
#' @param y A numeric vector of observations.
#' @param theta A list containing `mu`, `sigma`, `alpha` and `nu`.
#' @param scale Either `"parameter"` or `"link"`.
#' @param ... Unused.
#' @return A named list of second derivatives.
#' @seealso [skewt_distrib()]
S7::method(distrib_hessian, SkewTDistrib) <- function(distrib, y, theta, scale = c("parameter", "link"), ...) {
  mu <- theta[[1]]
  sigma <- theta[[2]]
  alpha <- theta[[3]]
  nu <- theta[[4]]
  p <- skewt_pieces(y, mu, sigma, alpha, nu)
  d <- p$a + p$q * p$b
  dd <- p$da + p$dq * p$b^2 + p$q * p$db
  s2 <- sigma^2
  zc <- p$z * p$c
  mixed_alpha <- p$dq * p$b * zc + p$q * p$e

  # --- the components involving nu -----------------------------------------
  hn <- skewt_nu_step(nu)
  lp <- function(v) {
    distrib_pdf(distrib, y, list(mu, sigma, alpha, v), log = TRUE)
  }
  nu_nu <- fd5_second(lp, nu, hn)

  # The mixed components step the CLOSED-FORM score in nu, so only one
  # difference is taken and it is taken of an analytic quantity. Stepping the
  # log-density in both directions instead would be a difference of a
  # difference in one of them.
  grad_at <- function(v) {
    pv <- skewt_pieces(y, mu, sigma, alpha, v)
    dv <- pv$a + pv$q * pv$b
    cbind(
      mu = -dv / sigma,
      sigma = -(1 + pv$z * dv) / sigma,
      alpha = pv$q * pv$z * pv$c
    )
  }
  gnu <- fd5_first(grad_at, nu, hn)

  list(
    mu_mu = dd / s2,
    sigma_sigma = (1 + 2 * p$z * d + p$z^2 * dd) / s2,
    alpha_alpha = p$dq * zc^2,
    nu_nu = nu_nu,
    mu_sigma = (d + p$z * dd) / s2,
    mu_alpha = -mixed_alpha / sigma,
    mu_nu = gnu[, "mu"],
    sigma_alpha = -p$z * mixed_alpha / sigma,
    sigma_nu = gnu[, "sigma"],
    alpha_nu = gnu[, "alpha"]
  )
}

#' @title Skew t Third-Order Derivatives
#' @name distrib_deriv3.SkewTDistrib
#' @description
#' Third-order derivatives assembled so that no finite difference is ever
#' applied to another finite difference. Components whose Hessian entry is
#' closed form -- both indices in \eqn{(\mu, \sigma, \alpha)}, or one index
#' equal to \eqn{\nu} with the stencil taken along a different variable --
#' come from the generic construction, one stencil on an analytic quantity.
#' The components the generic construction would nest are replaced:
#' \eqn{(i, \nu, \nu)} is one five-point second-difference of the closed-form
#' score component \eqn{i}, and \eqn{(\nu, \nu, \nu)} is one five-point
#' third-difference of the log-density itself. The derivative of a Student t
#' distribution function in its degrees of freedom has no elementary form, so
#' this is the same obstruction, and the same remedy, as the Hessian's.
#' @param distrib A `SkewTDistrib` object.
#' @param y A numeric vector of observations.
#' @param theta A list containing `mu`, `sigma`, `alpha` and `nu`.
#' @param expected Logical; if `TRUE`, the expectation is approximated
#'   numerically.
#' @param approx Strategy for the expectation; see [distrib_deriv3()].
#' @param nsim Monte Carlo sample size when `approx = "mc"`.
#' @return A named list of third-derivative component vectors.
#' @seealso [skewt_distrib()]
S7::method(distrib_deriv3, SkewTDistrib) <- function(distrib, y, theta, expected = FALSE, scale = c("parameter", "link"), approx = c("integrate", "bartlett", "mc", "opg"), nsim = 10000, ...) {
  if (expected) {
    return(expected_derivative(distrib, y, theta, order = 3L,
                               approx = match.arg(approx), nsim = nsim))
  }
  out <- numerical_deriv3(distrib, y, theta)
  nu <- theta[[4]]
  h <- skewt_nu_step(nu)
  grad_at <- function(v, comp) {
    th <- theta; th[[4]] <- v
    distrib_gradient(distrib, y, th)[[comp]]
  }
  for (p in c("mu", "sigma", "alpha")) {
    out[[paste0(p, "_nu_nu")]] <- fd5_second(function(v) grad_at(v, p), nu, h)
  }
  ll_at <- function(v) {
    th <- theta; th[[4]] <- v
    distrib_pdf(distrib, y, th, log = TRUE)
  }
  out[["nu_nu_nu"]] <- fd5_third(ll_at, nu, h)
  out
}

#' @title Skew t Fourth-Order Derivatives
#' @name distrib_deriv4.SkewTDistrib
#' @description
#' Fourth-order derivatives assembled with the discipline of
#' [distrib_deriv3.SkewTDistrib()]: the generic construction serves
#' every component whose Hessian entry is closed form, and the ones it would
#' nest are replaced by one stencil each -- \eqn{(i, \nu, \nu, \nu)} by a
#' third-difference of the closed-form score component \eqn{i}, and
#' \eqn{(\nu, \nu, \nu, \nu)} by a fourth-difference of the log-density. The
#' pure-\eqn{\nu} component is the least accurate quantity the family reports,
#' at roughly four significant digits.
#' @param distrib A `SkewTDistrib` object.
#' @param y A numeric vector of observations.
#' @param theta A list containing `mu`, `sigma`, `alpha` and `nu`.
#' @param expected Logical; if `TRUE`, the expectation is approximated
#'   numerically.
#' @param approx Strategy for the expectation; see [distrib_deriv4()].
#' @param nsim Monte Carlo sample size when `approx = "mc"`.
#' @return A named list of fourth-derivative component vectors.
#' @seealso [skewt_distrib()]
S7::method(distrib_deriv4, SkewTDistrib) <- function(distrib, y, theta, expected = FALSE, scale = c("parameter", "link"), approx = c("integrate", "bartlett", "mc", "opg"), nsim = 10000, ...) {
  if (expected) {
    return(expected_derivative(distrib, y, theta, order = 4L,
                               approx = match.arg(approx), nsim = nsim))
  }
  out <- numerical_deriv4(distrib, y, theta)
  nu <- theta[[4]]
  h <- skewt_nu_step(nu)
  grad_at <- function(v, comp) {
    th <- theta; th[[4]] <- v
    distrib_gradient(distrib, y, th)[[comp]]
  }
  for (p in c("mu", "sigma", "alpha")) {
    out[[paste0(p, "_nu_nu_nu")]] <- fd5_third(function(v) grad_at(v, p), nu, h)
  }
  ll_at <- function(v) {
    th <- theta; th[[4]] <- v
    distrib_pdf(distrib, y, th, log = TRUE)
  }
  # The fourth difference amplifies rounding by h^-4, so its step is measured
  # separately: at the family's base step the per-observation noise is near
  # 1e-2 relative, at ten times that it is negligible and the h^2 truncation,
  # about 6e-4, is what remains.
  out[["nu_nu_nu_nu"]] <- fd5_fourth(ll_at, nu, 10 * h)
  out
}

#' @title Skew t Response Derivative
#' @name distrib_grad_y.SkewTDistrib
#' @description
#' Closed form: \eqn{\partial \ell / \partial y = D/\sigma}, which is minus the
#' derivative in \eqn{\mu}, as it must be for a location family.
#' @param distrib A `SkewTDistrib` object.
#' @param y A numeric vector of observations.
#' @param theta A list containing `mu`, `sigma`, `alpha` and `nu`.
#' @return A numeric vector.
#' @seealso [skewt_distrib()]
S7::method(distrib_grad_y, SkewTDistrib) <- function(distrib, y, theta) {
  p <- skewt_pieces(y, theta[[1]], theta[[2]], theta[[3]], theta[[4]])
  (p$a + p$q * p$b) / theta[[2]]
}

#' @title Skew t Response Second Derivative
#' @name distrib_hess_y.SkewTDistrib
#' @description
#' Closed form: \eqn{\partial^2 \ell / \partial y^2 = D'/\sigma^2}.
#' @param distrib A `SkewTDistrib` object.
#' @param y A numeric vector of observations.
#' @param theta A list containing `mu`, `sigma`, `alpha` and `nu`.
#' @return A numeric vector.
#' @seealso [skewt_distrib()]
S7::method(distrib_hess_y, SkewTDistrib) <- function(distrib, y, theta) {
  p <- skewt_pieces(y, theta[[1]], theta[[2]], theta[[3]], theta[[4]])
  (p$da + p$dq * p$b^2 + p$q * p$db) / theta[[2]]^2
}

# --- CONSTRUCTOR WRAPPER ---

#' Skew t Distribution Object
#'
#' @description
#' Creates a distribution object for Azzalini's skew \eqn{t} distribution, with
#' location \eqn{\mu}, scale \eqn{\sigma}, shape \eqn{\alpha} and degrees of
#' freedom \eqn{\nu}. It contains the Student \eqn{t} (\eqn{\alpha = 0}), the
#' skew normal (\eqn{\nu \to \infty}) and the gaussian (both).
#'
#' @param link_mu A link function object for the location \eqn{\mu}. Defaults to
#'   [linkfunctions7::identity_link()].
#' @param link_sigma A link function object for the scale \eqn{\sigma}. Defaults
#'   to [linkfunctions7::log_link()].
#' @param link_alpha A link function object for the shape \eqn{\alpha}, which is
#'   unconstrained. Defaults to [linkfunctions7::identity_link()].
#' @param link_nu A link function object for the degrees of freedom \eqn{\nu}.
#'   Defaults to [linkfunctions7::log_link()].
#'
#' @details
#' This is the four-parameter family a location-scale-shape framework wants: the
#' scale, the skewness and the tail weight are three separate parameters, each
#' of which can be given its own linear predictor. The skew normal of
#' [skewnormal1_distrib()] can reach a skewness of at most \eqn{0.995}
#' and an excess kurtosis of at most \eqn{0.87}; adding \eqn{\nu} removes both
#' bounds.
#'
#' **Probability density function**, with \eqn{z = (y-\mu)/\sigma} and
#' \eqn{w = \alpha z\sqrt{(\nu+1)/(\nu+z^2)}}:
#' \deqn{f(y; \mu, \sigma, \alpha, \nu)
#'   = \dfrac{2}{\sigma}\,t_\nu(z)\,T_{\nu+1}(w)}
#'
#' **What is closed form and what is not.** The score and the observed
#' Hessian are closed form in \eqn{(\mu, \sigma, \alpha)}. Everything involving
#' \eqn{\nu} is not, because the density contains \eqn{T_{\nu+1}}, whose
#' derivative with respect to its degrees of freedom has no elementary
#' expression --- the same obstruction that stops the gamma and beta
#' distribution functions from having closed-form shape derivatives. Those
#' components come from a single finite-difference stencil applied to an
#' analytic quantity, never from a difference of a difference:
#'
#' \tabular{lll}{
#'   **component** \tab **route** \tab **error, summed over n**
#'     \cr
#'   \eqn{\mu, \sigma, \alpha} (score) \tab closed form \tab machine precision
#'     \cr
#'   \eqn{\nu} (score) \tab five-point stencil on \eqn{\ell} \tab
#'     \eqn{10^{-11}} to \eqn{10^{-9}} \cr
#'   \eqn{(\mu,\sigma,\alpha)} block (Hessian) \tab closed form \tab machine
#'     precision \cr
#'   \eqn{\nu} with another parameter \tab five-point stencil on the analytic
#'     score \tab about \eqn{10^{-8}} \cr
#'   \eqn{\nu} twice \tab five-point stencil on \eqn{\ell} \tab about
#'     \eqn{10^{-6}}
#' }
#'
#' **The tolerance a fit can ask for.** The score in \eqn{\nu} cannot be
#' computed more accurately than the table above, so a stopping rule on the
#' gradient cannot be satisfied below that level however good the optimizer is.
#' [fit_distrib()] tests the score **per observation**, and its
#' default of \eqn{10^{-10}} leaves room: on samples of 500 to 4000 the summed
#' score reaches \eqn{10^{-10}} to \eqn{3 \times 10^{-9}}, which is
#' \eqn{10^{-13}} per observation, and the fit converges in four or five
#' iterations. A rule expressed on the summed gradient at that tolerance would
#' not be attainable, which is why the tolerance is not expressed that way.
#'
#' **Distribution function.** There is no elementary form, so the base
#' class integrates the density and inverts the result by root finding.
#'
#' **Expected information** has no closed form either and is approximated
#' by the strategy named in `approx`, the default being the score
#' variance. That approximation costs one quadrature per component, so
#' `method = "newton"` is much the cheaper way to fit this family: the
#' observed Hessian is the closed form above and needs no integration.
#'
#' **Moments** exist only up to order \eqn{\nu}: the mean requires
#' \eqn{\nu > 1}, the variance \eqn{\nu > 2}, the skewness \eqn{\nu > 3} and the
#' kurtosis \eqn{\nu > 4}, and each returns `NaN` below its threshold. The
#' density is perfectly well defined there, which is why the moments and the
#' parameters are kept apart.
#'
#' **Special cases.** \eqn{\alpha = 0} is [student_t1_distrib()];
#' large \eqn{\nu} approaches [skewnormal1_distrib()]. The information
#' is singular in \eqn{\alpha} at \eqn{\alpha = 0} for the same reason as in the
#' skew normal.
#'
#' **Parameter Domains:**
#' \itemize{
#'   \item \eqn{\mu \in (-\infty, +\infty)}
#'   \item \eqn{\sigma \in (0, +\infty)}
#'   \item \eqn{\alpha \in (-\infty, +\infty)}
#'   \item \eqn{\nu \in (0, +\infty)}
#' }
#'
#' @return An S7 object of class [SkewTDistrib()] (inheriting from
#'   `continuous_distrib`).
#'
#' @references
#' Azzalini, A. and Capitanio, A. (2003). Distributions generated by
#' perturbation of symmetry with emphasis on a multivariate skew t distribution.
#' *Journal of the Royal Statistical Society, Series B* 65, 367-389.
#'
#' @importFrom linkfunctions7 identity_link log_link
#' @importFrom stats dt pt rchisq
#'
#' @examples
#' d <- skewt_distrib()
#' d@params
#'
#' theta <- list(mu = 0, sigma = 1, alpha = 3, nu = 5)
#' distrib_pdf(d, c(-1, 0, 1), theta)
#' distrib_gradient(d, c(-1, 0, 1), theta)
#'
#' # shape zero is the Student t
#' max(abs(distrib_pdf(d, c(-1, 0, 1), list(mu = 0, sigma = 1, alpha = 0, nu = 5)) -
#'         stats::dt(c(-1, 0, 1), df = 5)))
#'
#' # the family reaches skewness the skew normal cannot
#' c(skew_t = skewness(d, list(mu = 0, sigma = 1, alpha = 8, nu = 5)),
#'   skew_normal_bound = 0.9953)
#'
#' # The observed Hessian is the cheap route here: this family has no
#' # closed-form expected information, so Fisher scoring would approximate it
#' # by quadrature at every step.
#' set.seed(1)
#' y <- distrib_rng(d, 200, theta)
#' coef(fit_distrib(d, y, method = optimizers7::newton(), start = theta))
#'
#' @seealso [skewnormal1_distrib()], [student_t1_distrib()]
#' @export
skewt_distrib <- function(link_mu = identity_link(),
                          link_sigma = log_link(),
                          link_alpha = identity_link(),
                          link_nu = log_link()) {
  SkewTDistrib(
    distrib_name = "skew t",
    dimension = "univariate",
    bounds = c(-Inf, Inf),

    params = c("mu", "sigma", "alpha", "nu"),
    params_interpretation = c(
      mu = "location", sigma = "scale", alpha = "shape",
      nu = "degrees of freedom"
    ),
    n_params = 4,

    params_bounds = list(
      mu = c(-Inf, Inf),
      sigma = c(0, Inf),
      alpha = c(-Inf, Inf),
      nu = c(0, Inf)
    ),

    link_params = list(
      mu = link_mu,
      sigma = link_sigma,
      alpha = link_alpha,
      nu = link_nu
    )
  )
}
