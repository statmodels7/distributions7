#' @include distrib.R generics.R skewnormal1_distrib.R expected_loc_scale.R cdf_compiled.R
NULL

#' @title Skew t Distribution Class
#' @name SkewTDistrib
#'
#' @description
#' The S7 class of Azzalini's skew \eqn{t}: a Student \eqn{t} carrying a shape
#' parameter that tilts it, so that the tail weight and the asymmetry are
#' modeled by two parameters instead of one. With \eqn{z = (y-\mu)/\sigma} and
#' \eqn{w = \alpha z\sqrt{(\nu+1)/(\nu+z^2)}} the density is
#' \eqn{2 t_\nu(z)T_{\nu+1}(w)/\sigma}.
#'
#' It contains three families as limits. At \eqn{\alpha = 0} it is the Student
#' \eqn{t}; as \eqn{\nu \to \infty} it is the skew normal, approached at
#' \eqn{O(1/\nu)}; and with both it is the Gaussian. Its reason to exist is the
#' skewness: the skew normal cannot pass 0.9953, and this family reaches 2.05
#' at \eqn{\nu = 6} and 4.00 at \eqn{\nu = 4}.
#'
#' Build one with [skewt_distrib()], which supplies the four link functions.
#' This page documents the raw S7 constructor, which validates none of the
#' relationships between its properties.
#'
#' @inheritParams distrib
#'
#' @return An S7 object of class `SkewTDistrib`, inheriting from
#'   `continuous_distrib` and from `distrib`. For an object built by
#'   [skewt_distrib()] the properties hold `"skew t"`, `"univariate"`,
#'   `c(-Inf, Inf)`, `c("mu", "sigma", "alpha", "nu")`, the interpretations
#'   `c(mu = "location", sigma = "scale", alpha = "shape", nu = "degrees of
#'   freedom")`, `4`, and the domains \eqn{(-\infty,\infty)},
#'   \eqn{(0,\infty)}, \eqn{(-\infty,\infty)} and \eqn{(0,\infty)}.
#'
#' @section Methods:
#' Registered in this file:
#'   [`distrib_pdf()`][distrib_pdf.SkewTDistrib],
#'   [`distrib_rng()`][distrib_rng.SkewTDistrib],
#'   [`distrib_gradient()`][distrib_gradient.SkewTDistrib],
#'   [`distrib_hessian()`][distrib_hessian.SkewTDistrib],
#'   [`distrib_deriv3()`][distrib_deriv3.SkewTDistrib],
#'   [`distrib_deriv4()`][distrib_deriv4.SkewTDistrib],
#'   [`distrib_deriv5()`][distrib_deriv5.SkewTDistrib],
#'   [`distrib_grad_y()`][distrib_grad_y.SkewTDistrib],
#'   [`distrib_hess_y()`][distrib_hess_y.SkewTDistrib], and the expected
#'   information with its two derivatives,
#'   [`distrib_expected_hessian()`][distrib_expected_hessian.SkewTDistrib].
#'
#' Registered elsewhere: all four moments in `moments.R`
#' ([`mean()`][mean.SkewTDistrib], [`variance()`][variance.SkewTDistrib],
#' [`skewness()`][skewness.SkewTDistrib],
#' [`kurtosis()`][kurtosis.SkewTDistrib]); the mixed response-parameter
#' derivative [`distrib_cross_y()`][distrib_cross_y] in
#' `cross_derivatives_families.R`; and
#' [`distrib_grad_cdf()`][distrib_grad_cdf] in `cdf_derivatives_families.R`.
#'
#' The **distribution function** and the **quantile function** come from
#' [continuous_distrib()], by quadrature and by root finding on it. The
#' **expected information** has no elementary form and is computed by one
#' quadrature per distinct \eqn{(\alpha, \nu)}; see
#' [distrib_expected_hessian.SkewTDistrib()].
#'
#' @section How the derivatives are computed:
#' Every derivative, to order five, comes from a compiled kernel in closed
#' form. The density carries \eqn{T_{\nu+1}}, and the derivatives of a Student
#' \eqn{t} distribution function in its degrees of freedom have no elementary
#' expression; they are integrals of the derivatives of the \eqn{t} density and
#' are taken by quadrature. See [distrib_gradient.SkewTDistrib()].
#'
#' @seealso [skewt_distrib()] to build one;
#'   [skewnormal1_distrib()] for the \eqn{\nu \to \infty} limit;
#'   [student_t1_distrib()] for the \eqn{\alpha = 0} case;
#'   [skewt_pieces()] for the scalar functions the derivatives are built from.
#'
#' @examples
#' d <- skewt_distrib()
#' S7::S7_inherits(d, continuous_distrib)
#'
#' d@params
#' d@params_interpretation
#'
#' # Two parameters for two departures from the Gaussian, which is why the
#' # family exists.
#' vapply(d@link_params, function(l) l@link_name, character(1))
#'
#' # It passes the skew normal's ceiling of 0.9953.
#' vapply(c(4, 6, 20, 1e6),
#'        function(v) skewness(d, list(mu = 0, sigma = 1, alpha = 50, nu = v)), 0)
SkewTDistrib <- S7::new_class("SkewTDistrib", parent = continuous_distrib)

# --- S7 METHODS IMPLEMENTATION ---

#' @title The Pieces a Skew t Evaluates From
#'
#' @description
#' Assembles the standardized variable, the argument of the tilting
#' distribution function and the six scalar functions every closed-form
#' derivative of the log-density is a combination of. Every method in this
#' family calls it once and then writes its own formula in terms of the result.
#'
#' @details
#' With \eqn{z = (y-\mu)/\sigma}, \eqn{m = \nu + 1}, \eqn{s = \nu + z^2} and
#' \eqn{c = \sqrt{m/s}}, the tilting argument is \eqn{w = \alpha z c}. The
#' functions returned are
#' \deqn{A = \dfrac{\partial}{\partial z}\log t_\nu(z) = -\dfrac{m z}{s},
#'       \qquad
#'       A' = -\dfrac{m(\nu - z^2)}{s^2},}
#' \deqn{E = \dfrac{\nu\sqrt{m}}{s^{3/2}}, \qquad
#'       B = \dfrac{\partial w}{\partial z} = \alpha E, \qquad
#'       B' = -\dfrac{3\alpha\nu\sqrt{m}\,z}{s^{5/2}},}
#' and \eqn{Q = t_m(w)/T_m(w)} with
#' \eqn{Q' = Q\{-(m+1)w/(m + w^2) - Q\}}, the last from differentiating the
#' quotient.
#'
#' \eqn{Q} is formed as `exp(dt(log = TRUE) - pt(log.p = TRUE))` because both
#' factors underflow together in the far left tail while the ratio stays
#' finite. It matters as the degrees of freedom grow and the \eqn{t} tail
#' approaches the Gaussian's: measured at \eqn{w = -60} with
#' \eqn{m = 2000}, the log route returns 21.4345 and `dt(w, m)/pt(w, m)`
#' returns `NaN`.
#'
#' Nothing here differentiates in \eqn{\nu}. The derivative kernels evaluate
#' the same quantities in bounded variables; see
#' [distrib_gradient.SkewTDistrib()].
#'
#' @param y A numeric vector of observations.
#' @param mu,sigma,alpha,nu The four parameters, numeric vectors of length 1 or
#'   of the length of `y`. Nothing is validated: `sigma` and `nu` must be
#'   strictly positive.
#'
#' @return A named list of ten numeric vectors, each of the length of the
#'   recycled inputs: `z` the standardized variable, `w` the tilting argument,
#'   `c` the factor relating them, `a` and `da` for \eqn{A} and \eqn{A'}, `e`
#'   for \eqn{E}, `b` and `db` for \eqn{B} and \eqn{B'}, and `q` and `dq` for
#'   \eqn{Q} and \eqn{Q'}.
#'
#' @section Notation:
#' \eqn{t_\nu} and \eqn{T_\nu} are the standard Student \eqn{t} density and
#' distribution function on \eqn{\nu} degrees of freedom, \eqn{\mu} the
#' location, \eqn{\sigma} the scale, \eqn{\alpha} the shape and \eqn{\nu} the
#' degrees of freedom.
#'
#' @seealso [distrib_gradient.SkewTDistrib()], which writes the score in these
#'   terms, and [skewt_distrib()] for the family.
#'
#' @examples
#' p <- distributions7:::skewt_pieces(c(-1.5, 0.4, 2.1), 0, 1, 3, 6)
#' names(p)
#'
#' # The score in the location is -(A + QB)/sigma.
#' d <- skewt_distrib()
#' all.equal(-(p$a + p$q * p$b) / 1,
#'           distrib_gradient(d, c(-1.5, 0.4, 2.1),
#'                            list(mu = 0, sigma = 1, alpha = 3, nu = 6))$mu)
#'
#' # Q survives a tail where the direct quotient does not.
#' c(log_route = exp(dt(-60, df = 2000, log = TRUE) -
#'                   pt(-60, df = 2000, log.p = TRUE)),
#'   direct = dt(-60, df = 2000) / pt(-60, df = 2000))
#'
#' @keywords internal
skewt_pieces <- function(y, mu, sigma, alpha, nu) {
  z <- (y - mu) / sigma
  m <- nu + 1
  s <- nu + z^2
  cc <- sqrt(m / s)
  w <- alpha * z * cc
  e <- nu * sqrt(m) / s^1.5
  a <- -m * z / s
  da <- -m * (nu - z^2) / s^2
  db <- -3 * alpha * nu * sqrt(m) * z / s^2.5
  # s overflows near |z| = 1e154; far out the pieces are written in
  # rs = sqrt(s), z/rs and nu/s, as the registry's skewt_pieces() writes them
  big <- abs(z) > 1e100
  if (any(big)) {
    n <- length(z)
    rs <- abs(z) * sqrt(1 + nu / (z * z))
    az <- z / rs
    irs <- 1 / rs
    b2 <- nu * irs * irs
    rm <- sqrt(m)
    pick <- function(old, new) { old <- rep_len(old, n); old[big] <- rep_len(new, n)[big]; old }
    cc <- pick(cc, rm * irs)
    w <- pick(w, alpha * rm * az)
    a <- pick(a, -m * az * irs)
    da <- pick(da, -m * (b2 - az * az) * irs * irs)
    e <- pick(e, rm * b2 * irs)
    db <- pick(db, -3 * alpha * rm * b2 * az * irs * irs)
  }
  # Q is formed on the log scale for the same reason numericals7::mills_ratio() is: both the
  # density and the distribution function underflow in the far left tail while
  # their ratio stays finite.
  q <- exp(stats::dt(w, df = m, log = TRUE) - stats::pt(w, df = m, log.p = TRUE))
  list(
    z = z, w = w, c = cc,
    a = a,
    da = da,
    e = e,
    b = alpha * e,
    db = db,
    q = q,
    dq = q * (-(m + 1) * w / (m + w^2) - q)
  )
}

#' @title Skew t Density
#' @name distrib_pdf.SkewTDistrib
#'
#' @description
#' Computes the skew \eqn{t} density, with \eqn{z = (y-\mu)/\sigma} and
#' \eqn{w = \alpha z\sqrt{(\nu+1)/(\nu+z^2)}}:
#' \deqn{f(y; \mu, \sigma, \alpha, \nu) = \dfrac{2}{\sigma}\,
#'       t_\nu(z)\,T_{\nu+1}(w),}
#' with \eqn{t_\nu} the standard Student \eqn{t} density and \eqn{T_{\nu+1}}
#' its distribution function on **one more** degree of freedom. The extra
#' degree of freedom is deliberate: without it the tilted density would not
#' integrate to one.
#'
#' The two logarithms are taken separately and added, `dt(log = TRUE)` beside
#' `pt(log.p = TRUE)`, so the light tail of the skewed side returns a large
#' negative number instead of `-Inf`.
#'
#' @param distrib A `SkewTDistrib` object, from [skewt_distrib()].
#' @param y A numeric vector of observations, anywhere on the real line.
#' @param theta A named list with components `mu`, `sigma`, `alpha` and `nu`,
#'   each a numeric vector of length 1 or of the length of `y`. `sigma` and
#'   `nu` must be strictly positive; `mu` and `alpha` may take any finite
#'   value.
#' @param log Logical of length 1. When `TRUE` the log-density is returned.
#'   Defaults to `FALSE`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of densities, of the length of the recycled inputs.
#'
#' @section Notation:
#' \eqn{\mu} is the location, \eqn{\sigma > 0} the scale, \eqn{\alpha} the
#' shape and \eqn{\nu > 0} the degrees of freedom. Neither \eqn{\mu} nor
#' \eqn{\sigma} is a moment.
#'
#' @seealso [distrib_gradient.SkewTDistrib()] for the score,
#'   [skewnormal1_distrib()] for the \eqn{\nu \to \infty} limit,
#'   [student_t1_distrib()] for the \eqn{\alpha = 0} case, and [distrib_pdf()]
#'   for the generic.
#'
#' @examples
#' d <- skewt_distrib()
#' y <- c(-1.5, -0.3, 0.4, 2.1)
#' th <- list(mu = 0, sigma = 1, alpha = 3, nu = 6)
#'
#' # The formula written out.
#' w <- 3 * y * sqrt(7 / (6 + y^2))
#' all.equal(distrib_pdf(d, y, th), 2 * dt(y, 6) * pt(w, 7))
#'
#' # It integrates to one.
#' integrate(function(v) distrib_pdf(d, v, th), -Inf, Inf)$value
#'
#' # Shape zero is the Student t.
#' all.equal(distrib_pdf(d, y, list(mu = 0, sigma = 1, alpha = 0, nu = 6)),
#'           dt(y, 6))
#'
#' # Large degrees of freedom give the skew normal, at O(1/nu).
#' sn <- skewnormal1_distrib()
#' vapply(c(1e2, 1e4, 1e6), function(v)
#'   max(abs(distrib_pdf(d, y, list(mu = 0, sigma = 1, alpha = 3, nu = v)) -
#'           distrib_pdf(sn, y, list(mu = 0, sigma = 1, alpha = 3)))), 0)
S7::method(distrib_pdf, SkewTDistrib) <- function(distrib, y, theta, log = FALSE, ...) {
  mu <- theta[[1]]
  sigma <- theta[[2]]
  alpha <- theta[[3]]
  nu <- theta[[4]]
  z <- (y - mu) / sigma
  w <- alpha * z * sqrt((nu + 1) / (nu + z^2))
  # z^2 overflows near |z| = 1e154; far out w is formed from rs = sqrt(nu + z^2)
  # without it, as the registry's skewt_logpdf() forms it
  big <- abs(z) > 1e100
  if (any(big)) {
    rs <- abs(z) * sqrt(1 + nu / (z * z))
    wb <- alpha * sqrt(nu + 1) * (z / rs)
    w[big] <- rep_len(wb, length(w))[big]
  }
  log_d <- log(2) - log(sigma) + stats::dt(z, df = nu, log = TRUE) +
    stats::pt(w, df = nu + 1, log.p = TRUE)
  if (log) log_d else exp(log_d)
}

#' @title Skew t Random Generation
#' @name distrib_rng.SkewTDistrib
#'
#' @description
#' Draws from the skew \eqn{t} exactly, from its scale-mixture representation:
#' with \eqn{Z} standard skew normal of shape \eqn{\alpha} and
#' \eqn{V \sim \chi^2_\nu} independent of it,
#' \deqn{Y = \mu + \sigma\,\dfrac{Z}{\sqrt{V/\nu}}.}
#' \eqn{Z} itself is drawn from [distrib_rng.SkewNormal1Distrib()]'s
#' representation, so the whole draw is three `rnorm`/`rchisq` calls and no
#' inversion or rejection.
#'
#' The representation reads off both limits. As \eqn{\nu \to \infty} the mixing
#' factor tends to one and \eqn{Y} is skew normal; at \eqn{\alpha = 0} the
#' numerator is Gaussian and \eqn{Y} is Student \eqn{t}.
#'
#' @param distrib A `SkewTDistrib` object, from [skewt_distrib()].
#' @param n A single positive integer, the number of draws.
#' @param theta A named list with components `mu`, `sigma`, `alpha` and `nu`,
#'   each a numeric vector of length 1 or of length `n`; a component of length
#'   1 is recycled.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of `n` draws.
#'
#' @section Notation:
#' \eqn{\nu} is the degrees of freedom of the mixing chi-squared and
#' \eqn{\delta = \alpha/\sqrt{1+\alpha^2}} the weight of the half-normal
#' component of \eqn{Z}.
#'
#' @seealso [distrib_pdf.SkewTDistrib()] for the density the draws follow,
#'   [distrib_rng.SkewNormal1Distrib()] for the inner representation, and
#'   [distrib_rng()] for the generic.
#'
#' @examples
#' d <- skewt_distrib()
#' th <- list(mu = 1, sigma = 2, alpha = 3, nu = 8)
#'
#' set.seed(21)
#' x <- distrib_rng(d, 1e5, th)
#'
#' # Two moments against their closed forms. Both exist here; at nu <= 2 the
#' # variance does not.
#' rbind(sample = c(mean(x), var(x)),
#'       theory = c(mean(d, th), variance(d, th)))
#'
#' # The mixture, written out at the same seed.
#' set.seed(21)
#' delta <- 3 / sqrt(1 + 9)
#' z <- delta * abs(rnorm(1e5)) + sqrt(1 - delta^2) * rnorm(1e5)
#' all.equal(x, 1 + 2 * z / sqrt(rchisq(1e5, df = 8) / 8))
S7::method(distrib_rng, SkewTDistrib) <- function(distrib, n, theta, ...) {
  alpha <- theta[[3]]
  nu <- theta[[4]]
  delta <- alpha / sqrt(1 + alpha^2)
  z <- delta * abs(stats::rnorm(n)) + sqrt(1 - delta^2) * stats::rnorm(n)
  theta[[1]] + theta[[2]] * z / sqrt(stats::rchisq(n, df = nu) / nu)
}

#' @title Skew t Score
#' @name distrib_gradient.SkewTDistrib
#'
#' @description
#' Computes the four first derivatives of the log-density. Three are closed
#' form: with \eqn{D = A + QB} in the notation of [skewt_pieces()],
#' \deqn{\dfrac{\partial \ell}{\partial \mu} = -\dfrac{D}{\sigma},
#'       \qquad
#'       \dfrac{\partial \ell}{\partial \sigma} = -\dfrac{1 + zD}{\sigma},
#'       \qquad
#'       \dfrac{\partial \ell}{\partial \alpha} = Q z c.}
#'
#' The fourth contains \eqn{\partial_m \log T_m(w)}, \eqn{m = \nu + 1}, the
#' derivative of a Student \eqn{t} distribution function in its degrees of
#' freedom, which has no elementary expression and is taken as an integral.
#'
#' @details
#' # The compiled kernel
#'
#' Every order of this family, from the score to the fifth derivatives, comes
#' from one construction in compiled code. With \eqn{r = \sqrt{\nu + z^2}},
#' \eqn{a = z/r} and \eqn{m = \nu + 1},
#' \deqn{\ell = \log 2 - \log\sigma + C(\nu) - \frac{\nu + 1}{2}
#'   \log\Big(1 + \frac{z^2}{\nu}\Big) + \log T_m(w), \qquad
#'   w = \alpha\sqrt{m}\,a,}
#' with \eqn{C(m) = \log\Gamma((m+1)/2) - \log\Gamma(m/2) - \log(m\pi)/2}.
#' Every derivative in \eqn{(\mu, \sigma, \alpha, \nu)} is written out by the
#' chain rule in the bounded variables \eqn{a}, \eqn{1/r}, \eqn{1/\sigma},
#' \eqn{\alpha}, \eqn{\nu} and \eqn{\sqrt m}, never in \eqn{z} itself, and
#' reduced with \eqn{a^2 = 1 - \nu/r^2}, so that the terms which cancel as
#' \eqn{|z|} grows cancel in the algebra. The expressions are generated by a
#' script kept with the package sources.
#'
#' The derivatives of \eqn{C} are \eqn{2^{-k} D^{(k)}(m/2)}, with
#' \eqn{D(x) = \log\Gamma(x + 1/2) - \log\Gamma(x) - \log(x)/2}, shifted to
#' \eqn{x \ge 20} by its recurrence and summed there by its asymptotic series.
#' The derivatives of \eqn{\log T_m(w)} are polynomials in
#' \eqn{q = t_m(w)/T_m(w)}, in \eqn{\Delta = \partial_w \log t_m(w) - q}, in the
#' partial derivatives of \eqn{\log t_m(w)} and in
#' \eqn{\delta_j = \partial_m^j T_m(w)/T_m(w) - \partial_m^j t_m(w)/t_m(w)}.
#' Since \eqn{T_m(0) = 1/2} for every \eqn{m},
#' \eqn{\partial_m^j T_m(w) = \int_0^w \partial_m^j t_m(u)\,du}. The integral is
#' taken by 20-point Gauss-Legendre for \eqn{|w| \le 2}; beyond, it is minus the
#' integral over the tail on the far side of \eqn{w}, taken by 12-point
#' Gauss-Legendre on panels of a variable in which the tail of the density
#' decays exponentially. On the left the integrand is
#' \eqn{t_m(u)\{P_j(u) - P_j(w)\}}, \eqn{P_j = \partial_m^j t_m/t_m}, whose
#' integral is \eqn{T_m(w)\,\delta_j}.
#'
#' Against `mpmath.diff` at 50 digits, on the closed form through the
#' incomplete beta function, every order from one to five agrees to
#' \eqn{2\times10^{-13}} relative to its largest component at five points
#' with \eqn{\nu} from 0.8 to 60.
#'
#' @param distrib A `SkewTDistrib` object, from [skewt_distrib()].
#' @param y A numeric vector of observations.
#' @param theta A named list with components `mu`, `sigma`, `alpha` and `nu`,
#'   each a numeric vector of length 1 or of the length of `y`. `sigma` and
#'   `nu` must be strictly positive.
#' @param scale Either `"parameter"`, the default, or `"link"`. The
#'   transformation to the link scale is applied in the generic's body, so this
#'   method always returns the parameter scale.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#' @param threads A single positive integer, how many threads the kernel may
#'   use.
#'
#' @return A named list of four numeric vectors, `mu`, `sigma`, `alpha` and
#'   `nu`, each of the length of the recycled inputs.
#'
#' @section Notation:
#' \eqn{z = (y-\mu)/\sigma}, \eqn{c = \sqrt{(\nu+1)/(\nu+z^2)}},
#' \eqn{w = \alpha z c}, and \eqn{A}, \eqn{B}, \eqn{Q} are as
#' [skewt_pieces()] defines them.
#'
#' @seealso [skewt_pieces()] for the scalar functions,
#'   [distrib_hessian.SkewTDistrib()] for the second derivatives, and
#'   [distrib_gradient()] for the generic.
#'
#' @examples
#' d <- skewt_distrib()
#' y <- c(-1.5, -0.3, 0.4, 2.1)
#' th <- list(mu = 0, sigma = 1, alpha = 3, nu = 6)
#' g <- distrib_gradient(d, y, th)
#'
#' # The location component written out in skewt_pieces()' terms.
#' p <- distributions7:::skewt_pieces(y, 0, 1, 3, 6)
#' all.equal(g$mu, -(p$a + p$q * p$b))
#'
#' # All four against numerical differentiation of the log-likelihood.
#' f <- function(v) sum(distrib_pdf(d, y, list(mu = v[1], sigma = v[2],
#'                                             alpha = v[3], nu = v[4]),
#'                                  log = TRUE))
#' rbind(analytic = vapply(g, sum, 0),
#'       numeric = numDeriv::grad(f, c(0, 1, 3, 6)))
#'
#' # At shape zero the score in alpha is not zero: the tilting factor is at
#' # its inflection, so alpha is still identified.
#' distrib_gradient(d, y, list(mu = 0, sigma = 1, alpha = 0, nu = 6))$alpha
S7::method(distrib_gradient, SkewTDistrib) <- function(distrib, y, theta, scale = c("parameter", "link"), ..., threads = 1L) {
  skewt_gradient_cpp(y, theta[[1]], theta[[2]], theta[[3]], theta[[4]], threads)
}

#' @title Skew t Observed Hessian
#' @name distrib_hessian.SkewTDistrib
#'
#' @description
#' Computes the ten second derivatives of the log-density with the compiled
#' kernel described on [distrib_gradient.SkewTDistrib()]. The block in
#' \eqn{(\mu, \sigma, \alpha)} has the closed forms below; the components
#' involving \eqn{\nu} carry the derivatives of \eqn{\log T_{\nu+1}(w)} in the
#' degrees of freedom, which are integrals.
#'
#' @details
#' # The closed-form block
#'
#' With \eqn{D = A + QB} and \eqn{D' = A' + Q'B^2 + QB'} in the notation of
#' [skewt_pieces()],
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
#'
#' @param distrib A `SkewTDistrib` object, from [skewt_distrib()].
#' @param y A numeric vector of observations.
#' @param theta A named list with components `mu`, `sigma`, `alpha` and `nu`,
#'   each a numeric vector of length 1 or of the length of `y`. `sigma` and
#'   `nu` must be strictly positive.
#' @param scale Either `"parameter"`, the default, or `"link"`. The
#'   transformation is applied in the generic's body, so this method always
#'   returns the parameter scale.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#' @param threads A single positive integer, how many threads the kernel may
#'   use.
#'
#' @return A named list of ten numeric vectors in [hess_names()]'s order:
#'   `mu_mu`, `sigma_sigma`, `alpha_alpha`, `nu_nu`, then `mu_sigma`,
#'   `mu_alpha`, `mu_nu`, `sigma_alpha`, `sigma_nu`, `alpha_nu`.
#'
#' @section Notation:
#' \eqn{z = (y-\mu)/\sigma}, \eqn{c = \sqrt{(\nu+1)/(\nu+z^2)}}, and
#' \eqn{A}, \eqn{B}, \eqn{E}, \eqn{Q} are as [skewt_pieces()] defines them.
#'
#' @seealso [distrib_gradient.SkewTDistrib()] for the order below,
#'   [distrib_deriv3.SkewTDistrib()] for the order above, and
#'   [distrib_hessian()] for the generic.
#'
#' @examples
#' d <- skewt_distrib()
#' y <- c(-1.5, -0.3, 0.4, 2.1)
#' th <- list(mu = 0, sigma = 1, alpha = 3, nu = 6)
#' h <- distrib_hessian(d, y, th)
#' names(h)
#'
#' # Against numerical differentiation of the log-likelihood: one component
#' # of the closed-form block, and two involving nu.
#' f <- function(v) sum(distrib_pdf(d, y, list(mu = v[1], sigma = v[2],
#'                                             alpha = v[3], nu = v[4]),
#'                                  log = TRUE))
#' H <- numDeriv::hessian(f, c(0, 1, 3, 6))
#' rbind(analytic = c(sum(h$mu_mu), sum(h$mu_nu), sum(h$nu_nu)),
#'       numeric = c(H[1, 1], H[1, 4], H[4, 4]))
#'
#' # The curvature in the location turns positive far out, as a Student t's
#' # does and a Gaussian's does not.
#' range(distrib_hess_y(d, seq(-20, 20, by = 0.5), th))
S7::method(distrib_hessian, SkewTDistrib) <- function(distrib, y, theta, scale = c("parameter", "link"), ..., threads = 1L) {
  skewt_hessian_cpp(y, theta[[1]], theta[[2]], theta[[3]], theta[[4]], threads)
}

#' @title The Skew t Tower in the Location, Scale and Shape
#' @name skewt_msa_tower
#'
#' @description
#' Builds the table of \eqn{\partial^i_z\,\partial^c_\alpha \ell}, for every
#' pair with \eqn{c + i \le 4}. Together with [skewt_msa_component()] it gives
#' in closed form every derivative of the log-density that does not involve
#' \eqn{\nu}, by a derivation in \eqn{z} independent of the compiled kernels
#' of [distrib_gradient.SkewTDistrib()], which are tested against it.
#'
#' @details
#' # Why the block closes
#'
#' With \eqn{z = (y-\mu)/\sigma} and \eqn{u(z) = z\sqrt{(\nu+1)/(\nu+z^2)}},
#' the log-density is
#' \deqn{\ell = \log 2 - \log\sigma + g(z) + \Lambda(\alpha\,u(z)),}
#' where \eqn{g = \log t_\nu} and \eqn{\Lambda = \log T_{\nu+1}}. Neither
#' \eqn{\mu} nor \eqn{\sigma} appears anywhere but inside \eqn{z}, and
#' \eqn{\alpha} nowhere but inside \eqn{w = \alpha u}, so two observations
#' settle the whole block.
#'
#' The shape enters the argument of \eqn{\Lambda} linearly, so differentiating
#' in it never leaves the table:
#' \deqn{\Phi_c(z) \;=\; \dfrac{\partial^c \ell}{\partial\alpha^c}
#'   \;=\; u(z)^c\,\Lambda^{(c)}(w), \qquad c \ge 1 .}
#'
#' The location and the scale then act on a function of \eqn{z} alone, and for
#' any such \eqn{F},
#' \deqn{\dfrac{\partial^{a+b} F}{\partial\mu^a\,\partial\sigma^b}
#'       = \dfrac{(-1)^{a+b}}{\sigma^{a+b}}\;P_{a,b}(z), \qquad
#'       P_{a,b}(z) = \sum_{i=0}^{b}\binom{b}{i}
#'         \left[\prod_{k=i}^{b-1}(a+k)\right] z^i\,F^{(a+i)}(z),}
#' which follows by induction from \eqn{P_{a,0} = F^{(a)}} and
#' \eqn{P_{a,b+1} = (a+b)\,P_{a,b} + z\,P_{a,b}'}. The explicit
#' \eqn{-\log\sigma} contributes only to the pure-\eqn{\sigma} components,
#' where it adds \eqn{(-1)^b (b-1)!\,\sigma^{-b}}.
#'
#' # Nothing in it needs a difference
#'
#' \eqn{u} and \eqn{g} are rational in \eqn{z} up to one square root, so their
#' four derivatives are written out. \eqn{\Lambda^{(k)}} follows from the
#' Riccati recursion \eqn{Q' = Q(G-Q)} for \eqn{Q = t_{\nu+1}/T_{\nu+1}} and
#' \eqn{G = \partial_w \log t_{\nu+1}}, differentiated twice more, exactly as a
#' skew normal's inverse Mills ratio is handled. The two chain rules that
#' remain, \eqn{(u^c)^{(j)}} and \eqn{\partial_z^p \Lambda^{(c)}(\alpha u)},
#' are Faa di Bruno sums over the partial Bell polynomials of \eqn{u}.
#'
#' What is left outside the table is \eqn{\nu}, which enters the degrees of
#' freedom of \eqn{T_{\nu+1}} and therefore carries the obstruction
#' [distrib_gradient.SkewTDistrib()] records.
#'
#' @param y A numeric vector of observations.
#' @param mu,sigma,alpha,nu The four parameters, each of length 1 or of the
#'   length of `y`.
#'
#' @return A named list. `z` is the standardized residual, `rs` is
#'   \eqn{\sqrt{\nu + z^2}} and `a` is `z / rs`; the remaining fifteen
#'   elements are named `"c_i"` and hold
#'   \eqn{\partial_z^i \Phi_c(z)} times `rs^i`, which bounds them, for every
#'   \eqn{c + i \le 4}, with
#'   \eqn{\Phi_0 = \ell} up to the terms free of \eqn{z}, so that `"0_0"` is
#'   `NA_real_` and is never read.
#'
#' @seealso [skewt_msa_component()], which reads one derivative off the table,
#'   and [skewt_pieces()], whose `a`, `e` and `q` are this function's
#'   \eqn{g'}, \eqn{u'} and \eqn{Q}.
#'
#' @examples
#' tw <- distributions7:::skewt_msa_tower(c(-0.4, 1.2), 0, 1, 0.7, 8)
#' names(tw)
#'
#' @keywords internal
skewt_msa_tower <- function(y, mu, sigma, alpha, nu) {
  z <- (y - mu) / sigma
  m <- nu + 1
  rm <- sqrt(m)
  # s = nu + z^2 through rs = sqrt(s), formed so that z^2 never overflows,
  # with a = z/rs and b2 = nu/s, both bounded (a^2 + b2 = 1). Every
  # z-derivative of order j below is held times rs^j, which bounds it: the
  # Bell polynomials are homogeneous of weight j, so the table is the
  # derivatives times rs^i, and skewt_msa_component() divides the power
  # back out. Written in z and s, the fourth derivatives were NaN from |y|
  # of about 1e80.
  az <- abs(z)
  rs <- ifelse(az > 1, az * sqrt(1 + nu / (z * z)), sqrt(nu + z * z))
  a <- z / rs
  b2 <- nu / (rs * rs)

  # u = z sqrt((nu+1)/(nu+z^2)) and its four z-derivatives, written out.
  u <- rm * a
  u1 <- rm * b2
  u2 <- -3 * rm * b2 * a
  u3 <- -3 * rm * b2 * (b2 - 4 * a^2)
  u4 <- 15 * rm * b2 * a * (3 * b2 - 4 * a^2)

  # g = log t_nu and its four z-derivatives.
  gz <- list(
    -m * a,
    -m * (b2 - a^2),
    2 * m * a * (3 * b2 - a^2),
    6 * m * (b2^2 - 6 * b2 * a^2 + a^4)
  )

  # Lam^(k)(w) = Q^(k-1)(w) by the Riccati recursion, Q formed on the log scale
  # for the reason skewt_pieces() gives.
  w <- alpha * u
  q <- exp(stats::dt(w, df = m, log = TRUE) - stats::pt(w, df = m, log.p = TRUE))
  mw <- m + w^2
  g0 <- -(m + 1) * w / mw
  g1 <- -(m + 1) * (m - w^2) / mw^2
  g2 <- 2 * (m + 1) * w * (3 * m - w^2) / mw^3
  q1 <- q * (g0 - q)
  q2 <- q1 * (g0 - q) + q * (g1 - q1)
  q3 <- q2 * (g0 - q) + 2 * q1 * (g1 - q1) + q * (g2 - q2)
  lam <- list(q, q1, q2, q3)

  hu <- list(u1, u2, u3, u4)
  bell <- function(mm, j) {
    switch(mm,
      hu[[1]],
      switch(j, hu[[2]], hu[[1]]^2),
      switch(j, hu[[3]], 3 * hu[[1]] * hu[[2]], hu[[1]]^3),
      switch(j, hu[[4]], 4 * hu[[1]] * hu[[3]] + 3 * hu[[2]]^2,
        6 * hu[[1]]^2 * hu[[2]], hu[[1]]^4)
    )
  }
  upow <- function(cc, j) {
    if (j == 0L) {
      return(u^cc)
    }
    acc <- 0
    for (r in seq_len(min(j, cc))) {
      acc <- acc + prod(cc - seq_len(r) + 1) * u^(cc - r) * bell(j, r)
    }
    acc
  }
  lamchain <- function(cc, p) {
    if (p == 0L) {
      return(lam[[cc]])
    }
    acc <- 0
    for (r in seq_len(p)) acc <- acc + lam[[cc + r]] * alpha^r * bell(p, r)
    acc
  }

  out <- vector("list", 15L)
  nms <- character(15L)
  k <- 0L
  for (cc in 0:4) {
    for (i in 0:(4L - cc)) {
      k <- k + 1L
      nms[k] <- paste0(cc, "_", i)
      out[[k]] <- if (cc == 0L) {
        if (i == 0L) NA_real_ else gz[[i]] + lamchain(0L, i)
      } else {
        acc <- 0
        for (p in 0:i) acc <- acc + choose(i, p) * upow(cc, i - p) * lamchain(cc, p)
        acc
      }
    }
  }
  names(out) <- nms
  c(list(z = z, a = a, rs = rs), out)
}

#' @title One Skew t Derivative in the Location, Scale and Shape
#' @name skewt_msa_component
#'
#' @description
#' Reads \eqn{\partial^a_\mu \partial^b_\sigma \partial^c_\alpha \ell} off the
#' table [skewt_msa_tower()] builds, by the \eqn{P_{a,b}} sum that function's
#' page derives.
#'
#' @details
#' A term whose coefficient \eqn{\prod_{k=i}^{b-1}(a+k)} vanishes is dropped
#' rather than multiplied, because the factor it would multiply is the entry
#' `"0_0"`, which the table does not hold: in R `0 * NA` is `NA`.
#'
#' @param tw A list from [skewt_msa_tower()].
#' @param sigma The scale, of length 1 or of the length of the response.
#' @param a,b,c Non-negative integers, the number of times the derivative is
#'   taken in \eqn{\mu}, in \eqn{\sigma} and in \eqn{\alpha}. Their sum is the
#'   order.
#'
#' @return A numeric vector, one value per observation.
#'
#' @seealso [skewt_msa_tower()] for the table and the derivation, and
#'   [skewt_msa_derivs()], which loops this over a whole order.
#'
#' @examples
#' tw <- distributions7:::skewt_msa_tower(c(-0.4, 1.2), 0, 1, 0.7, 8)
#' distributions7:::skewt_msa_component(tw, 1, a = 2L, b = 1L, c = 0L)
#'
#' @keywords internal
skewt_msa_component <- function(tw, sigma, a, b, c) {
  # the table holds the z-derivatives times rs^j, so z^i times the entry of
  # order a + i is (z/rs)^i times it, over rs^a
  za <- tw[["a"]]
  acc <- 0
  for (i in 0:b) {
    ris <- if (i > b - 1L) 1 else prod(a + (i:(b - 1L)))
    if (ris == 0) next
    acc <- acc + choose(b, i) * ris * za^i * tw[[paste0(c, "_", a + i)]]
  }
  val <- (-1)^(a + b) * acc / (tw[["rs"]]^a * sigma^(a + b))
  # The explicit -log(sigma) survives only where nothing differentiates z.
  if (a == 0L && c == 0L && b >= 1L) {
    val <- val + (-1)^b * factorial(b - 1L) / sigma^b
  }
  val
}

#' @title Every Skew t Derivative of an Order That Avoids the Degrees of Freedom
#' @name skewt_msa_derivs
#'
#' @description
#' Returns, in closed form, the components of a given order whose indices are
#' all drawn from \eqn{(\mu, \sigma, \alpha)}: ten of the twenty at order
#' three and fifteen of the thirty-five at order four.
#'
#' @details
#' The multi-indices and the names come from [deriv_indices()] and
#' [deriv_names()], the same enumeration, so a name is never recovered by
#' splitting a string.
#'
#' @param distrib A `SkewTDistrib` object, from [skewt_distrib()].
#' @param y A numeric vector of observations.
#' @param theta A named list with components `mu`, `sigma`, `alpha` and `nu`.
#' @param order A single integer, the order of differentiation.
#'
#' @return A named list of numeric vectors, a subset of [deriv_names()] at that
#'   order, holding every component free of \eqn{\nu}.
#'
#' @seealso [skewt_msa_tower()] for the derivation and
#'   [distrib_deriv3.SkewTDistrib()] for the compiled kernel tested against
#'   this.
#'
#' @examples
#' d <- skewt_distrib()
#' m <- distributions7:::skewt_msa_derivs(
#'   d, c(-0.4, 1.2), list(mu = 0, sigma = 1, alpha = 0.7, nu = 8), 3L)
#' names(m)
#'
#' @keywords internal
skewt_msa_derivs <- function(distrib, y, theta, order) {
  params <- distrib@params
  idx <- deriv_indices(params, order)
  nms <- deriv_names(params, order)
  keep <- vapply(idx, function(i) !any(i == 4L), logical(1))
  tw <- skewt_msa_tower(y, theta[[1]], theta[[2]], theta[[3]], theta[[4]])
  out <- lapply(idx[keep], function(i) {
    skewt_msa_component(tw, theta[[2]], sum(i == 1L), sum(i == 2L), sum(i == 3L))
  })
  names(out) <- nms[keep]
  out
}

#' @title Skew t Expected Hessian and Its Derivatives
#' @name distrib_expected_hessian.SkewTDistrib
#' @aliases distrib_dexpected_hessian.SkewTDistrib
#'   distrib_d2expected_hessian.SkewTDistrib
#' @description
#' Returns the expected information, and through [distrib_dexpected_hessian()]
#' and [distrib_d2expected_hessian()] its first and second derivatives in the
#' parameters.
#'
#' @details
#' No component has an elementary form. The family is a location-scale family,
#' so every component equals its value at \eqn{\mu = 0}, \eqn{\sigma = 1} times
#' \eqn{\sigma^{-k}}, with \eqn{k} the number of indices on \eqn{\mu} or
#' \eqn{\sigma}. The value there depends on \eqn{(\alpha, \nu)} alone and is one
#' integral over \eqn{z} of the observed derivatives against the density, taken
#' once per distinct pair by the exp-sinh rule of [loc_scale_expected()].
#'
#' The integral is exact to the rule's accuracy, and what it integrates are the
#' observed derivatives of [distrib_gradient.SkewTDistrib()], exact in every
#' parameter.
#'
#' The tail is integrated to \eqn{|z| = 10^{60}}, which leaves a relative
#' \eqn{10^{-60\nu}/\nu} and is negligible for \eqn{\nu \ge 0.3}.
#'
#' `approx` and `nsim` are accepted for the generic's sake and ignored.
#'
#' @param distrib A `SkewTDistrib` object, from [skewt_distrib()].
#' @param y A numeric vector of observations. Its length sets the length of
#'   each returned component; the values themselves are not read.
#' @param theta A named list with components `mu`, `sigma`, `alpha` and `nu`,
#'   each a numeric vector of length 1 or of the length of `y`.
#' @param scale One of `"parameter"` (the default) or `"link"`, matched by
#'   [base::match.arg()].
#' @param approx,nsim Ignored.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#' @param threads The thread count passed to the family's kernels.
#'
#' @return A named list of numeric vectors of length `length(y)`: for
#'   [distrib_expected_hessian()] the components keyed as [hess_names()], for
#'   the two derivatives those keyed as [dexpected_names()] and
#'   [d2expected_names()].
#'
#' @seealso [loc_scale_expected()] for the construction,
#'   [distrib_hessian.SkewTDistrib()] for the quantity this is the expectation
#'   of, and [distrib_expected_hessian()] for the generic.
#'
#' @examples
#' d <- skewt_distrib()
#' th <- list(mu = 1, sigma = 2, alpha = 2, nu = 5)
#' e <- distrib_expected_hessian(d, 0, th)
#' vapply(e, function(v) v[1], numeric(1))
#'
#' # One quadrature serves every observation sharing a shape: the location and
#' # the scale only rescale the result.
#' e2 <- distrib_expected_hessian(d, 0, list(mu = -3, sigma = 1, alpha = 2, nu = 5))
#' e$alpha_alpha / e2$alpha_alpha
NULL

register_loc_scale_expected(SkewTDistrib)

#' @title Skew t Third Derivatives
#' @name distrib_deriv3.SkewTDistrib
#'
#' @description
#' Computes the twenty third derivatives of the log-density with the compiled
#' kernel described on [distrib_gradient.SkewTDistrib()].
#'
#' @details
#' With `expected = TRUE` the whole order is an expectation and comes from
#' `expected_derivative()`; the family has no closed-form expected
#' information, so `approx` and `nsim` are read.
#'
#' @param distrib A `SkewTDistrib` object, from [skewt_distrib()].
#' @param y A numeric vector of observations. With `expected = TRUE` only its
#'   length matters.
#' @param theta A named list with components `mu`, `sigma`, `alpha` and `nu`.
#' @param expected Logical of length 1. When `FALSE`, the default, the observed
#'   derivatives at `y` are returned.
#' @param scale Either `"parameter"`, the default, or `"link"`. The
#'   transformation is applied in the generic's body.
#' @param approx One of `"integrate"`, `"bartlett"`, `"mc"` or `"opg"`, read
#'   only when `expected = TRUE`.
#' @param nsim A single positive integer, the Monte Carlo sample size used when
#'   `approx = "mc"`. Defaults to `10000`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#' @param threads A single positive integer, how many threads the kernel may
#'   use.
#'
#' @return A named list of twenty numeric vectors, one per distinct third-order
#'   component, from `mu_mu_mu` to `nu_nu_nu` as [deriv_names()] names them.
#'
#' @seealso [distrib_hessian.SkewTDistrib()] for the order below,
#'   [distrib_deriv4.SkewTDistrib()] for the order above, and
#'   [distrib_deriv3()] for the generic.
#'
#' @examples
#' d <- skewt_distrib()
#' y <- c(-1.5, -0.3, 0.4, 2.1)
#' th <- list(mu = 0, sigma = 1, alpha = 3, nu = 6)
#' d3 <- distrib_deriv3(d, y, th)
#' c(components = length(d3), involving_nu = sum(grepl("nu", names(d3))))
#'
#' # A closed-form-block component against a difference of the Hessian.
#' eps <- 1e-5
#' rbind(analytic = d3$mu_mu_alpha,
#'       numeric = (distrib_hessian(d, y, list(mu = 0, sigma = 1,
#'                                             alpha = 3 + eps, nu = 6))$mu_mu -
#'                  distrib_hessian(d, y, list(mu = 0, sigma = 1,
#'                                             alpha = 3 - eps, nu = 6))$mu_mu) /
#'                 (2 * eps))
#'
#' # The pure-nu component against a single stencil on the log-density.
#' ld <- function(v) sum(distrib_pdf(d, y, list(mu = 0, sigma = 1,
#'                                              alpha = 3, nu = v), log = TRUE))
#' c(ours = sum(d3$nu_nu_nu),
#'   stencil = numericals7::fd_derivative(ld, 6, 3L, h = 0.05))
S7::method(distrib_deriv3, SkewTDistrib) <- function(distrib, y, theta, expected = FALSE, scale = c("parameter", "link"), approx = c("integrate", "bartlett", "mc", "opg"), nsim = 10000, ...,
                                                  threads = 1L) {
  if (expected) {
    return(expected_derivative(distrib, y, theta, order = 3L,
                               approx = match.arg(approx), nsim = nsim))
  }
  skewt_deriv3_cpp(y, theta[[1]], theta[[2]], theta[[3]], theta[[4]], threads)
}

#' @title Skew t Fourth Derivatives
#' @name distrib_deriv4.SkewTDistrib
#'
#' @description
#' Computes the thirty-five fourth derivatives of the log-density with the
#' compiled kernel described on [distrib_gradient.SkewTDistrib()].
#'
#' @details
#' With `expected = TRUE` the whole order is an expectation and comes from
#' `expected_derivative()`; `approx` and `nsim` are then read.
#'
#' @param distrib A `SkewTDistrib` object, from [skewt_distrib()].
#' @param y A numeric vector of observations. With `expected = TRUE` only its
#'   length matters.
#' @param theta A named list with components `mu`, `sigma`, `alpha` and `nu`.
#' @param expected Logical of length 1. When `FALSE`, the default, the observed
#'   derivatives at `y` are returned.
#' @param scale Either `"parameter"`, the default, or `"link"`. The
#'   transformation is applied in the generic's body.
#' @param approx One of `"integrate"`, `"bartlett"`, `"mc"` or `"opg"`, read
#'   only when `expected = TRUE`.
#' @param nsim A single positive integer, the Monte Carlo sample size used when
#'   `approx = "mc"`. Defaults to `10000`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#' @param threads A single positive integer, how many threads the kernel may
#'   use.
#'
#' @return A named list of thirty-five numeric vectors, one per distinct
#'   fourth-order component, from `mu_mu_mu_mu` to `nu_nu_nu_nu`.
#'
#' @seealso [distrib_deriv3.SkewTDistrib()] for the order below,
#'   [distrib_deriv5.SkewTDistrib()] for the order above, and
#'   [distrib_deriv4()] for the generic.
#'
#' @examples
#' d <- skewt_distrib()
#' y <- c(-1.5, -0.3, 0.4, 2.1)
#' th <- list(mu = 0, sigma = 1, alpha = 3, nu = 6)
#' d4 <- distrib_deriv4(d, y, th)
#' c(components = length(d4), involving_nu = sum(grepl("nu", names(d4))))
#'
#' # The components span five orders of magnitude, so a relative comparison
#' # on the smallest of them says nothing about the largest.
#' s <- sort(vapply(d4, function(v) sum(abs(v)), 0), decreasing = TRUE)
#' s[c(1, 2, length(s) - 1, length(s))]
#'
#' # A closed-form-block component against a difference of the third order.
#' eps <- 1e-5
#' rbind(analytic = d4$mu_mu_alpha_alpha,
#'       numeric = (distrib_deriv3(d, y, list(mu = 0, sigma = 1,
#'                                            alpha = 3 + eps, nu = 6))$mu_mu_alpha -
#'                  distrib_deriv3(d, y, list(mu = 0, sigma = 1,
#'                                            alpha = 3 - eps, nu = 6))$mu_mu_alpha) /
#'                 (2 * eps))
S7::method(distrib_deriv4, SkewTDistrib) <- function(distrib, y, theta, expected = FALSE, scale = c("parameter", "link"), approx = c("integrate", "bartlett", "mc", "opg"), nsim = 10000, ...,
                                                  threads = 1L) {
  if (expected) {
    return(expected_derivative(distrib, y, theta, order = 4L,
                               approx = match.arg(approx), nsim = nsim))
  }
  skewt_deriv4_cpp(y, theta[[1]], theta[[2]], theta[[3]], theta[[4]], threads)
}

#' @title Skew t Fifth Derivatives
#' @name distrib_deriv5.SkewTDistrib
#'
#' @description
#' Computes the fifty-six fifth derivatives of the log-density with the
#' compiled kernel described on [distrib_gradient.SkewTDistrib()].
#'
#' @param distrib A `SkewTDistrib` object, from [skewt_distrib()].
#' @param y A numeric vector of observations.
#' @param theta A named list with components `mu`, `sigma`, `alpha` and `nu`.
#' @param scale Either `"parameter"`, the default, or `"link"`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#' @param threads A single positive integer, how many threads the kernel may
#'   use.
#'
#' @return A named list of fifty-six numeric vectors, one per distinct
#'   fifth-order component, as [deriv_names()] names them.
#'
#' @seealso [distrib_deriv4.SkewTDistrib()] for the order below and
#'   [distrib_deriv5()] for the generic.
#'
#' @examples
#' d <- skewt_distrib()
#' d5 <- distrib_deriv5(d, c(-1.5, 0.4, 2.1),
#'                      list(mu = 0, sigma = 1, alpha = 3, nu = 6))
#' length(d5)
S7::method(distrib_deriv5, SkewTDistrib) <- function(
    distrib, y, theta, scale = c("parameter", "link"), ..., threads = 1L) {
  deriv5_scale(distrib, y, theta,
               skewt_deriv5_cpp(y, theta[[1]], theta[[2]], theta[[3]], theta[[4]], threads),
               match.arg(scale))
}

#' @title Skew t Response Derivative
#' @name distrib_grad_y.SkewTDistrib
#'
#' @description
#' Computes \eqn{\partial\ell/\partial y = D/\sigma}, with \eqn{D = A + QB} in
#' the notation of [skewt_pieces()]. It is minus the derivative in \eqn{\mu},
#' exactly, because the response and the location enter the density only
#' through their difference; the identity holds for every location family.
#'
#' It is closed form at every parameter value, \eqn{\nu} included: nothing here
#' differentiates the degrees of freedom.
#'
#' @param distrib A `SkewTDistrib` object, from [skewt_distrib()].
#' @param y A numeric vector of observations.
#' @param theta A named list with components `mu`, `sigma`, `alpha` and `nu`,
#'   each a numeric vector of length 1 or of the length of `y`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of the length of the recycled inputs.
#'
#' @section Notation:
#' \eqn{\ell} is the log-density of one observation and \eqn{A}, \eqn{B},
#' \eqn{Q} are as [skewt_pieces()] defines them.
#'
#' @seealso [distrib_hess_y.SkewTDistrib()] for the second derivative,
#'   [distrib_gradient.SkewTDistrib()] for the parameter derivatives, and
#'   [distrib_grad_y()] for the generic.
#'
#' @examples
#' d <- skewt_distrib()
#' y <- c(-1.5, -0.3, 0.4, 2.1)
#' th <- list(mu = 0, sigma = 1, alpha = 3, nu = 6)
#'
#' # It is minus the score in the location, exactly.
#' all.equal(distrib_grad_y(d, y, th), -distrib_gradient(d, y, th)$mu)
#'
#' # Against a central difference of the log-density.
#' eps <- 1e-6
#' rbind(analytic = distrib_grad_y(d, y, th),
#'       numeric = (distrib_pdf(d, y + eps, th, log = TRUE) -
#'                  distrib_pdf(d, y - eps, th, log = TRUE)) / (2 * eps))
#'
#' # The score redescends in the heavy tail, as a Student t's does.
#' distrib_grad_y(d, c(2, 8, 32, 128), th)
S7::method(distrib_grad_y, SkewTDistrib) <- function(distrib, y, theta, ...) {
  p <- skewt_pieces(y, theta[[1]], theta[[2]], theta[[3]], theta[[4]])
  (p$a + p$q * p$b) / theta[[2]]
}

#' @title Skew t Second Response Derivative
#' @name distrib_hess_y.SkewTDistrib
#'
#' @description
#' Computes \eqn{\partial^2\ell/\partial y^2 = D'/\sigma^2}, with
#' \eqn{D' = A' + Q'B^2 + QB'} in the notation of [skewt_pieces()]. It is the
#' same expression as \eqn{\partial^2\ell/\partial\mu^2}: two signs cancel
#' where one did not at first order.
#'
#' Unlike the skew normal's, the value **changes sign**. The Student
#' \eqn{t}'s curvature is positive wherever \eqn{|z| > \sqrt{\nu}}, so the
#' log-density is convex in the far tail. That convexity makes the family's
#' score redescend, and its estimates resistant to an outlier.
#'
#' @param distrib A `SkewTDistrib` object, from [skewt_distrib()].
#' @param y A numeric vector of observations.
#' @param theta A named list with components `mu`, `sigma`, `alpha` and `nu`,
#'   each a numeric vector of length 1 or of the length of `y`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of the length of the recycled inputs. It is not of
#'   one sign.
#'
#' @section Notation:
#' \eqn{z = (y-\mu)/\sigma}, \eqn{\nu} the degrees of freedom, and \eqn{A'},
#' \eqn{B}, \eqn{B'}, \eqn{Q}, \eqn{Q'} are as [skewt_pieces()] defines them.
#'
#' @seealso [distrib_grad_y.SkewTDistrib()] for the first derivative,
#'   [distrib_hessian.SkewTDistrib()] for the parameter curvature it shares an
#'   expression with, [distrib_hess_y.SkewNormal1Distrib()] for the case that
#'   does not change sign, and [distrib_hess_y()] for the generic.
#'
#' @examples
#' d <- skewt_distrib()
#' y <- c(-1.5, -0.3, 0.4, 2.1)
#' th <- list(mu = 0, sigma = 1, alpha = 3, nu = 6)
#'
#' # It equals the curvature in the location, without a sign change.
#' all.equal(distrib_hess_y(d, y, th), distrib_hessian(d, y, th)$mu_mu)
#'
#' # Against a central difference of the response derivative.
#' eps <- 1e-5
#' rbind(analytic = distrib_hess_y(d, y, th),
#'       numeric = (distrib_grad_y(d, y + eps, th) -
#'                  distrib_grad_y(d, y - eps, th)) / (2 * eps))
#'
#' # Convex in the far tail, where the skew normal is concave everywhere.
#' sn <- skewnormal1_distrib()
#' far <- c(-20, -8, 8, 20)
#' rbind(skew_t = distrib_hess_y(d, far, th),
#'       skew_normal = distrib_hess_y(sn, far, list(mu = 0, sigma = 1, alpha = 3)))
S7::method(distrib_hess_y, SkewTDistrib) <- function(distrib, y, theta, ...) {
  p <- skewt_pieces(y, theta[[1]], theta[[2]], theta[[3]], theta[[4]])
  (p$da + p$dq * p$b^2 + p$q * p$db) / theta[[2]]^2
}

# --- CONSTRUCTOR WRAPPER ---

#' @title Skew t Distribution Object
#'
#' @description
#' Builds a skew \eqn{t} distribution object with location \eqn{\mu}, scale
#' \eqn{\sigma}, shape \eqn{\alpha} and degrees of freedom \eqn{\nu}. It is
#' the four-parameter family a location-scale-shape framework wants: the
#' scale, the asymmetry and the tail weight are three separate parameters,
#' each of which can carry its own linear predictor.
#'
#' Three families sit inside it. \eqn{\alpha = 0} is
#' [student_t1_distrib()], large \eqn{\nu} approaches
#' [skewnormal1_distrib()], and both together give the Gaussian.
#'
#' @param link_mu A `linkfunctions7` link object for the location \eqn{\mu},
#'   which is unconstrained. Defaults to [linkfunctions7::identity_link()].
#' @param link_sigma A link object for the scale \eqn{\sigma}, which must be
#'   strictly positive. Defaults to [linkfunctions7::log_link()].
#' @param link_alpha A link object for the shape \eqn{\alpha}, which is
#'   unconstrained. Defaults to [linkfunctions7::identity_link()].
#' @param link_nu A link object for the degrees of freedom \eqn{\nu}, which
#'   must be strictly positive. Defaults to [linkfunctions7::log_link()].
#'
#' @details
#' # Why the fourth parameter
#'
#' The skew normal of [skewnormal1_distrib()] reaches a skewness of at most
#' 0.99527 and an excess kurtosis of at most 0.86918, both approached only as
#' \eqn{|\alpha| \to \infty}. Adding \eqn{\nu} removes both ceilings: measured
#' at \eqn{\alpha = 50}, the skewness is 1.190 at \eqn{\nu = 30}, 2.050 at
#' \eqn{\nu = 6} and 3.998 at \eqn{\nu = 4}.
#'
#' # How the derivatives are computed
#'
#' The density contains \eqn{T_{\nu+1}}, and the derivatives of a Student
#' \eqn{t} distribution function in its degrees of freedom have no elementary
#' expression. They are integrals of the derivatives of the \eqn{t} density and
#' are taken by quadrature, so that every derivative this family reports, to
#' order five, is exact to rounding; see [distrib_gradient.SkewTDistrib()].
#'
#' # Fitting
#'
#' The expected information has no closed form and is computed by one
#' quadrature over \eqn{z} per distinct \eqn{(\alpha, \nu)}, so Fisher scoring
#' on an intercept-only fit pays for one integral per iteration.
#'
#' The distribution function and the quantile function likewise have no
#' elementary form; the base class integrates the density and inverts the
#' result by root finding.
#'
#' # Moments
#'
#' They exist only up to order \eqn{\nu}: the mean requires \eqn{\nu > 1}, the
#' variance \eqn{\nu > 2}, the skewness \eqn{\nu > 3} and the excess kurtosis
#' \eqn{\nu > 4}, and each returns `NaN` below its threshold. The density is
#' well defined at every positive \eqn{\nu}, which is why the moments and the
#' parameters are kept apart.
#'
#' # Identification near symmetry
#'
#' The skew normal's expected information loses a rank at \eqn{\alpha = 0},
#' because its shape score is a fixed multiple of its location score there.
#' **The skew \eqn{t} does not inherit that at finite \eqn{\nu.}** The tilting
#' argument carries \eqn{c = \sqrt{(\nu+1)/(\nu+z^2)}}, which depends on the
#' observation, so the two scores are not proportional: measured at
#' \eqn{\alpha = 0}, \eqn{\nu = 6}, their ratio spreads over 0.76 across
#' observations, where the skew normal's is constant to the last digit. The
#' spread falls as \eqn{O(1/\nu)}, so the singularity is inherited only in the
#' limit: 0.068 at \eqn{\nu = 100}, \eqn{7\times10^{-4}} at \eqn{10^4}.
#'
#' What is weakly identified here is \eqn{\nu} itself. The smallest eigenvalue
#' of the information belongs to that direction and falls with \eqn{\nu}:
#' measured at \eqn{\alpha = 0}, it is \eqn{3.6\times10^{-3}} at \eqn{\nu = 4},
#' \eqn{8.4\times10^{-4}} at 6 and \eqn{1.7\times10^{-6}} at 30. A nearly
#' Gaussian tail carries little information about how heavy it is.
#'
#' # Parameter domains
#'
#' - \eqn{\mu \in (-\infty, \infty)}
#' - \eqn{\sigma \in (0, \infty)}
#' - \eqn{\alpha \in (-\infty, \infty)}
#' - \eqn{\nu \in (0, \infty)}
#'
#' @return An S7 object of class [SkewTDistrib], inheriting from
#'   `continuous_distrib`. Its `params` are
#'   `c("mu", "sigma", "alpha", "nu")`, its `bounds` `c(-Inf, Inf)`, and its
#'   `link_params` the four links given here.
#'
#' @references
#' Azzalini, A. and Capitanio, A. (2003). Distributions generated by
#' perturbation of symmetry with emphasis on a multivariate skew t distribution.
#' *Journal of the Royal Statistical Society, Series B* 65, 367-389.
#'
#' @importFrom linkfunctions7 identity_link log_link
#' @importFrom stats dt pt rchisq
#'
#' @seealso [skewnormal1_distrib()] for the \eqn{\nu \to \infty} limit,
#'   [student_t1_distrib()] for the \eqn{\alpha = 0} case,
#'   [skewt_pieces()] for the scalar functions the derivatives are built from,
#'   and [SkewTDistrib] for the class and its method list.
#'
#' @examples
#' d <- skewt_distrib()
#' d@params
#' th <- list(mu = 0, sigma = 1, alpha = 3, nu = 5)
#'
#' distrib_pdf(d, c(-1, 0, 1), th)
#'
#' # Shape zero is the Student t.
#' all.equal(distrib_pdf(d, c(-1, 0, 1),
#'                       list(mu = 0, sigma = 1, alpha = 0, nu = 5)),
#'           dt(c(-1, 0, 1), df = 5))
#'
#' # It passes both of the skew normal's ceilings.
#' rbind(skew_t = c(skewness(d, list(mu = 0, sigma = 1, alpha = 8, nu = 5)),
#'                  kurtosis(d, list(mu = 0, sigma = 1, alpha = 8, nu = 5))),
#'       skew_normal_bound = c(0.99527, 0.86918))
#'
#' # Moments exist only up to order nu.
#' t(vapply(c(1.5, 2.5, 3.5, 4.5), function(v) {
#'   p <- list(mu = 0, sigma = 1, alpha = 3, nu = v)
#'   c(nu = v, mean = mean(d, p), var = variance(d, p),
#'     skew = skewness(d, p), kurt = kurtosis(d, p))
#' }, numeric(5)))
#'
#' # Fitting by Fisher scoring, the default, and by Newton on the observed
#' # Hessian reaches the same estimate.
#' set.seed(1)
#' x <- distrib_rng(d, 200, th)
#' rbind(fisher = coef(fit_distrib(d, x, start = th)),
#'       newton = coef(fit_distrib(d, x, method = optimizers7::newton(), start = th)))
#'
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

#' @title Skew t Distribution Function
#' @name distrib_cdf.SkewTDistrib
#'
#' @description
#' Computes \eqn{P(Y \le q)} by numerical integration of the density, taken
#' by the compiled rule of [compiled_cdf()] over the tail on the side of
#' \eqn{q} away from \eqn{\mu}.
#'
#' @param distrib A `SkewTDistrib` object, from [skewt_distrib()].
#' @param q A numeric vector of quantiles.
#' @param theta A named list with components `mu`, `sigma`, `alpha` and `nu`.
#' @param lower.tail Logical; if `FALSE`, \eqn{P(Y > q)} is returned.
#' @param log.p Logical; if `TRUE`, the logarithm is returned.
#' @param ... Unused.
#'
#' @return A numeric vector of probabilities.
#'
#' @examples
#' d <- skewt_distrib()
#' th <- list(mu = 0.3, sigma = 1.2, alpha = 2, nu = 6)
#' distrib_cdf(d, c(-1, 0.5, 2), th)
#' c(method = distrib_cdf(d, 2, th),
#'   integral = integrate(function(v) distrib_pdf(d, v, th), -Inf, 2)$value)
#'
#' @keywords internal
S7::method(distrib_cdf, SkewTDistrib) <- compiled_cdf
