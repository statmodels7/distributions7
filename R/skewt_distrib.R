#' @include distrib.R generics.R skewnormal1_distrib.R
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
#'   [`distrib_grad_y()`][distrib_grad_y.SkewTDistrib],
#'   [`distrib_hess_y()`][distrib_hess_y.SkewTDistrib].
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
#' [continuous_distrib()], by quadrature and by root finding on it. So does the
#' **expected information**: this family has none in elementary form, so
#' [distrib_expected_hessian()] approximates it, and
#' `method = "newton"` is much the cheaper way to fit it.
#'
#' @section What is closed form and what is not:
#' Every derivative in \eqn{\mu}, \eqn{\sigma} and \eqn{\alpha} is closed form.
#' Every derivative involving \eqn{\nu} is not, and the obstruction is
#' mathematical: the density carries \eqn{T_{\nu+1}}, and the derivative of a
#' Student \eqn{t} distribution function in its degrees of freedom has no
#' elementary expression. Those components come from **one** stencil applied to
#' an analytic quantity, never from a difference of a difference.
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
#' Nothing here involves \eqn{\nu} by differentiation; the components in
#' \eqn{\nu} are obtained separately, by a stencil.
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

#' @title The Step a Skew t Differences the Degrees of Freedom With
#'
#' @description
#' Returns the finite-difference step used for the derivatives in \eqn{\nu},
#' `pmax(1e-3 * abs(nu), 1e-6)`: relative to \eqn{\nu} itself, so that the same
#' number of significant digits is differenced at every scale, and floored so
#' that it stays a number for a degrees of freedom near zero.
#'
#' @details
#' The relative step \eqn{10^{-3}} is measured rather than assumed. Swept over
#' \eqn{\nu} from 2 to 30 and sample sizes from 500 to 4000, it is where the
#' truncation error of [fd5_first()]'s five-point stencil has fallen to the
#' level of the rounding error and the two are balanced. A smaller step is
#' dominated by rounding, which the stencil amplifies by \eqn{18/(12h)}, and a
#' larger one by truncation, which grows as \eqn{h^4}.
#'
#' The choice shows up in the stopping rule of a fit. [fit_distrib()] takes
#' its rule from the method it is given, and `optimizers7::crit_grad()` tests
#' the score **per observation** at a default tolerance of \eqn{10^{-6}}.
#' Measured on samples of 500 to 4000 the run converges at between
#' \eqn{3\times10^{-10}} and \eqn{9\times10^{-9}} per observation, so this
#' step leaves two orders of room. A three-point stencil would not: its bias
#' does not cancel over the sum, and a run would spend its whole budget
#' reporting failure at the maximum.
#'
#' @param nu A numeric vector of degrees of freedom.
#'
#' @return A numeric vector of steps, of the length of `nu`.
#'
#' @seealso [fd5_first()] and its siblings, which consume the step, and
#'   [distrib_gradient.SkewTDistrib()], the first method that needs it.
#'
#' @examples
#' distributions7:::skewt_nu_step(c(2, 6, 30, 1e-5))
#'
#' # Relative above the floor, absolute below it.
#' nu <- c(1e-4, 1e-3, 1e-2, 1, 100)
#' rbind(nu = nu, step = distributions7:::skewt_nu_step(nu))
#'
#' @keywords internal
skewt_nu_step <- function(nu) {
  pmax(1e-3 * abs(nu), 1e-6)
}

#' @title A Five-Point First Derivative
#'
#' @description
#' Returns \eqn{\{f(x-2h) - 8f(x-h) + 8f(x+h) - f(x+2h)\}/(12h)}, the central
#' stencil on five nodes that is exact on polynomials up to degree four, so its
#' truncation error is \eqn{O(h^4)}. The weights come from
#' [numericals7::fd_derivative()] at `accuracy = 4`, which builds them from the
#' Vandermonde system on the offsets \eqn{-2, -1, 0, 1, 2}.
#'
#' @details
#' The three-point stencil is not accurate enough for the derivative in
#' \eqn{\nu} of a fitted likelihood. Its truncation error is \eqn{O(h^2)} per
#' observation and does **not** cancel when the observations are summed,
#' because it is a bias: on a sample of a few thousand it leaves
#' the summed score at about \eqn{10^{-8}}, an order of magnitude worse than
#' the five-point stencil, which reaches the level of rounding at the cost of
#' two more evaluations of `f`.
#'
#' This is one stencil applied to an analytic quantity. Nothing in this family
#' differences a differenced value.
#'
#' @param f A function of one scalar, returning a numeric vector. It is called
#'   four times, at \eqn{x \pm h} and \eqn{x \pm 2h}.
#' @param x A single number, the point to differentiate at.
#' @param h A single positive number, the step. For this family it comes from
#'   [skewt_nu_step()].
#'
#' @return A numeric vector, of whatever length `f` returns.
#'
#' @seealso [fd5_second()], [fd5_third()] and [fd5_fourth()] for the other
#'   orders, [skewt_nu_step()] for the step, and
#'   [numericals7::fd_derivative()] for the stencil library.
#'
#' @examples
#' f <- function(x) exp(x) * sin(x)
#' truth <- exp(0.7) * (sin(0.7) + cos(0.7))
#' c(stencil = distributions7:::fd5_first(f, 0.7, 1e-3), truth = truth)
#'
#' # The error falls as h^4 until rounding takes over.
#' vapply(c(0.1, 0.05, 0.025),
#'        function(h) abs(distributions7:::fd5_first(f, 0.7, h) - truth), 0)
#'
#' @keywords internal
fd5_first <- function(f, x, h) {
  # numericals7's shared weights at accuracy four: the displayed formula is
  # exactly what the Vandermonde construction produces on five nodes.
  numericals7::fd_derivative(f, x, 1L, h = h, accuracy = 4L)
}

#' @title A Five-Point Second Derivative
#'
#' @description
#' Returns \eqn{\{-f(x-2h) + 16f(x-h) - 30f(x) + 16f(x+h) - f(x+2h)\}/(12h^2)},
#' the central stencil on five nodes for the second derivative, with truncation
#' error \eqn{O(h^4)}. Like [fd5_first()] it comes from
#' [numericals7::fd_derivative()] at `accuracy = 4`.
#'
#' Rounding is amplified by \eqn{h^{-2}} here, one power more than in the first
#' derivative, so the attainable accuracy at the same step is one or two digits
#' lower.
#'
#' @param f A function of one scalar, returning a numeric vector. It is called
#'   five times, at \eqn{x}, \eqn{x \pm h} and \eqn{x \pm 2h}.
#' @param x A single number, the point to differentiate at.
#' @param h A single positive number, the step.
#'
#' @return A numeric vector, of whatever length `f` returns.
#'
#' @seealso [fd5_first()] for the order below, [fd5_third()] for the order
#'   above, and [distrib_hessian.SkewTDistrib()], which uses this for the
#'   \eqn{\nu} components.
#'
#' @examples
#' f <- function(x) exp(x) * sin(x)
#' c(stencil = distributions7:::fd5_second(f, 0.7, 1e-3),
#'   truth = 2 * exp(0.7) * cos(0.7))
#'
#' @keywords internal
fd5_second <- function(f, x, h) {
  numericals7::fd_derivative(f, x, 2L, h = h, accuracy = 4L)
}

#' @title A Five-Point Third Derivative
#'
#' @description
#' Returns \eqn{\{-f(x-2h)/2 + f(x-h) - f(x+h) + f(x+2h)/2\}/h^3}, the central
#' stencil on five nodes for the third derivative. Five nodes is the smallest
#' number that carries a third derivative at all, so the accuracy here is
#' \eqn{O(h^2)} where the first two orders get \eqn{O(h^4)} from the same
#' offsets.
#'
#' One stencil applied to an analytic quantity, never a difference of
#' differences.
#'
#' @param f A function of one scalar, returning a numeric vector. It is called
#'   four times, at \eqn{x \pm h} and \eqn{x \pm 2h}.
#' @param x A single number, the point to differentiate at.
#' @param h A single positive number, the step.
#'
#' @return A numeric vector, of whatever length `f` returns.
#'
#' @seealso [fd5_second()] for the order below, [fd5_fourth()] for the order
#'   above, and [distrib_deriv3.SkewTDistrib()], which uses this.
#'
#' @examples
#' f <- function(x) exp(x) * sin(x)
#' c(stencil = distributions7:::fd5_third(f, 0.7, 1e-2),
#'   truth = 2 * exp(0.7) * (cos(0.7) - sin(0.7)))
#'
#' @keywords internal
fd5_third <- function(f, x, h) {
  numericals7::fd_derivative(f, x, 3L, h = h)
}

#' @title A Five-Point Fourth Derivative
#'
#' @description
#' Returns \eqn{\{f(x-2h) - 4f(x-h) + 6f(x) - 4f(x+h) + f(x+2h)\}/h^4}, the
#' central stencil on five nodes for the fourth derivative, accurate to
#' \eqn{O(h^2)}.
#'
#' Rounding is amplified by \eqn{h^{-4}}, which sets what this can deliver:
#' measured on \eqn{e^x\sin x} at \eqn{x = 0.7} with \eqn{h = 10^{-2}} it
#' returns \eqn{-5.18939} against a true \eqn{-5.18918}, four significant
#' digits. On a component whose value is itself small the surviving digits are
#' fewer, and [distrib_deriv4.SkewTDistrib()] says so where a reader meets one.
#'
#' @param f A function of one scalar, returning a numeric vector. It is called
#'   five times, at \eqn{x}, \eqn{x \pm h} and \eqn{x \pm 2h}.
#' @param x A single number, the point to differentiate at.
#' @param h A single positive number, the step. A step too small is worse than
#'   one too large here, the rounding growing four times faster than the
#'   truncation falls.
#'
#' @return A numeric vector, of whatever length `f` returns.
#'
#' @seealso [fd5_third()] for the order below and
#'   [distrib_deriv4.SkewTDistrib()] for the method that uses this.
#'
#' @examples
#' f <- function(x) exp(x) * sin(x)
#' truth <- -4 * exp(0.7) * sin(0.7)
#' c(stencil = distributions7:::fd5_fourth(f, 0.7, 1e-2), truth = truth)
#'
#' # Too small a step is worse than too large: rounding grows as h^-4.
#' vapply(c(1e-1, 1e-2, 1e-3, 1e-4),
#'        function(h) abs(distributions7:::fd5_fourth(f, 0.7, h) - truth), 0)
#'
#' @keywords internal
fd5_fourth <- function(f, x, h) {
  numericals7::fd_derivative(f, x, 4L, h = h)
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
#' The fourth is not. \eqn{\partial\log T_{\nu+1}(w)/\partial\nu} is a
#' derivative of a Student \eqn{t} distribution function with respect to its
#' degrees of freedom, which has no elementary expression, the same obstruction
#' the gamma and beta distribution functions meet in their shape. That one
#' component is a single central difference of the **log-density**, taken with
#' [fd5_first()] at the step of [skewt_nu_step()].
#'
#' @details
#' # Accuracy
#'
#' Measured at \eqn{\mu = 0}, \eqn{\sigma = 1}, \eqn{\alpha = 3},
#' \eqn{\nu = 6} on four observations, the summed score agrees with
#' `numDeriv::grad` on the log-likelihood to \eqn{6\times10^{-12}} in
#' \eqn{\mu}, \eqn{2\times10^{-11}} in \eqn{\sigma},
#' \eqn{3\times10^{-12}} in \eqn{\alpha} and \eqn{5\times10^{-11}} in
#' \eqn{\nu}. The \eqn{\nu} component is therefore the loosest of the four,
#' and it is what sets the accuracy a fit of this family can stop at; see
#' [skewt_nu_step()].
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
#'   [distrib_hessian.SkewTDistrib()] for the second derivatives,
#'   [fd5_first()] for the stencil, and [distrib_gradient()] for the generic.
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

#' @title Skew t Observed Hessian
#' @name distrib_hessian.SkewTDistrib
#'
#' @description
#' Computes the ten second derivatives of the log-density. The block in
#' \eqn{(\mu, \sigma, \alpha)} is closed form; every component involving
#' \eqn{\nu} comes from one stencil applied to an analytic quantity, for the
#' reason [distrib_gradient.SkewTDistrib()] gives.
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
#' # The four components in the degrees of freedom
#'
#' \eqn{\partial^2\ell/\partial\nu^2} is one five-point second difference of
#' the log-density, [fd5_second()]. The three mixed ones step the
#' **closed-form score** in \eqn{\nu} with [fd5_first()], so only one
#' difference is taken and it is taken of an analytic quantity. Stepping the
#' log-density in both variables instead would be a difference of a difference
#' in \eqn{\nu}, which this family never does.
#'
#' Measured at \eqn{\mu = 0}, \eqn{\sigma = 1}, \eqn{\alpha = 3},
#' \eqn{\nu = 6}, the summed `nu_nu` agrees with `numDeriv::hessian` on the
#' log-likelihood to \eqn{2\times10^{-9}} relative and `mu_nu` to the printed
#' digit.
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
#'   [distrib_deriv3.SkewTDistrib()] for the order above,
#'   [fd5_second()] for the stencil, and [distrib_hessian()] for the generic.
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

#' @title Skew t Third Derivatives
#' @name distrib_deriv3.SkewTDistrib
#'
#' @description
#' Computes the twenty third derivatives of the log-density, assembled so that
#' no stencil is ever applied to another stencil's output. Ten of the twenty
#' involve \eqn{\nu}.
#'
#' @details
#' # How the twenty are obtained
#'
#' A component whose Hessian entry is closed form goes through the generic
#' construction of [numerical_deriv3()], which is one stencil on an analytic
#' quantity. That covers both indices in \eqn{(\mu, \sigma, \alpha)}, and also
#' one index equal to \eqn{\nu} where the stencil runs along a different
#' variable.
#'
#' The components the generic construction would nest are replaced:
#' \eqn{(i, \nu, \nu)} for \eqn{i} in \eqn{(\mu, \sigma, \alpha)} is one
#' five-point second difference of the **closed-form** score component
#' \eqn{i}, through [fd5_second()]; and \eqn{(\nu, \nu, \nu)} is one
#' five-point third difference of the log-density itself, through
#' [fd5_third()].
#'
#' # Accuracy
#'
#' The pure-\eqn{\nu} component is the loosest at this order. Measured at
#' \eqn{\mu = 0}, \eqn{\sigma = 1}, \eqn{\alpha = 3}, \eqn{\nu = 6} on four
#' observations, `nu_nu_nu` is \eqn{-0.0061280} against \eqn{-0.0061270} from
#' an independent single stencil on the log-density at \eqn{h = 0.05}, so
#' about four significant digits. The package's own thirteen-check battery
#' reports order 3 against finite differences at \eqn{2\times10^{-6}} for this
#' family.
#'
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
#'
#' @return A named list of twenty numeric vectors, one per distinct third-order
#'   component, from `mu_mu_mu` to `nu_nu_nu` as [deriv_names()] names them.
#'
#' @seealso [distrib_hessian.SkewTDistrib()] for the order below,
#'   [distrib_deriv4.SkewTDistrib()] for the order above,
#'   [fd5_third()] for the stencil, and [distrib_deriv3()] for the generic.
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
#' # The pure-nu component against an independent single stencil on the
#' # log-density, which shares no arithmetic with the route above.
#' ld <- function(v) sum(distrib_pdf(d, y, list(mu = 0, sigma = 1,
#'                                              alpha = 3, nu = v), log = TRUE))
#' c(ours = sum(d3$nu_nu_nu),
#'   stencil = numericals7::fd_derivative(ld, 6, 3L, h = 0.05))
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

#' @title Skew t Fourth Derivatives
#' @name distrib_deriv4.SkewTDistrib
#'
#' @description
#' Computes the thirty-five fourth derivatives of the log-density, with the
#' discipline of [distrib_deriv3.SkewTDistrib()]: the generic construction
#' serves every component whose Hessian entry is closed form, and the ones it
#' would nest are replaced by one stencil each. \eqn{(i, \nu, \nu, \nu)}
#' becomes a third difference of the closed-form score component \eqn{i}, and
#' \eqn{(\nu, \nu, \nu, \nu)} a fourth difference of the log-density. Twenty of
#' the thirty-five involve \eqn{\nu}.
#'
#' @details
#' # The step for the pure-nu component
#'
#' A fourth difference amplifies rounding by \eqn{h^{-4}}, so
#' [fd5_fourth()] is called at **ten times** [skewt_nu_step()]'s step rather
#' than at it. The choice is measured: at the family's base step the
#' per-observation noise is near \eqn{10^{-2}} relative, and at ten times that
#' it is negligible while the \eqn{O(h^2)} truncation, about
#' \eqn{6\times10^{-4}}, is what remains.
#'
#' `nu_nu_nu_nu` is the least accurate quantity this family reports, and it is
#' also the smallest: measured on four observations at \eqn{\nu = 6} it is
#' \eqn{-4.5\times10^{-6}} while `sigma_sigma_sigma_sigma` is 127. A
#' **relative** comparison on a component that small is not informative, which
#' is why the package's battery scales its order-4 check by the size of the
#' whole array and reports \eqn{3.5\times10^{-4}} here.
#'
#' # Cost
#'
#' This is the dearest method in the family: the twenty \eqn{\nu} components
#' each cost four or five evaluations of an analytic quantity over the whole
#' vector. Measured at \eqn{n = 20{,}000} it takes about sixteen seconds,
#' against sixty milliseconds for the score.
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
#'
#' @return A named list of thirty-five numeric vectors, one per distinct
#'   fourth-order component, from `mu_mu_mu_mu` to `nu_nu_nu_nu`.
#'
#' @seealso [distrib_deriv3.SkewTDistrib()] for the order below,
#'   [fd5_fourth()] for the stencil and its \eqn{h^{-4}} behavior, and
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
#' # What is closed form and what is not
#'
#' Every derivative in \eqn{(\mu, \sigma, \alpha)} is closed form. Everything
#' involving \eqn{\nu} is not, because the density contains \eqn{T_{\nu+1}} and
#' the derivative of a Student \eqn{t} distribution function with respect to
#' its degrees of freedom has no elementary expression. It is the same
#' obstruction that stops the gamma and beta distribution functions from having
#' closed-form shape derivatives.
#'
#' Those components come from a single stencil applied to an analytic quantity,
#' never from a difference of a difference:
#'
#' \tabular{lll}{
#'   **component** \tab **route** \tab **agreement with an independent route**
#'     \cr
#'   \eqn{\mu, \sigma, \alpha} (score) \tab closed form \tab \eqn{10^{-12}}
#'     \cr
#'   \eqn{\nu} (score) \tab [fd5_first()] on \eqn{\ell} \tab \eqn{5\times10^{-11}}
#'     \cr
#'   \eqn{(\mu,\sigma,\alpha)} block (Hessian) \tab closed form \tab
#'     \eqn{10^{-12}} \cr
#'   \eqn{\nu} with another parameter \tab [fd5_first()] on the analytic score
#'     \tab to the printed digit \cr
#'   \eqn{\nu} twice \tab [fd5_second()] on \eqn{\ell} \tab \eqn{2\times10^{-9}}
#'     \cr
#'   \eqn{\nu} three times \tab [fd5_third()] on \eqn{\ell} \tab
#'     \eqn{10^{-4}} \cr
#'   \eqn{\nu} four times \tab [fd5_fourth()] on \eqn{\ell} \tab about one
#'     figure
#' }
#'
#' # The tolerance a fit can ask for
#'
#' The score in \eqn{\nu} cannot be computed more accurately than that table,
#' so no stopping rule on the gradient can be satisfied below it however good
#' the optimizer is. [fit_distrib()] takes its rule from the method it is
#' given, and `crit_grad()`'s default tolerance of \eqn{10^{-6}} is tested on
#' the score **per observation**, which leaves room: measured on samples of
#' 500 to 4000, the run converges with a summed score between
#' \eqn{1.5\times10^{-7}} and \eqn{2.4\times10^{-5}}, i.e. between
#' \eqn{3\times10^{-10}} and \eqn{9\times10^{-9}} per observation.
#'
#' # Fitting
#'
#' The expected information has no closed form, so it is approximated by the
#' strategy named in `approx`, at one quadrature per component.
#' `method = optimizers7::newton()` is much the cheaper route: the observed
#' Hessian is the closed form above and needs no integration.
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
#' # The observed Hessian is the cheap route: this family has no closed-form
#' # expected information, so Fisher scoring would quadrature it at every step.
#' set.seed(1)
#' x <- distrib_rng(d, 200, th)
#' coef(fit_distrib(d, x, method = optimizers7::newton(), start = th))
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
