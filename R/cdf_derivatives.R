#' @include distrib.R generics.R utility_functions.R numerical_functions.R
#' @include gaussian1_distrib.R logistic_distrib.R cauchy_distrib.R laplace_distrib.R
#' @include laplace2_distrib.R
#' @include lognormal1_distrib.R invgauss1_distrib.R poisson_distrib.R binomial_distrib.R
#' @include bernoulli_distrib.R negbin2_distrib.R student_t1_distrib.R pseudohuber_distrib.R
NULL

# ===========================================================================
# Derivatives of the distribution function with respect to the parameters.
#
# The identity that governs everything here is one exchange of derivative and
# integral. The region of integration does not depend on theta, so
#
#     d^I F(q) / F(q) = E[ d^I f / f  |  Y <= q ],
#
# a partial expectation of exactly the quantity the Bartlett lemma expands: at
# first order the conditional mean of the score, at second order the conditional
# mean of l^(ij) + l^(i) l^(j).
#
# Two consequences shape the implementation. For a *discrete* distribution the
# expectation is a finite sum, so the identity is not an approximation and the
# default method uses it directly. For a *continuous* one it is an integral over
# a semi-infinite region, and evaluating it by quadrature is both slower and less
# accurate than differencing the cdf, which for every distribution in the
# catalog is an analytic function; the default there differences the cdf, and
# distributions with a closed form register it.
# ===========================================================================

#' Put CDF Derivatives on the Requested Tail and Scale
#'
#' @description
#' Converts derivatives of the distribution function \eqn{F} into derivatives
#' of whichever tail was asked for, on the natural or the logarithmic scale.
#' Every route to a cdf derivative in this package produces derivatives of
#' \eqn{F} itself, so the `lower.tail` and `log` arguments are handled once,
#' here, and no method has to implement four cases.
#'
#' @details
#' # The two conversions
#'
#' Switching to the upper tail flips the sign, \eqn{S = 1 - F} giving
#' \eqn{\partial^I S = -\partial^I F} at every order. Switching to the log
#' scale divides by the probability, which at second order brings in the
#' familiar correction
#' \deqn{\partial^2 \log P = \frac{\partial^2 P}{P}
#'       - \frac{\partial P}{P}\,\frac{\partial P}{P}.}
#' The two commute, and both are applied to whatever was handed in.
#'
#' # What the caller must supply
#'
#' `dF1` is always read, at second order as well: the correction above needs
#' the first derivatives of the same tail. A caller asking for a Hessian
#' therefore passes both lists, and passes them as derivatives of \eqn{F} on
#' the natural scale, whichever tail the result is wanted on.
#'
#' @section Notation:
#' \eqn{F} is the distribution function, \eqn{S = 1 - F} the survival
#' function, \eqn{P} whichever of the two was asked for, and \eqn{\partial^I}
#' a derivative with respect to a multi-index of parameters.
#'
#' @param distrib An object inheriting from `distrib`. Only its `params` are
#'   read, to name and pair the components.
#' @param Fq The distribution function at the quantile, a numeric vector.
#' @param dF1 A named list of first derivatives of \eqn{F}, one component per
#'   parameter, in the parameter order.
#' @param dF2 A named list of second derivatives of \eqn{F}, keyed as
#'   [hess_names()], or `NULL` (the default) when only the gradient is wanted.
#' @param lower.tail Is the lower tail wanted? A single logical. `TRUE` leaves
#'   the signs alone and `FALSE` flips every one of them.
#' @param log Are derivatives of the log probability wanted? A single logical.
#'   `TRUE` divides by the probability, which returns `-Inf` or `NaN` in a tail
#'   where that probability has underflowed to zero.
#'
#' @return A named list of numeric vectors: one per parameter when `dF2` is
#'   `NULL`, otherwise one per [hess_names()] component. The gradient is not
#'   returned alongside the Hessian.
#'
#' @seealso [distrib_grad_cdf()] and [distrib_hess_cdf()], the generics whose
#'   methods all end here; [loc_scale_grad_cdf()] and [discrete_cdf_deriv()]
#'   for two of the routes that feed it.
#'
#' @keywords internal
cdf_tail_scale <- function(distrib, Fq, dF1, dF2 = NULL, lower.tail, log) {
  params <- distrib@params
  P <- if (lower.tail) Fq else 1 - Fq
  sgn <- if (lower.tail) 1 else -1

  d1 <- lapply(dF1, function(v) sgn * v)
  names(d1) <- params
  if (is.null(dF2)) {
    if (!log) return(d1)
    return(stats::setNames(lapply(params, function(p) d1[[p]] / P), params))
  }

  d2 <- lapply(dF2, function(v) sgn * v)
  names(d2) <- names(dF2)
  if (!log) return(d2)

  pairs <- hess_pairs(params)
  stats::setNames(lapply(names(pairs), function(nm) {
    pr <- pairs[[nm]]
    d2[[nm]] / P - (d1[[pr[1]]] / P) * (d1[[pr[2]]] / P)
  }), names(pairs))
}

# --- the exact route, for discrete distributions ----------------------------
#
# d^I F(q) = sum over the support up to q of  f(y) * (d^I f / f)(y).
# The sum is finite whenever the support has a finite lower bound, which every
# discrete distribution in the package has and the class requires.

#' CDF Derivatives of a Discrete Distribution
#'
#' @description
#' Evaluates \eqn{\partial^I F(q)} for a discrete family as the finite sum
#' \deqn{\partial^I F(q) = \sum_{y \le q} f(y)\,
#'       \frac{\partial^I f}{f}(y).}
#' Nothing is differenced and nothing is integrated: the sum is exact.
#'
#' @details
#' # Where the identity comes from
#'
#' The region of integration does not depend on \eqn{\theta}, so derivative and
#' integral exchange and
#' \deqn{\frac{\partial^I F(q)}{F(q)}
#'       = \mathbb{E}\!\left[\frac{\partial^I f}{f} \;\middle|\; Y \le q\right],}
#' a partial expectation of exactly the quantity the Bartlett lemma expands: the
#' score at first order, and \eqn{\ell^{(ij)} + \ell^{(i)}\ell^{(j)}} at second.
#' For a discrete family that conditional expectation is a finite sum whenever
#' the support has a finite lower bound, which the discrete class requires, so
#' the identity is used as it stands.
#'
#' # A trap for whoever writes a test
#'
#' Checking this against the partial-expectation sum proves nothing: it is the
#' same sum computed twice. A discrete implementation has to be checked against
#' finite differences of the cdf, and a continuous one against the partial
#' expectation, so that the two sides of the comparison share no arithmetic.
#'
#' @section Notation:
#' \eqn{f} is the mass function, \eqn{F} the distribution function,
#' \eqn{\ell = \log f} and \eqn{\partial^I} a derivative with respect to a
#' multi-index of parameters.
#'
#' @param distrib An object inheriting from `discrete_distrib`.
#' @param q A numeric vector of quantiles. Each is handled separately, the
#'   support being walked from the lower bound up to it, so the cost grows with
#'   the largest quantile asked for.
#' @param theta A named list of parameters on the parameter scale.
#' @param order The derivative order, 1 or 2. At order 1 the components are the
#'   parameters and at order 2 they are [hess_names()].
#'
#' @return A named list of numeric vectors, derivatives of \eqn{F} itself on
#'   the natural scale. [cdf_tail_scale()] puts them on the requested tail.
#'
#' @seealso [distrib_grad_cdf.discrete_distrib()], its caller;
#'   [numerical_cdf_deriv()], the continuous route;
#'   [cdf_tail_scale()].
#'
#' @keywords internal
discrete_cdf_deriv <- function(distrib, q, theta, order) {
  params <- distrib@params
  lo <- distrib@bounds[1]
  n <- length(q)
  nms <- if (order == 1L) params else hess_names(params)
  pairs <- if (order == 2L) hess_pairs(params) else NULL

  out <- lapply(nms, function(nm) numeric(n))
  names(out) <- nms

  for (k in seq_len(n)) {
    th_k <- lapply(theta[seq_along(params)],
                   function(v) if (length(v) > 1) v[k] else v)
    if (q[k] < lo) next
    grid <- seq(lo, floor(q[k] + 1e-9))
    fy <- distrib_pdf(distrib, grid, th_k)
    g <- distrib_gradient(distrib, grid, th_k)
    if (order == 1L) {
      for (nm in nms) out[[nm]][k] <- sum(fy * g[[nm]])
    } else {
      h <- distrib_hessian(distrib, grid, th_k)
      for (nm in nms) {
        pr <- pairs[[nm]]
        out[[nm]][k] <- sum(fy * (h[[nm]] + g[[pr[1]]] * g[[pr[2]]]))
      }
    }
  }
  out
}

# --- the numerical route, for continuous distributions ---------------------

#' Finite-Difference CDF Derivatives
#'
#' @description
#' Differentiates the distribution function with respect to the parameters by
#' central differences of [distrib_cdf()]. This is the route a continuous
#' family takes when it registers no closed form, and it is deliberately
#' preferred there to the exact partial-expectation identity: for a continuous
#' family that identity is an integral over a semi-infinite region, while the
#' cdf itself is analytic for every family in the catalog, so differencing it
#' is both cheaper and more accurate.
#'
#' @details
#' # The stencils
#'
#' At order 1 each component is one central difference,
#' \eqn{\{F(\theta + h) - F(\theta - h)\}/(2h)}. At order 2 a diagonal
#' component is the three-point second difference and an off-diagonal one is
#' the four-point mixed stencil
#' \eqn{\{F(+,+) - F(+,-) - F(-,+) + F(-,-)\}/(4h_ih_j)}, which differences
#' two different variables and is therefore one stencil; the package never
#' composes two first differences in the same variable.
#'
#' # The step
#'
#' The relative step is \eqn{\varepsilon^{1/(k+2)}}, which is
#' \eqn{6.1\times10^{-6}} at order 1 and \eqn{1.2\times10^{-4}} at order 2, and
#' it is scaled per component by `pmax(1, abs(theta[[j]]))`. **One step is
#' chosen per observation and not per parameter**: `theta` may be vectorized,
#' and a step read off its first element would be the wrong size everywhere
#' else.
#'
#' # What it costs and what it delivers
#'
#' Measured on a Gaussian at 1000 quantiles, against the closed form the family
#' registers: the relative error is \eqn{6.1\times10^{-11}} at order 1 and
#' \eqn{1.7\times10^{-7}} at order 2, and the calls take 0.29 ms and 0.71 ms
#' against the closed form's 0.16 ms. Of the 42 univariate families, 8 reach
#' this function for their gradient (beta1, beta2, chisq, gamma1, gamma2,
#' gengamma1 and the two von Mises), where the derivative of an incomplete
#' gamma or beta in its shape is hypergeometric or the cdf is itself a
#' quadrature.
#'
#' @section Notation:
#' \eqn{F} is the distribution function, \eqn{\theta} the parameter on its own
#' scale, \eqn{h} the step and \eqn{\varepsilon} the machine epsilon.
#'
#' @param distrib An object inheriting from `distrib`.
#' @param q A numeric vector of quantiles.
#' @param theta A named list of parameters on the parameter scale. Components
#'   may be vectors; the step is then chosen elementwise.
#' @param order The derivative order, 1 (the default) or 2.
#' @param h_rel The relative step. A single number, defaulting to
#'   \eqn{\varepsilon^{1/(\mathrm{order}+2)}}, which balances the stencil's
#'   truncation against its rounding. A step much smaller than the default is
#'   worse, the rounding growing as \eqn{h^{-\mathrm{order}}}.
#' @param which A character vector naming the components to differentiate, or
#'   `NULL` (the default) for all of them: parameter names at order 1 and
#'   [hess_names()] components at order 2. A family with a closed form for part
#'   of the order passes the rest here, so that only those cost cdf
#'   evaluations.
#'
#' @return A named list of numeric vectors, derivatives of \eqn{F} itself and
#'   not of its logarithm, in the order `which` gives or in the enumeration
#'   order when it is `NULL`.
#'
#' @seealso [discrete_cdf_deriv()], the exact route for a discrete family;
#'   [cdf_tail_scale()], which puts the result on the requested tail;
#'   [distrib_grad_cdf()] and [distrib_hess_cdf()].
#'
#' @examples
#' d <- gaussian1_distrib()
#' th <- list(mu = 0.3, sigma = 1.2)
#' q <- c(-1, 0.5, 2)
#'
#' # One central difference per parameter.
#' numerical_cdf_deriv(d, q, th, order = 1)
#'
#' # Against the closed form the Gaussian registers, on the same scale.
#' exact <- distrib_grad_cdf(d, q, th, log = FALSE)
#' max(abs(numerical_cdf_deriv(d, q, th, 1)$mu / exact$mu - 1))
#'
#' # Only the components asked for are differenced.
#' names(numerical_cdf_deriv(d, q, th, order = 2, which = "mu_mu"))
#'
#' @export
numerical_cdf_deriv <- function(distrib, q, theta, order = 1L,
                                h_rel = .Machine$double.eps^(1 / (order + 2)),
                                which = NULL) {
  params <- distrib@params
  p <- length(params)
  keep <- if (is.null(which)) params else which
  bump <- function(j, s, hj) {
    t2 <- theta
    t2[[j]] <- theta[[j]] + s * hj
    t2
  }
  # one step per observation, not per parameter: theta may be vectorized, and a
  # step chosen from its first element would be wrong everywhere else.
  hs <- lapply(seq_len(p), function(j) h_rel * pmax(1, abs(theta[[j]])))

  if (order == 1L) {
    js <- match(keep, params)
    out <- lapply(js, function(j) {
      (distrib_cdf(distrib, q, bump(j, 1, hs[[j]])) -
         distrib_cdf(distrib, q, bump(j, -1, hs[[j]]))) / (2 * hs[[j]])
    })
    return(stats::setNames(out, keep))
  }

  pairs <- hess_pairs(params)
  if (!is.null(which)) pairs <- pairs[intersect(names(pairs), which)]
  stats::setNames(lapply(names(pairs), function(nm) {
    pr <- pairs[[nm]]
    i <- match(pr[1], params); j <- match(pr[2], params)
    if (i == j) {
      (distrib_cdf(distrib, q, bump(i, 1, hs[[i]])) -
         2 * distrib_cdf(distrib, q, theta) +
         distrib_cdf(distrib, q, bump(i, -1, hs[[i]]))) / hs[[i]]^2
    } else {
      shift <- function(a, b) {
        t2 <- theta
        t2[[i]] <- theta[[i]] + a * hs[[i]]
        t2[[j]] <- theta[[j]] + b * hs[[j]]
        t2
      }
      (distrib_cdf(distrib, q, shift(1, 1)) - distrib_cdf(distrib, q, shift(1, -1)) -
         distrib_cdf(distrib, q, shift(-1, 1)) + distrib_cdf(distrib, q, shift(-1, -1))) /
        (4 * hs[[i]] * hs[[j]])
    }
  }), names(pairs))
}


# --- default methods -------------------------------------------------------

#' @title Default Log-CDF Gradient for Continuous Distributions
#' @name distrib_grad_cdf.continuous_distrib
#'
#' @description
#' The fallback for a continuous family that registers no closed form: one
#' central difference of [distrib_cdf()] per parameter, through
#' [numerical_cdf_deriv()], put on the requested tail and scale by
#' [cdf_tail_scale()]. The step is \eqn{\varepsilon^{1/3}} relative, about
#' \eqn{6.1\times10^{-6}}, and the accuracy measured against a family's own
#' closed form is \eqn{6.1\times10^{-11}} relative.
#'
#' @details
#' Eight of the 42 univariate families reach this method: beta1, beta2, chisq,
#' gamma1, gamma2, gengamma1 and the two von Mises. In the first six the
#' derivative wanted is that of an incomplete gamma or beta function with
#' respect to its shape, which is hypergeometric and has no elementary form; in
#' the last two the distribution function is itself a quadrature. Every other
#' continuous family registers a formula.
#'
#' @param distrib A `continuous_distrib` object.
#' @param q A numeric vector of quantiles.
#' @param theta A named list of parameters on the parameter scale. Components
#'   may be vectors, in which case the step is chosen elementwise.
#' @param lower.tail Is the lower tail wanted? A single logical, `TRUE` by
#'   default.
#' @param log Are derivatives of the log probability wanted? A single logical,
#'   `TRUE` by default. Far into a tail the probability underflows and the
#'   result is `-Inf` or `NaN`.
#'
#' @return A named list of numeric vectors, one per parameter, each the length
#'   of `q` recycled against `theta`.
#'
#' @seealso [numerical_cdf_deriv()] for the stencil and its step;
#'   [distrib_grad_cdf.discrete_distrib()], which is exact;
#'   [distrib_hess_cdf.continuous_distrib()] for the second order.
#'
#' @examples
#' # A gamma reaches this method: its cdf derivative in the shape has no
#' # elementary form, so the cdf is differenced instead.
#' d <- gamma2_distrib()
#' distrib_grad_cdf(d, c(1, 2, 4), list(mu = 2, sigma2 = 1))
#'
#' @keywords internal
S7::method(distrib_grad_cdf, continuous_distrib) <- function(distrib, q, theta,
                                                             lower.tail = TRUE, log = TRUE) {
  cdf_tail_scale(distrib, distrib_cdf(distrib, q, theta),
                 numerical_cdf_deriv(distrib, q, theta, 1L),
                 NULL, lower.tail, log)
}

#' @title Default Log-CDF Hessian for Continuous Distributions
#' @name distrib_hess_cdf.continuous_distrib
#'
#' @description
#' The fallback for a continuous family that registers no closed second
#' derivative: [numerical_cdf_deriv()] at order 2, with the first-order part
#' taken from [distrib_grad_cdf()], so a closed gradient the family does
#' register is used here too. The step is
#' \eqn{\varepsilon^{1/4}} relative, about \eqn{1.2\times10^{-4}}, and the
#' accuracy measured against a family's own closed form is
#' \eqn{1.7\times10^{-7}} relative.
#'
#' @details
#' A diagonal component is the three-point second difference and an
#' off-diagonal one the four-point mixed stencil, which differences two
#' different variables and is therefore a single stencil. The log-scale
#' correction \eqn{\partial^2 P/P - (\partial P/P)^2} is applied afterwards by
#' [cdf_tail_scale()], which is why the gradient is fetched even when only the
#' Hessian was asked for.
#'
#' @param distrib A `continuous_distrib` object.
#' @param q A numeric vector of quantiles.
#' @param theta A named list of parameters on the parameter scale.
#' @param lower.tail Is the lower tail wanted? A single logical, `TRUE` by
#'   default.
#' @param log Are derivatives of the log probability wanted? A single logical,
#'   `TRUE` by default.
#'
#' @return A named list of numeric vectors keyed as [hess_names()], each the
#'   length of `q` recycled against `theta`. The gradient is not returned
#'   alongside.
#'
#' @seealso [distrib_grad_cdf.continuous_distrib()] for the first order;
#'   [numerical_cdf_deriv()] for the stencils; [cdf_tail_scale()].
#'
#' @examples
#' # A gamma reaches this method at both orders.
#' d <- gamma2_distrib()
#' distrib_hess_cdf(d, c(1, 2), list(mu = 2, sigma2 = 1))
#'
#' @keywords internal
S7::method(distrib_hess_cdf, continuous_distrib) <- function(distrib, q, theta,
                                                             lower.tail = TRUE, log = TRUE) {
  # the first-order part comes from distrib_grad_cdf, so a closed form registered
  # for this family is used here too rather than differenced a second time
  cdf_tail_scale(distrib, distrib_cdf(distrib, q, theta),
                 distrib_grad_cdf(distrib, q, theta, lower.tail = TRUE, log = FALSE),
                 numerical_cdf_deriv(distrib, q, theta, 2L),
                 lower.tail, log)
}

#' @title Log-CDF Gradient for Discrete Distributions
#' @name distrib_grad_cdf.discrete_distrib
#'
#' @description
#' Exact, and nothing is differenced. The partial expectation of the score is a
#' finite sum over the support up to \eqn{q},
#' \eqn{\partial^I F(q) = \sum_{y \le q} f(y)\,\partial^I f/f\,(y)}, evaluated
#' by [discrete_cdf_deriv()] and put on the requested tail by
#' [cdf_tail_scale()].
#'
#' @details
#' The sum is finite because the discrete class requires a finite lower bound
#' on the support. The cost therefore grows with the largest quantile asked
#' for and not with the number of parameters: on a Poisson of mean 30 at 200
#' quantiles the whole gradient takes 0.14 ms.
#'
#' Six of the shipped discrete families use this method; the rest register a
#' formula, several of which are exact simplifications of this same sum. A
#' Poisson's, for instance, telescopes to \eqn{\partial F(k)/\partial\mu =
#' -f(k)}.
#'
#' @param distrib A `discrete_distrib` object.
#' @param q A numeric vector of quantiles. Values below the support give a
#'   derivative of zero, the sum being empty.
#' @param theta A named list of parameters on the parameter scale.
#' @param lower.tail Is the lower tail wanted? A single logical, `TRUE` by
#'   default.
#' @param log Are derivatives of the log probability wanted? A single logical,
#'   `TRUE` by default.
#'
#' @return A named list of numeric vectors, one per parameter, each the length
#'   of `q` recycled against `theta`.
#'
#' @seealso [discrete_cdf_deriv()] for the sum;
#'   [distrib_grad_cdf.PoissonDistrib()] for a family that simplifies it;
#'   [distrib_hess_cdf.discrete_distrib()] for the second order.
#'
#' @examples
#' # A beta-binomial reaches this method, and the sum is exact.
#' d <- betabinom1_distrib(size = 10)
#' distrib_grad_cdf(d, c(2, 5, 8), list(mu = 0.3, sigma = 0.5))
#'
#' @keywords internal
S7::method(distrib_grad_cdf, discrete_distrib) <- function(distrib, q, theta,
                                                           lower.tail = TRUE, log = TRUE) {
  cdf_tail_scale(distrib, distrib_cdf(distrib, q, theta),
                 discrete_cdf_deriv(distrib, q, theta, 1L),
                 NULL, lower.tail, log)
}

#' @title Log-CDF Hessian for Discrete Distributions
#' @name distrib_hess_cdf.discrete_distrib
#'
#' @description
#' Exact, by the same finite sum as the gradient: at second order the summand
#' is \eqn{\ell^{(ij)} + \ell^{(i)}\ell^{(j)}}, the Bartlett lemma's expansion
#' of \eqn{\partial^2 f/f}. The first-order part comes from
#' [distrib_grad_cdf()], so a family's own simplification is used where it has
#' one.
#'
#' @param distrib A `discrete_distrib` object.
#' @param q A numeric vector of quantiles. Values below the support give zero.
#' @param theta A named list of parameters on the parameter scale.
#' @param lower.tail Is the lower tail wanted? A single logical, `TRUE` by
#'   default.
#' @param log Are derivatives of the log probability wanted? A single logical,
#'   `TRUE` by default.
#'
#' @return A named list of numeric vectors keyed as [hess_names()], each the
#'   length of `q` recycled against `theta`. The gradient is not returned
#'   alongside.
#'
#' @seealso [distrib_grad_cdf.discrete_distrib()] for the first order;
#'   [discrete_cdf_deriv()] for the sum; [cdf_tail_scale()].
#'
#' @examples
#' # A negative binomial: two parameters, so three Hessian components.
#' d <- negbin2_distrib()
#' distrib_hess_cdf(d, c(2, 5), list(mu = 3, theta = 2))
#'
#' @keywords internal
S7::method(distrib_hess_cdf, discrete_distrib) <- function(distrib, q, theta,
                                                           lower.tail = TRUE, log = TRUE) {
  cdf_tail_scale(distrib, distrib_cdf(distrib, q, theta),
                 distrib_grad_cdf(distrib, q, theta, lower.tail = TRUE, log = FALSE),
                 discrete_cdf_deriv(distrib, q, theta, 2L),
                 lower.tail, log)
}


# --- closed forms for location-scale families ------------------------------
#
# If F(y) = G(z) with z = (y - mu)/sigma, then differentiating the composition
# gives everything in terms of the density and its derivative in y:
#
#   dF/dmu     = -f
#   dF/dsigma  = -z f
#   d2F/dmu2   = f * l_y                        l_y = d log f / dy
#   d2F/dmu dsigma = f * (z l_y + 1/sigma)
#   d2F/dsigma2    = f * (z^2 l_y + 2z/sigma)
#
# so a family only has to say that it is location-scale. This covers the
# censored-regression workhorses -- Gaussian, logistic, Cauchy, Laplace -- and
# needs nothing from the family beyond its density and response derivative,
# both of which every distribution already provides.

#' CDF Derivatives of a Location-Scale Family
#'
#' @description
#' Closed-form derivatives of \eqn{F} for a family that is location-scale in
#' its first two parameters. Writing \eqn{F(q) = G(z)} with
#' \eqn{z = (q-\mu)/\sigma} and differentiating the composition puts everything
#' in terms of the density and its response derivative, both of which every
#' distribution already supplies, so a family has only to declare that it is
#' location-scale.
#'
#' @details
#' With \eqn{\ell_y = \partial \log f/\partial y},
#' \deqn{\frac{\partial F}{\partial \mu} = -f, \qquad
#'       \frac{\partial F}{\partial \sigma} = -z f,}
#' \deqn{\frac{\partial^2 F}{\partial \mu^2} = f\,\ell_y, \qquad
#'       \frac{\partial^2 F}{\partial \mu \partial \sigma}
#'         = f\left(z\,\ell_y + \frac{1}{\sigma}\right), \qquad
#'       \frac{\partial^2 F}{\partial \sigma^2}
#'         = f\left(z^2 \ell_y + \frac{2z}{\sigma}\right).}
#'
#' This covers the censored-regression workhorses: the Gaussian, the logistic,
#' the Cauchy and the Laplace all register these bodies unchanged, and the
#' Student t and the pseudo-Huber use them for their first two components and
#' difference the shape.
#'
#' @section Notation:
#' \eqn{\mu} is the location, \eqn{\sigma > 0} the scale,
#' \eqn{z = (q-\mu)/\sigma}, \eqn{f} the density and
#' \eqn{\ell_y = \partial\log f/\partial y} its response derivative.
#'
#' @param distrib An object inheriting from `distrib` whose first two
#'   parameters are a location and a scale. Nothing checks that; a family that
#'   is not location-scale and registers this body gets wrong numbers.
#' @param q A numeric vector of quantiles.
#' @param theta A named list of parameters, the location first and the scale
#'   second. Any further parameters are read by [distrib_pdf()] and
#'   [distrib_grad_y()] but get no component here.
#' @param order The derivative order, 1 or 2.
#'
#' @return A named list of numeric vectors, derivatives of \eqn{F} itself on
#'   the natural scale: two components at order 1, and the three of
#'   [hess_names()] at order 2.
#'
#' @seealso [loc_scale_grad_cdf()] and [loc_scale_hess_cdf()], the two bodies
#'   families register; [partial_loc_scale_grad_cdf()], the variant for a
#'   family with a shape parameter as well; [cdf_tail_scale()].
#'
#' @keywords internal
loc_scale_cdf_deriv <- function(distrib, q, theta, order) {
  mu <- theta[[1]]
  s <- theta[[2]]
  z <- (q - mu) / s
  f <- distrib_pdf(distrib, q, theta)
  if (order == 1L) {
    return(stats::setNames(list(-f, -z * f), distrib@params))
  }
  ly <- distrib_grad_y(distrib, q, theta)
  vals <- list(f * ly, f * (z^2 * ly + 2 * z / s), f * (z * ly + 1 / s))
  stats::setNames(vals, hess_names(distrib@params))
}

#' Location-Scale CDF Gradient
#'
#' @description
#' The [distrib_grad_cdf()] body the location-scale families share:
#' [loc_scale_cdf_deriv()] at order 1, put on the requested tail and scale by
#' [cdf_tail_scale()]. The Gaussian, the logistic, the Cauchy and the Laplace
#' register this function itself, so their four pages document one piece of
#' arithmetic.
#'
#' @param distrib An object inheriting from `distrib` whose first two
#'   parameters are a location and a scale.
#' @param q A numeric vector of quantiles.
#' @param theta A named list of parameters, the location first and the scale
#'   second.
#' @param lower.tail Is the lower tail wanted? A single logical, `TRUE` by
#'   default.
#' @param log Are derivatives of the log probability wanted? A single logical,
#'   `TRUE` by default.
#'
#' @return A named list of two numeric vectors, one per parameter, each the
#'   length of `q` recycled against `theta`.
#'
#' @seealso [loc_scale_cdf_deriv()] for the formulas;
#'   [loc_scale_hess_cdf()] for the second order;
#'   [distrib_grad_cdf.Gaussian1Distrib()] for one of the four registrations.
#'
#' @keywords internal
loc_scale_grad_cdf <- function(distrib, q, theta, lower.tail = TRUE, log = TRUE) {
  cdf_tail_scale(distrib, distrib_cdf(distrib, q, theta),
                 loc_scale_cdf_deriv(distrib, q, theta, 1L), NULL, lower.tail, log)
}

#' Location-Scale CDF Hessian
#'
#' @description
#' The [distrib_hess_cdf()] body the location-scale families share:
#' [loc_scale_cdf_deriv()] at both orders, put on the requested tail and scale
#' by [cdf_tail_scale()]. Both orders are needed even when only the Hessian was
#' asked for, the log-scale correction reading the first derivatives.
#'
#' @param distrib An object inheriting from `distrib` whose first two
#'   parameters are a location and a scale.
#' @param q A numeric vector of quantiles.
#' @param theta A named list of parameters, the location first and the scale
#'   second.
#' @param lower.tail Is the lower tail wanted? A single logical, `TRUE` by
#'   default.
#' @param log Are derivatives of the log probability wanted? A single logical,
#'   `TRUE` by default.
#'
#' @return A named list of three numeric vectors keyed as [hess_names()], each
#'   the length of `q` recycled against `theta`. The gradient is not returned
#'   alongside.
#'
#' @seealso [loc_scale_cdf_deriv()] for the formulas;
#'   [loc_scale_grad_cdf()] for the first order;
#'   [distrib_hess_cdf.Gaussian1Distrib()] for one of the four registrations.
#'
#' @keywords internal
loc_scale_hess_cdf <- function(distrib, q, theta, lower.tail = TRUE, log = TRUE) {
  cdf_tail_scale(distrib, distrib_cdf(distrib, q, theta),
                 loc_scale_cdf_deriv(distrib, q, theta, 1L),
                 loc_scale_cdf_deriv(distrib, q, theta, 2L), lower.tail, log)
}

#' @title Gaussian Log-CDF Gradient
#' @name distrib_grad_cdf.Gaussian1Distrib
#'
#' @description
#' Closed form, from the location-scale structure:
#' \eqn{\partial F/\partial\mu = -f(q)} and
#' \eqn{\partial F/\partial\sigma = -z f(q)} with \eqn{z = (q-\mu)/\sigma}. The
#' method is [loc_scale_grad_cdf()] itself, shared with the logistic, the
#' Cauchy and the Laplace, so nothing here is particular to the Gaussian beyond
#' its density.
#'
#' @details
#' A censored Gaussian likelihood is cheap for this reason: a right-censored
#' observation contributes \eqn{\log S(q)}, whose score in \eqn{\mu} is
#' \eqn{f(q)/S(q)}, the inverse Mills ratio, and the whole of it comes from the
#' density the family already computes.
#'
#' @section Notation:
#' \eqn{\mu} is the mean, \eqn{\sigma > 0} the standard deviation,
#' \eqn{z = (q-\mu)/\sigma}, \eqn{f} the density and \eqn{F} the distribution
#' function.
#'
#' @param distrib A `Gaussian1Distrib` object, from [gaussian1_distrib()].
#' @param q A numeric vector of quantiles.
#' @param theta A named list with components `mu` and `sigma` (positive), each
#'   a numeric vector of length 1 or `n`.
#' @param lower.tail Is the lower tail wanted? A single logical, `TRUE` by
#'   default; `FALSE` flips the sign of every component.
#' @param log Are derivatives of the log probability wanted? A single logical,
#'   `TRUE` by default. Far into a tail the probability underflows to zero and
#'   the result is `-Inf` or `NaN`.
#'
#' @return A named list of two numeric vectors, `mu` and `sigma`, each the
#'   length of `q` recycled against `theta`.
#'
#' @seealso [loc_scale_grad_cdf()] for the shared body;
#'   [distrib_hess_cdf.Gaussian1Distrib()] for the second order;
#'   [gaussian1_distrib()].
#'
#' @examples
#' d <- gaussian1_distrib()
#' th <- list(mu = 0.3, sigma = 1.2)
#' q <- c(-1, 0.5, 2)
#'
#' # On the natural scale the mean component is minus the density.
#' all.equal(distrib_grad_cdf(d, q, th, log = FALSE)$mu,
#'           -distrib_pdf(d, q, th))
#'
#' # On the log scale it is -f/F, and the upper tail flips to +f/S.
#' distrib_grad_cdf(d, q, th)$mu
#' distrib_grad_cdf(d, q, th, lower.tail = FALSE)$mu
S7::method(distrib_grad_cdf, Gaussian1Distrib) <- loc_scale_grad_cdf

#' @title Gaussian Log-CDF Hessian
#' @name distrib_hess_cdf.Gaussian1Distrib
#'
#' @description
#' Closed form, from the same location-scale structure. With
#' \eqn{\ell_y = \partial\log f/\partial y}, which for a Gaussian is
#' \eqn{-z/\sigma},
#' \deqn{\frac{\partial^2 F}{\partial\mu^2} = f\,\ell_y, \qquad
#'       \frac{\partial^2 F}{\partial\mu\,\partial\sigma}
#'         = f\left(z\ell_y + \frac{1}{\sigma}\right), \qquad
#'       \frac{\partial^2 F}{\partial\sigma^2}
#'         = f\left(z^2\ell_y + \frac{2z}{\sigma}\right).}
#' The method is [loc_scale_hess_cdf()] itself, shared with the logistic, the
#' Cauchy and the Laplace.
#'
#' @section Notation:
#' \eqn{\mu} is the mean, \eqn{\sigma > 0} the standard deviation,
#' \eqn{z = (q-\mu)/\sigma}, \eqn{f} the density and
#' \eqn{\ell_y = \partial\log f/\partial y} its response derivative.
#'
#' @param distrib A `Gaussian1Distrib` object, from [gaussian1_distrib()].
#' @param q A numeric vector of quantiles.
#' @param theta A named list with components `mu` and `sigma` (positive), each
#'   a numeric vector of length 1 or `n`.
#' @param lower.tail Is the lower tail wanted? A single logical, `TRUE` by
#'   default.
#' @param log Are derivatives of the log probability wanted? A single logical,
#'   `TRUE` by default.
#'
#' @return A named list of three numeric vectors keyed as [hess_names()],
#'   `mu_mu`, `sigma_sigma` and `mu_sigma`, each the length of `q` recycled
#'   against `theta`. The gradient is not returned alongside.
#'
#' @seealso [loc_scale_hess_cdf()] for the shared body;
#'   [distrib_grad_cdf.Gaussian1Distrib()] for the first order;
#'   [gaussian1_distrib()].
#'
#' @examples
#' d <- gaussian1_distrib()
#' th <- list(mu = 0.3, sigma = 1.2)
#' q <- c(-1, 0.5, 2)
#'
#' distrib_hess_cdf(d, q, th)
#'
#' # Against a central difference of the cdf, which shares no arithmetic.
#' exact <- distrib_hess_cdf(d, q, th, log = FALSE)
#' fd <- numerical_cdf_deriv(d, q, th, order = 2)
#' max(abs(unlist(exact[names(fd)]) / unlist(fd) - 1))
S7::method(distrib_hess_cdf, Gaussian1Distrib) <- loc_scale_hess_cdf

#' @title Logistic Log-CDF Gradient
#' @name distrib_grad_cdf.LogisticDistrib
#'
#' @description
#' Closed form, from the location-scale structure:
#' \eqn{\partial F/\partial\mu = -f(q)} and
#' \eqn{\partial F/\partial\sigma = -z f(q)} with \eqn{z = (q-\mu)/\sigma}. The
#' method is [loc_scale_grad_cdf()] itself. For this family the density is
#' \eqn{F(1-F)/\sigma}, so the log-scale gradient is elementary in \eqn{F}
#' alone and no special evaluation is needed anywhere on the line.
#'
#' @section Notation:
#' \eqn{\mu} is the location, \eqn{\sigma > 0} the scale,
#' \eqn{z = (q-\mu)/\sigma}, \eqn{f} the density and \eqn{F} the distribution
#' function.
#'
#' @param distrib A `LogisticDistrib` object, from [logistic_distrib()].
#' @param q A numeric vector of quantiles.
#' @param theta A named list with components `mu` and `sigma` (positive), each
#'   a numeric vector of length 1 or `n`.
#' @param lower.tail Is the lower tail wanted? A single logical, `TRUE` by
#'   default.
#' @param log Are derivatives of the log probability wanted? A single logical,
#'   `TRUE` by default.
#'
#' @return A named list of two numeric vectors, `mu` and `sigma`, each the
#'   length of `q` recycled against `theta`.
#'
#' @seealso [loc_scale_grad_cdf()] for the shared body;
#'   [distrib_hess_cdf.LogisticDistrib()] for the second order;
#'   [logistic_distrib()].
#'
#' @examples
#' d <- logistic_distrib()
#' th <- list(mu = 0.3, sigma = 1.2)
#' q <- c(-1, 0.5, 2)
#'
#' # The mean component on the log scale is -f/F, which here is -(1-F)/sigma.
#' Fq <- distrib_cdf(d, q, th)
#' all.equal(distrib_grad_cdf(d, q, th)$mu, -(1 - Fq) / 1.2)
S7::method(distrib_grad_cdf, LogisticDistrib) <- loc_scale_grad_cdf

#' @title Logistic Log-CDF Hessian
#' @name distrib_hess_cdf.LogisticDistrib
#'
#' @description
#' Closed form, from the same location-scale structure, through
#' [loc_scale_hess_cdf()]. With \eqn{\ell_y = \partial\log f/\partial y}, which
#' for a logistic is \eqn{(1-2F)/\sigma},
#' \deqn{\frac{\partial^2 F}{\partial\mu^2} = f\,\ell_y, \qquad
#'       \frac{\partial^2 F}{\partial\mu\,\partial\sigma}
#'         = f\left(z\ell_y + \frac{1}{\sigma}\right), \qquad
#'       \frac{\partial^2 F}{\partial\sigma^2}
#'         = f\left(z^2\ell_y + \frac{2z}{\sigma}\right).}
#'
#' @section Notation:
#' \eqn{\mu} is the location, \eqn{\sigma > 0} the scale,
#' \eqn{z = (q-\mu)/\sigma}, \eqn{f} the density and
#' \eqn{\ell_y = \partial\log f/\partial y}.
#'
#' @param distrib A `LogisticDistrib` object, from [logistic_distrib()].
#' @param q A numeric vector of quantiles.
#' @param theta A named list with components `mu` and `sigma` (positive), each
#'   a numeric vector of length 1 or `n`.
#' @param lower.tail Is the lower tail wanted? A single logical, `TRUE` by
#'   default.
#' @param log Are derivatives of the log probability wanted? A single logical,
#'   `TRUE` by default.
#'
#' @return A named list of three numeric vectors keyed as [hess_names()], each
#'   the length of `q` recycled against `theta`.
#'
#' @seealso [loc_scale_hess_cdf()] for the shared body;
#'   [distrib_grad_cdf.LogisticDistrib()] for the first order;
#'   [logistic_distrib()].
#'
#' @examples
#' d <- logistic_distrib()
#' th <- list(mu = 0.3, sigma = 1.2)
#' q <- c(-1, 0.5, 2)
#'
#' # Against a central difference of the cdf, which shares no arithmetic.
#' exact <- distrib_hess_cdf(d, q, th, log = FALSE)
#' fd <- numerical_cdf_deriv(d, q, th, order = 2)
#' max(abs(unlist(exact[names(fd)]) / unlist(fd) - 1))
S7::method(distrib_hess_cdf, LogisticDistrib) <- loc_scale_hess_cdf

#' @title Cauchy Log-CDF Gradient
#' @name distrib_grad_cdf.CauchyDistrib
#'
#' @description
#' Closed form, from the location-scale structure:
#' \eqn{\partial F/\partial\mu = -f(q)} and
#' \eqn{\partial F/\partial\sigma = -z f(q)} with \eqn{z = (q-\mu)/\sigma}. The
#' method is [loc_scale_grad_cdf()] itself.
#'
#' @details
#' No moment of the Cauchy exists, and none is needed here: the distribution
#' function is elementary, \eqn{F(q) = 1/2 + \arctan(z)/\pi}, and the identity
#' \eqn{\partial F/\partial\mu = -f} holds for every location-scale family
#' whether or not its moments converge. The quantities this page returns are
#' therefore exact where [mean.CauchyDistrib()] and its siblings return `NaN`.
#'
#' @section Notation:
#' \eqn{\mu} is the location, \eqn{\sigma > 0} the scale,
#' \eqn{z = (q-\mu)/\sigma}, \eqn{f} the density and \eqn{F} the distribution
#' function. Neither \eqn{\mu} nor \eqn{\sigma} is a moment: they are the
#' median and the half-interquartile range.
#'
#' @param distrib A `CauchyDistrib` object, from [cauchy_distrib()].
#' @param q A numeric vector of quantiles.
#' @param theta A named list with components `mu` and `sigma` (positive), each
#'   a numeric vector of length 1 or `n`.
#' @param lower.tail Is the lower tail wanted? A single logical, `TRUE` by
#'   default.
#' @param log Are derivatives of the log probability wanted? A single logical,
#'   `TRUE` by default. The Cauchy's tails are heavy, so the probability
#'   underflows far later here than for a Gaussian.
#'
#' @return A named list of two numeric vectors, `mu` and `sigma`, each the
#'   length of `q` recycled against `theta`.
#'
#' @seealso [loc_scale_grad_cdf()] for the shared body;
#'   [distrib_hess_cdf.CauchyDistrib()] for the second order;
#'   [mean.CauchyDistrib()] for the moments, which do not exist;
#'   [cauchy_distrib()].
#'
#' @examples
#' d <- cauchy_distrib()
#' th <- list(mu = 0.3, sigma = 1.2)
#' q <- c(-1, 0.5, 2)
#'
#' # Exact, even though no moment of this family exists.
#' all.equal(distrib_grad_cdf(d, q, th, log = FALSE)$mu,
#'           -distrib_pdf(d, q, th))
#'
#' # The heavy tail keeps the log-scale gradient finite far out.
#' distrib_grad_cdf(d, -1000, th)$mu
S7::method(distrib_grad_cdf, CauchyDistrib) <- loc_scale_grad_cdf

#' @title Cauchy Log-CDF Hessian
#' @name distrib_hess_cdf.CauchyDistrib
#'
#' @description
#' Closed form, from the same location-scale structure, through
#' [loc_scale_hess_cdf()]. With \eqn{\ell_y = \partial\log f/\partial y},
#' which for a Cauchy is \eqn{-2z/\{\sigma(1+z^2)\}},
#' \deqn{\frac{\partial^2 F}{\partial\mu^2} = f\,\ell_y, \qquad
#'       \frac{\partial^2 F}{\partial\mu\,\partial\sigma}
#'         = f\left(z\ell_y + \frac{1}{\sigma}\right), \qquad
#'       \frac{\partial^2 F}{\partial\sigma^2}
#'         = f\left(z^2\ell_y + \frac{2z}{\sigma}\right).}
#' As with the gradient, these are exact although the family has no moments.
#'
#' @section Notation:
#' \eqn{\mu} is the location, \eqn{\sigma > 0} the scale,
#' \eqn{z = (q-\mu)/\sigma}, \eqn{f} the density and
#' \eqn{\ell_y = \partial\log f/\partial y}.
#'
#' @param distrib A `CauchyDistrib` object, from [cauchy_distrib()].
#' @param q A numeric vector of quantiles.
#' @param theta A named list with components `mu` and `sigma` (positive), each
#'   a numeric vector of length 1 or `n`.
#' @param lower.tail Is the lower tail wanted? A single logical, `TRUE` by
#'   default.
#' @param log Are derivatives of the log probability wanted? A single logical,
#'   `TRUE` by default.
#'
#' @return A named list of three numeric vectors keyed as [hess_names()], each
#'   the length of `q` recycled against `theta`.
#'
#' @seealso [loc_scale_hess_cdf()] for the shared body;
#'   [distrib_grad_cdf.CauchyDistrib()] for the first order;
#'   [cauchy_distrib()].
#'
#' @examples
#' d <- cauchy_distrib()
#' th <- list(mu = 0.3, sigma = 1.2)
#' q <- c(-1, 0.5, 2)
#'
#' # Against a central difference of the cdf, which shares no arithmetic.
#' exact <- distrib_hess_cdf(d, q, th, log = FALSE)
#' fd <- numerical_cdf_deriv(d, q, th, order = 2)
#' max(abs(unlist(exact[names(fd)]) / unlist(fd) - 1))
S7::method(distrib_hess_cdf, CauchyDistrib) <- loc_scale_hess_cdf

#' @title Laplace Log-CDF Gradient
#' @name distrib_grad_cdf.LaplaceDistrib
#'
#' @description
#' Closed form, from the location-scale structure:
#' \eqn{\partial F/\partial\mu = -f(q)} and
#' \eqn{\partial F/\partial\sigma = -z f(q)} with \eqn{z = (q-\mu)/\sigma}. The
#' method is [loc_scale_grad_cdf()] itself. Both are continuous everywhere,
#' including at \eqn{q = \mu}: the density has a kink there but no jump, and
#' the first derivatives of \eqn{F} read the density itself.
#'
#' @section Notation:
#' \eqn{\mu} is the location, \eqn{\sigma > 0} the scale,
#' \eqn{z = (q-\mu)/\sigma}, \eqn{f} the density and \eqn{F} the distribution
#' function. The variance is \eqn{2\sigma^2}, so the scale is not a standard
#' deviation.
#'
#' @param distrib A `LaplaceDistrib` object, from [laplace_distrib()].
#' @param q A numeric vector of quantiles.
#' @param theta A named list with components `mu` and `sigma` (positive), each
#'   a numeric vector of length 1 or `n`.
#' @param lower.tail Is the lower tail wanted? A single logical, `TRUE` by
#'   default.
#' @param log Are derivatives of the log probability wanted? A single logical,
#'   `TRUE` by default.
#'
#' @return A named list of two numeric vectors, `mu` and `sigma`, each the
#'   length of `q` recycled against `theta`.
#'
#' @seealso [loc_scale_grad_cdf()] for the shared body;
#'   [distrib_hess_cdf.LaplaceDistrib()], where the kink does show;
#'   [distrib_grad_cdf.Laplace2Distrib()] for the rate parametrization;
#'   [laplace_distrib()].
#'
#' @examples
#' d <- laplace_distrib()
#' th <- list(mu = 0.3, sigma = 1.2)
#' q <- c(-1, 0.5, 2)
#'
#' all.equal(distrib_grad_cdf(d, q, th, log = FALSE)$mu,
#'           -distrib_pdf(d, q, th))
#'
#' # Continuous through the kink at q = mu.
#' distrib_grad_cdf(d, 0.3 + c(-1e-8, 0, 1e-8), th, log = FALSE)$mu
S7::method(distrib_grad_cdf, LaplaceDistrib) <- loc_scale_grad_cdf

#' @title Laplace Log-CDF Hessian
#' @name distrib_hess_cdf.LaplaceDistrib
#'
#' @description
#' Closed form, from the same location-scale structure, through
#' [loc_scale_hess_cdf()]. The response derivative it reads is
#' \eqn{\ell_y = -\mathrm{sign}(z)/\sigma}, which jumps at \eqn{q = \mu}, so
#' the second derivatives of \eqn{F} jump there too. That is a property of the
#' law and not a defect: the Laplace is the toolkit's non-regular family, and
#' its location has no second derivative at the kink.
#'
#' @details
#' Away from \eqn{q = \mu} the three components are the ordinary
#' location-scale ones. At \eqn{q = \mu} exactly, `sign(0)` is 0 and the
#' returned `mu_mu` is the average of the two one-sided limits, which is the
#' value the sign convention gives and is reported as it stands.
#'
#' @section Notation:
#' \eqn{\mu} is the location, \eqn{\sigma > 0} the scale,
#' \eqn{z = (q-\mu)/\sigma}, \eqn{f} the density and
#' \eqn{\ell_y = \partial\log f/\partial y}.
#'
#' @param distrib A `LaplaceDistrib` object, from [laplace_distrib()].
#' @param q A numeric vector of quantiles.
#' @param theta A named list with components `mu` and `sigma` (positive), each
#'   a numeric vector of length 1 or `n`.
#' @param lower.tail Is the lower tail wanted? A single logical, `TRUE` by
#'   default.
#' @param log Are derivatives of the log probability wanted? A single logical,
#'   `TRUE` by default.
#'
#' @return A named list of three numeric vectors keyed as [hess_names()], each
#'   the length of `q` recycled against `theta`.
#'
#' @seealso [loc_scale_hess_cdf()] for the shared body;
#'   [distrib_grad_cdf.LaplaceDistrib()], which is continuous at the kink;
#'   [laplace_distrib()].
#'
#' @examples
#' d <- laplace_distrib()
#' th <- list(mu = 0.3, sigma = 1.2)
#'
#' # The second derivative in the location jumps across q = mu.
#' distrib_hess_cdf(d, 0.3 + c(-1e-6, 1e-6), th, log = FALSE)$mu_mu
#'
#' # Away from the kink it agrees with a central difference of the cdf.
#' q <- c(-1, 2)
#' exact <- distrib_hess_cdf(d, q, th, log = FALSE)
#' fd <- numerical_cdf_deriv(d, q, th, order = 2)
#' max(abs(unlist(exact[names(fd)]) / unlist(fd) - 1))
S7::method(distrib_hess_cdf, LaplaceDistrib) <- loc_scale_hess_cdf

#' @title Laplace Log-CDF Gradient in Location and Rate
#' @name distrib_grad_cdf.Laplace2Distrib
#'
#' @description
#' Closed form: \eqn{\partial F/\partial\mu = -f(q)} and
#' \eqn{\partial F/\partial\lambda = (q-\mu)\,f(q)/\lambda}. This is the same
#' law as [laplace_distrib()] written by its rate \eqn{\lambda = 1/\sigma}, and
#' the chain rule turns the scale component into the one above.
#'
#' @section Notation:
#' \eqn{\mu} is the location, \eqn{\lambda > 0} the rate, \eqn{f} the density
#' and \eqn{F} the distribution function. The variance is \eqn{2/\lambda^2}.
#'
#' @param distrib A `Laplace2Distrib` object, from [laplace2_distrib()].
#' @param q A numeric vector of quantiles.
#' @param theta A named list with components `mu` and `lambda` (positive), each
#'   a numeric vector of length 1 or `n`.
#' @param lower.tail Is the lower tail wanted? A single logical, `TRUE` by
#'   default.
#' @param log Are derivatives of the log probability wanted? A single logical,
#'   `TRUE` by default.
#'
#' @return A named list of two numeric vectors, `mu` and `lambda`, each the
#'   length of `q` recycled against `theta`.
#'
#' @seealso [distrib_hess_cdf.Laplace2Distrib()], where the kink shows;
#'   [distrib_grad_cdf.LaplaceDistrib()] for the scale parametrization;
#'   [laplace2_distrib()].
#'
#' @examples
#' d <- laplace2_distrib()
#' th <- list(mu = 0.3, lambda = 1 / 1.2)
#' q <- c(-1, 0.5, 2)
#'
#' g <- distrib_grad_cdf(d, q, th, log = FALSE)
#' all.equal(g$mu, -distrib_pdf(d, q, th))
#' all.equal(g$lambda, (q - 0.3) * distrib_pdf(d, q, th) / (1 / 1.2))
S7::method(distrib_grad_cdf, Laplace2Distrib) <- function(distrib, q, theta,
                                                          lower.tail = TRUE, log = TRUE) {
  f <- distrib_pdf(distrib, q, theta)
  d1 <- stats::setNames(list(-f, (q - theta[[1]]) * f / theta[[2]]), distrib@params)
  cdf_tail_scale(distrib, distrib_cdf(distrib, q, theta), d1, NULL, lower.tail, log)
}

#' @title Laplace Log-CDF Hessian in Location and Rate
#' @name distrib_hess_cdf.Laplace2Distrib
#'
#' @description
#' Closed form. Writing \eqn{r = q - \mu} and \eqn{a = |r|},
#' \deqn{\frac{\partial^2 F}{\partial\mu^2} = -\lambda\,\mathrm{sign}(r)\,f,
#'       \qquad
#'       \frac{\partial^2 F}{\partial\lambda^2}
#'         = -\frac{\mathrm{sign}(r)\,a^2 f}{\lambda},
#'       \qquad
#'       \frac{\partial^2 F}{\partial\mu\,\partial\lambda}
#'         = -\left(\frac{1}{\lambda} - a\right) f.}
#' The sign function is what carries the kink: all three jump at
#' \eqn{q = \mu}, as they do in the scale parametrization.
#'
#' @section Notation:
#' \eqn{\mu} is the location, \eqn{\lambda > 0} the rate, \eqn{r = q - \mu},
#' \eqn{a = |r|}, \eqn{f} the density and \eqn{F} the distribution function.
#'
#' @param distrib A `Laplace2Distrib` object, from [laplace2_distrib()].
#' @param q A numeric vector of quantiles. At \eqn{q = \mu} exactly, `sign(0)`
#'   is 0 and the components are the average of their two one-sided limits.
#' @param theta A named list with components `mu` and `lambda` (positive), each
#'   a numeric vector of length 1 or `n`.
#' @param lower.tail Is the lower tail wanted? A single logical, `TRUE` by
#'   default.
#' @param log Are derivatives of the log probability wanted? A single logical,
#'   `TRUE` by default.
#'
#' @return A named list of three numeric vectors keyed as [hess_names()], each
#'   the length of `q` recycled against `theta`.
#'
#' @seealso [distrib_grad_cdf.Laplace2Distrib()] for the first order;
#'   [distrib_hess_cdf.LaplaceDistrib()] for the scale parametrization;
#'   [laplace2_distrib()].
#'
#' @examples
#' d <- laplace2_distrib()
#' th <- list(mu = 0.3, lambda = 1 / 1.2)
#'
#' # All three components jump across q = mu.
#' distrib_hess_cdf(d, 0.3 + c(-1e-6, 1e-6), th, log = FALSE)
#'
#' # Away from the kink they agree with a central difference of the cdf.
#' q <- c(-1, 2)
#' exact <- distrib_hess_cdf(d, q, th, log = FALSE)
#' fd <- numerical_cdf_deriv(d, q, th, order = 2)
#' max(abs(unlist(exact[names(fd)]) / unlist(fd) - 1))
S7::method(distrib_hess_cdf, Laplace2Distrib) <- function(distrib, q, theta,
                                                          lower.tail = TRUE, log = TRUE) {
  lam <- theta[[2]]
  r <- q - theta[[1]]
  s <- sign(r)
  a <- abs(r)
  f <- distrib_pdf(distrib, q, theta)
  d1 <- stats::setNames(list(-f, r * f / lam), distrib@params)
  d2 <- stats::setNames(list(-lam * s * f, -s * a^2 * f / lam, -(1 / lam - a) * f),
                        hess_names(distrib@params))
  cdf_tail_scale(distrib, distrib_cdf(distrib, q, theta), d1, d2, lower.tail, log)
}


# --- lognormal -------------------------------------------------------------
#
# F(y) = Phi(z) with z = (log y - mu)/sigma, so the family is location-scale on
# the log scale and the same two derivatives appear, expressed through the
# density of Y itself: phi(z) = y sigma f(y), whence
#
#   dF/dmu     = -y f(y)
#   dF/dsigma2 = -y f(y) z / (2 sigma)
#
# the second carrying the extra 1/(2 sigma) because the package parametrizes by
# the variance on the log scale rather than the standard deviation.

#' @title Lognormal Log-CDF Gradient
#' @name distrib_grad_cdf.Lognormal1Distrib
#'
#' @description
#' Closed form. On the log scale the lognormal is a location-scale family, so
#' \eqn{F(q) = \Phi(z)} with \eqn{z = (\log q - \mu)/\sigma}, and the two
#' derivatives are
#' \deqn{\frac{\partial F}{\partial\mu} = -q\,f(q), \qquad
#'       \frac{\partial F}{\partial\sigma^2} = -\frac{q\,f(q)\,z}{2\sigma}.}
#'
#' @details
#' The factor \eqn{q} is the Jacobian of the change of variable: the standard
#' normal density at \eqn{z} is \eqn{\varphi(z) = q\sigma f(q)}, so the
#' location-scale formulas are rewritten in the density of \eqn{Y} itself and
#' nothing on the log scale has to be evaluated. The extra \eqn{1/(2\sigma)} in
#' the second component is the chain rule onto the variance, this
#' parametrization carrying \eqn{\sigma^2} rather than \eqn{\sigma}.
#'
#' Only the gradient is registered; the second derivatives fall to
#' [distrib_hess_cdf.continuous_distrib()], which differences the cdf and
#' reuses this closed gradient for its first-order part.
#'
#' @section Notation:
#' \eqn{\mu} is the mean of \eqn{\log Y}, \eqn{\sigma^2 > 0} its variance,
#' \eqn{z = (\log q - \mu)/\sigma}, \eqn{f} the density of \eqn{Y} and
#' \eqn{\varphi} the standard normal density.
#'
#' @param distrib A `Lognormal1Distrib` object, from [lognormal1_distrib()].
#' @param q A numeric vector of quantiles, positive. At or below zero the
#'   distribution function is zero and its derivatives are too.
#' @param theta A named list with components `mu` (any real value) and `sigma2`
#'   (positive), each a numeric vector of length 1 or `n`.
#' @param lower.tail Is the lower tail wanted? A single logical, `TRUE` by
#'   default.
#' @param log Are derivatives of the log probability wanted? A single logical,
#'   `TRUE` by default.
#'
#' @return A named list of two numeric vectors, `mu` and `sigma2`, each the
#'   length of `q` recycled against `theta`.
#'
#' @seealso [distrib_hess_cdf.continuous_distrib()], the second order;
#'   [distrib_grad_cdf.Gaussian1Distrib()], the family on the log scale;
#'   [lognormal1_distrib()].
#'
#' @examples
#' d <- lognormal1_distrib()
#' th <- list(mu = 0, sigma2 = 1)
#' q <- c(0.5, 1, 3)
#'
#' # The mean component is -q f(q), the Jacobian factor included.
#' all.equal(distrib_grad_cdf(d, q, th, log = FALSE)$mu,
#'           -q * distrib_pdf(d, q, th))
#'
#' # Against a central difference of the cdf, which shares no arithmetic.
#' fd <- numerical_cdf_deriv(d, q, th, order = 1)
#' max(abs(distrib_grad_cdf(d, q, th, log = FALSE)$sigma2 / fd$sigma2 - 1))
S7::method(distrib_grad_cdf, Lognormal1Distrib) <- function(distrib, q, theta,
                                                           lower.tail = TRUE, log = TRUE) {
  s <- sqrt(theta[[2]])
  z <- (base::log(q) - theta[[1]]) / s
  f <- distrib_pdf(distrib, q, theta)
  d1 <- list(mu = -q * f, sigma2 = -q * f * z / (2 * s))
  cdf_tail_scale(distrib, distrib_cdf(distrib, q, theta), d1, NULL, lower.tail, log)
}

# --- inverse Gaussian ------------------------------------------------------
#
# Unusually for a positive family, the distribution function is elementary:
#
#   F(y) = Phi(a) + exp(2/(phi mu)) Phi(b),
#   a = (y/mu - 1)/sqrt(phi y),   b = -(y/mu + 1)/sqrt(phi y),
#
# so it can simply be differentiated. The exponential overflows for small
# phi*mu, so the product exp(2/(phi mu)) Phi(b) is formed on the log scale --
# the factor is huge exactly where Phi(b) is tiny.

# --- discrete families ------------------------------------------------------
#
# For the Poisson the sum defining F telescopes:
#   d/dmu sum_{j<=k} e^-mu mu^j/j!  =  sum_j [-f(j) + f(j-1)]  =  -f(k),
# so the sensitivity of the cdf to the mean is minus the mass at the last point
# retained. The binomial has the companion identity
#   d/dp P(X <= k) = -n dbinom(k, n-1, p),
# and the negative binomial, in the (mu, theta) parametrization,
#   dF(k)/dmu = -f(k) (k + theta)/(theta + mu),
# which tends to the Poisson result as theta grows. The theta direction is a
# derivative of the incomplete beta in its parameter and has no elementary form,
# so it keeps the exact summation.

#' @title Poisson Log-CDF Gradient
#' @name distrib_grad_cdf.PoissonDistrib
#'
#' @description
#' Closed form, and exact: the sum defining \eqn{F} telescopes, so
#' \deqn{\frac{\partial}{\partial\mu}\sum_{j \le k} \frac{e^{-\mu}\mu^j}{j!}
#'       = \sum_{j \le k} \{f(j-1) - f(j)\} = -f(k).}
#' The sensitivity of the distribution function to the mean is minus the mass
#' at the last point retained, so one density evaluation replaces the whole
#' sum.
#'
#' @section Notation:
#' \eqn{\mu > 0} is the mean, \eqn{f} the mass function, \eqn{F} the
#' distribution function and \eqn{k = \lfloor q \rfloor}.
#'
#' @param distrib A `PoissonDistrib` object, from [poisson_distrib()].
#' @param q A numeric vector of quantiles. Non-integer values are floored, as
#'   they are by the distribution function; values below zero give a derivative
#'   of zero.
#' @param theta A named list with one component, `mu` (positive), a numeric
#'   vector of length 1 or `n`.
#' @param lower.tail Is the lower tail wanted? A single logical, `TRUE` by
#'   default.
#' @param log Are derivatives of the log probability wanted? A single logical,
#'   `TRUE` by default.
#'
#' @return A named list with one numeric vector, `mu`, the length of `q`
#'   recycled against `theta`.
#'
#' @seealso [distrib_grad_cdf.discrete_distrib()], the general sum this
#'   replaces; [distrib_grad_cdf.NegBin2Distrib()], whose mean component tends
#'   to this one; [poisson_distrib()].
#'
#' @examples
#' d <- poisson_distrib()
#' q <- c(0, 2, 5)
#'
#' # Minus the mass at the last point retained.
#' all.equal(distrib_grad_cdf(d, q, list(mu = 3), log = FALSE)$mu,
#'           -distrib_pdf(d, q, list(mu = 3)))
#'
#' # On the log scale, and on the upper tail.
#' distrib_grad_cdf(d, q, list(mu = 3))$mu
#' distrib_grad_cdf(d, q, list(mu = 3), lower.tail = FALSE)$mu
S7::method(distrib_grad_cdf, PoissonDistrib) <- function(distrib, q, theta,
                                                         lower.tail = TRUE, log = TRUE) {
  d1 <- list(mu = -distrib_pdf(distrib, floor(q), theta) * (q >= 0))
  cdf_tail_scale(distrib, distrib_cdf(distrib, q, theta), d1, NULL, lower.tail, log)
}

#' @title Binomial Log-CDF Gradient
#' @name distrib_grad_cdf.BinomialDistrib
#'
#' @description
#' Closed form, and exact: the companion identity to the Poisson's,
#' \deqn{\frac{\partial}{\partial p} P(X \le k) = -n\,\mathrm{dbinom}(k, n-1, p).}
#' The whole sum collapses to one binomial mass at one fewer trial, so nothing
#' is summed and nothing is differenced.
#'
#' @section Notation:
#' \eqn{p \in (0,1)} is the success probability, \eqn{n} the number of trials
#' held on the object, and \eqn{k = \lfloor q \rfloor}.
#'
#' @param distrib A `BinomialDistrib` object, from [binomial_distrib()],
#'   carrying the trial count in its `size` property.
#' @param q A numeric vector of quantiles. Non-integer values are floored;
#'   values below zero give a derivative of zero.
#' @param theta A named list with one component, `mu` (the success
#'   probability, strictly between 0 and 1), a numeric vector of length 1 or
#'   `n`.
#' @param lower.tail Is the lower tail wanted? A single logical, `TRUE` by
#'   default.
#' @param log Are derivatives of the log probability wanted? A single logical,
#'   `TRUE` by default.
#'
#' @return A named list with one numeric vector, `mu`, the length of `q`
#'   recycled against `theta`.
#'
#' @seealso [distrib_grad_cdf.BernoulliDistrib()], this identity at
#'   \eqn{n = 1}; [distrib_grad_cdf.discrete_distrib()], the general sum;
#'   [binomial_distrib()].
#'
#' @examples
#' d <- binomial_distrib(size = 10)
#' q <- c(2, 5, 8)
#'
#' # Minus n times a binomial mass at one fewer trial.
#' all.equal(distrib_grad_cdf(d, q, list(mu = 0.3), log = FALSE)$mu,
#'           -10 * dbinom(q, 9, 0.3))
S7::method(distrib_grad_cdf, BinomialDistrib) <- function(distrib, q, theta,
                                                          lower.tail = TRUE, log = TRUE) {
  n <- distrib@size
  d1 <- list(mu = -n * stats::dbinom(floor(q), n - 1, theta[[1]]) * (q >= 0))
  cdf_tail_scale(distrib, distrib_cdf(distrib, q, theta), d1, NULL, lower.tail, log)
}

#' @title Bernoulli Log-CDF Gradient
#' @name distrib_grad_cdf.BernoulliDistrib
#'
#' @description
#' Closed form, the binomial identity at \eqn{n = 1}. The distribution function
#' takes two values, \eqn{F(0) = 1 - p} and \eqn{F(1) = 1}, so the derivative
#' is \eqn{-1} at \eqn{k = 0} and exactly zero at \eqn{k = 1}, the upper value
#' being 1 whatever \eqn{p} is.
#'
#' @section Notation:
#' \eqn{p \in (0,1)} is the success probability and
#' \eqn{k = \lfloor q \rfloor}.
#'
#' @param distrib A `BernoulliDistrib` object, from [bernoulli_distrib()].
#' @param q A numeric vector of quantiles. Non-integer values are floored;
#'   values below zero give a derivative of zero.
#' @param theta A named list with one component, `mu` (the success
#'   probability, strictly between 0 and 1), a numeric vector of length 1 or
#'   `n`.
#' @param lower.tail Is the lower tail wanted? A single logical, `TRUE` by
#'   default.
#' @param log Are derivatives of the log probability wanted? A single logical,
#'   `TRUE` by default. At \eqn{k = 1} the lower-tail probability is 1 and the
#'   log-scale derivative is 0; the upper-tail probability is 0 there, and the
#'   log-scale derivative is `NaN`.
#'
#' @return A named list with one numeric vector, `mu`, the length of `q`
#'   recycled against `theta`.
#'
#' @seealso [distrib_grad_cdf.BinomialDistrib()], of which this is the
#'   \eqn{n = 1} case; [bernoulli_distrib()].
#'
#' @examples
#' d <- bernoulli_distrib()
#'
#' # -1 at zero, and exactly 0 at one.
#' distrib_grad_cdf(d, c(0, 1), list(mu = 0.3), log = FALSE)$mu
#'
#' # On the log scale at k = 0 that is -1 / (1 - p).
#' all.equal(distrib_grad_cdf(d, 0, list(mu = 0.3))$mu, -1 / 0.7)
S7::method(distrib_grad_cdf, BernoulliDistrib) <- function(distrib, q, theta,
                                                           lower.tail = TRUE, log = TRUE) {
  d1 <- list(mu = -stats::dbinom(floor(q), 0, theta[[1]]) * (q >= 0))
  cdf_tail_scale(distrib, distrib_cdf(distrib, q, theta), d1, NULL, lower.tail, log)
}

#' @title Negative Binomial Log-CDF Gradient
#' @name distrib_grad_cdf.NegBin2Distrib
#'
#' @description
#' Closed form in the mean and an exact sum in the dispersion. The mean
#' component collapses like the Poisson's,
#' \deqn{\frac{\partial F(k)}{\partial\mu}
#'       = -f(k)\,\frac{k + \theta}{\theta + \mu},}
#' which tends to the Poisson's \eqn{-f(k)} as \eqn{\theta} grows. The
#' dispersion component is a derivative of the incomplete beta function with
#' respect to its parameter, which has no elementary form, so it keeps the
#' exact summation of [discrete_cdf_deriv()].
#'
#' @section Notation:
#' \eqn{\mu > 0} is the mean, \eqn{\theta > 0} the dispersion, \eqn{f} the mass
#' function, \eqn{F} the distribution function and \eqn{k = \lfloor q \rfloor}.
#'
#' @param distrib A `NegBin2Distrib` object, from [negbin2_distrib()].
#' @param q A numeric vector of quantiles. Non-integer values are floored;
#'   values below zero give a derivative of zero in both components.
#' @param theta A named list with components `mu` (positive) and `theta`
#'   (positive), each a numeric vector of length 1 or `n`. The cost of the
#'   dispersion component grows with the largest quantile, the sum running the
#'   support up to it.
#' @param lower.tail Is the lower tail wanted? A single logical, `TRUE` by
#'   default.
#' @param log Are derivatives of the log probability wanted? A single logical,
#'   `TRUE` by default.
#'
#' @return A named list of two numeric vectors, `mu` and `theta`, each the
#'   length of `q` recycled against the parameters.
#'
#' @seealso [discrete_cdf_deriv()] for the summation the dispersion uses;
#'   [distrib_grad_cdf.PoissonDistrib()] for the limit;
#'   [negbin2_distrib()].
#'
#' @examples
#' d <- negbin2_distrib()
#' q <- c(2, 5)
#' th <- list(mu = 3, theta = 2)
#'
#' # The mean component, written out.
#' all.equal(distrib_grad_cdf(d, q, th, log = FALSE)$mu,
#'           -distrib_pdf(d, q, th) * (q + 2) / (2 + 3))
#'
#' # It tends to the Poisson's -f(k) as the dispersion grows.
#' rbind(negbin = distrib_grad_cdf(d, q, list(mu = 3, theta = 1e8),
#'                                 log = FALSE)$mu,
#'       poisson = distrib_grad_cdf(poisson_distrib(), q, list(mu = 3),
#'                                  log = FALSE)$mu)
S7::method(distrib_grad_cdf, NegBin2Distrib) <- function(distrib, q, theta,
                                                        lower.tail = TRUE, log = TRUE) {
  k <- floor(q)
  f <- distrib_pdf(distrib, k, theta)
  d1 <- list(mu = -f * (k + theta[[2]]) / (theta[[2]] + theta[[1]]) * (q >= 0),
             theta = discrete_cdf_deriv(distrib, q, theta, 1L)[["theta"]])
  cdf_tail_scale(distrib, distrib_cdf(distrib, q, theta), d1, NULL, lower.tail, log)
}

# --- families that are location-scale in only some of their parameters -----
#
# The Student t and the pseudo-Huber are location-scale in (mu, sigma) with a
# further shape parameter nu. The two location-scale derivatives are closed
# form; nu is a derivative of a hypergeometric-type integral and is differenced.

#' CDF Gradient When Only Some Parameters Are Location-Scale
#'
#' @description
#' The [distrib_grad_cdf()] body for a family that is location-scale in
#' \eqn{(\mu, \sigma)} and carries a further shape parameter: the two
#' location-scale directions in closed form, \eqn{-f} and \eqn{-z f}, and the
#' shape directions from [numerical_cdf_deriv()]. The Student t and the
#' pseudo-Huber are the two families that register it.
#'
#' @details
#' # Why the shape is differenced
#'
#' In both families the shape direction is the derivative of a
#' hypergeometric-type integral with respect to its parameter, and has no
#' elementary form. Only those components are passed to the numerical route,
#' through its `which` argument, so the differencing costs cdf evaluations for
#' the shape alone.
#'
#' # What the split is worth
#'
#' It is a speed gain, and the size of it depends on how dear the family's cdf
#' is. On a pseudo-Huber, whose distribution function is itself a quadrature,
#' a gradient at 500 quantiles costs 0.08 s here against 0.18 s when every
#' component is differenced. On accuracy the two agree closely: measured over
#' shapes from 0.1 to 100 and quantiles out to eight scales, the differenced
#' mean component is within \eqn{2\times10^{-11}} to \eqn{3\times10^{-8}}
#' relative of the closed one.
#'
#' @section Notation:
#' \eqn{\mu} is the location, \eqn{\sigma > 0} the scale,
#' \eqn{z = (q-\mu)/\sigma} and \eqn{f} the density.
#'
#' @param distrib An object inheriting from `distrib` whose first two
#'   parameters are a location and a scale, and which has at least one more.
#' @param q A numeric vector of quantiles.
#' @param theta A named list of parameters, the location first and the scale
#'   second.
#' @param lower.tail Is the lower tail wanted? A single logical, `TRUE` by
#'   default.
#' @param log Are derivatives of the log probability wanted? A single logical,
#'   `TRUE` by default.
#'
#' @return A named list of numeric vectors, one per parameter, in the parameter
#'   order.
#'
#' @seealso [loc_scale_grad_cdf()], the body for a family with no shape;
#'   [numerical_cdf_deriv()] for the shape components;
#'   [distrib_grad_cdf.StudentT1Distrib()] and
#'   [distrib_grad_cdf.PseudoHuberDistrib()], the two registrations.
#'
#' @keywords internal
partial_loc_scale_grad_cdf <- function(distrib, q, theta, lower.tail = TRUE, log = TRUE) {
  params <- distrib@params
  z <- (q - theta[[1]]) / theta[[2]]
  f <- distrib_pdf(distrib, q, theta)
  rest <- params[-(1:2)]
  d1 <- c(list(-f, -z * f), numerical_cdf_deriv(distrib, q, theta, 1L, which = rest))
  names(d1) <- params
  cdf_tail_scale(distrib, distrib_cdf(distrib, q, theta), d1, NULL, lower.tail, log)
}

#' @title Student t Log-CDF Gradient
#' @name distrib_grad_cdf.StudentT1Distrib
#'
#' @description
#' Closed form in the location and the scale, \eqn{-f(q)} and \eqn{-z f(q)}
#' with \eqn{z = (q-\mu)/\sigma}; the degrees of freedom are differenced. The
#' method is [partial_loc_scale_grad_cdf()] itself, shared with the
#' pseudo-Huber.
#'
#' @details
#' The derivative of a Student t distribution function with respect to its
#' degrees of freedom has no elementary form, which is the same obstruction the
#' skew t meets in its own \eqn{\nu} components. One central difference of the
#' analytic cdf covers it, and only that component pays for the evaluations.
#'
#' @section Notation:
#' \eqn{\mu} is the location, \eqn{\sigma > 0} the scale, \eqn{\nu > 0} the
#' degrees of freedom, \eqn{z = (q-\mu)/\sigma} and \eqn{f} the density.
#'
#' @param distrib A `StudentT1Distrib` object, from [student_t1_distrib()].
#' @param q A numeric vector of quantiles.
#' @param theta A named list with components `mu`, `sigma` (positive) and `nu`
#'   (positive), each a numeric vector of length 1 or `n`.
#' @param lower.tail Is the lower tail wanted? A single logical, `TRUE` by
#'   default.
#' @param log Are derivatives of the log probability wanted? A single logical,
#'   `TRUE` by default.
#'
#' @return A named list of three numeric vectors, `mu`, `sigma` and `nu`, each
#'   the length of `q` recycled against `theta`.
#'
#' @seealso [partial_loc_scale_grad_cdf()] for the shared body;
#'   [distrib_hess_cdf.StudentT1Distrib()] for the second order;
#'   [student_t1_distrib()].
#'
#' @examples
#' d <- student_t1_distrib()
#' th <- list(mu = 0.3, sigma = 1.2, nu = 6)
#' q <- c(-1, 0.5, 2)
#'
#' # The location component is exact, the density itself.
#' all.equal(distrib_grad_cdf(d, q, th, log = FALSE)$mu,
#'           -distrib_pdf(d, q, th))
#'
#' # The degrees of freedom are differenced; the component is small and negative
#' # in the lower tail, heavier tails putting more mass below a low quantile.
#' distrib_grad_cdf(d, q, th, log = FALSE)$nu
S7::method(distrib_grad_cdf, StudentT1Distrib) <- partial_loc_scale_grad_cdf

#' @title Pseudo-Huber Log-CDF Gradient
#' @name distrib_grad_cdf.PseudoHuberDistrib
#'
#' @description
#' Closed form in the location and the scale, \eqn{-f(q)} and \eqn{-z f(q)}
#' with \eqn{z = (q-\mu)/\sigma}; the shape \eqn{\nu} is differenced. The
#' method is [partial_loc_scale_grad_cdf()] itself, shared with the Student t.
#'
#' @details
#' This family's distribution function is itself a quadrature, so an evaluation
#' of it is dear and the split is worth more here than elsewhere: a gradient at
#' 500 quantiles costs 0.08 s against 0.18 s when all three components are
#' differenced, since the closed pair needs the density alone.
#'
#' @section Notation:
#' \eqn{\mu} is the location, \eqn{\sigma > 0} the scale, \eqn{\nu > 0} the
#' shape, \eqn{z = (q-\mu)/\sigma} and \eqn{f} the density.
#'
#' @param distrib A `PseudoHuberDistrib` object, from [pseudohuber_distrib()].
#' @param q A numeric vector of quantiles.
#' @param theta A named list with components `mu`, `sigma` (positive) and `nu`
#'   (positive), each a numeric vector of length 1 or `n`.
#' @param lower.tail Is the lower tail wanted? A single logical, `TRUE` by
#'   default.
#' @param log Are derivatives of the log probability wanted? A single logical,
#'   `TRUE` by default.
#'
#' @return A named list of three numeric vectors, `mu`, `sigma` and `nu`, each
#'   the length of `q` recycled against `theta`.
#'
#' @seealso [partial_loc_scale_grad_cdf()] for the shared body;
#'   [distrib_hess_cdf.PseudoHuberDistrib()] for the second order;
#'   [pseudohuber_distrib()].
#'
#' @examples
#' d <- pseudohuber_distrib()
#' th <- list(mu = 0.3, sigma = 1.2, nu = 4)
#' q <- c(-1, 0.5, 2)
#'
#' # The location component is exact, the density itself.
#' all.equal(distrib_grad_cdf(d, q, th, log = FALSE)$mu,
#'           -distrib_pdf(d, q, th))
#'
#' # Differencing the quadrature agrees, and costs more.
#' fd <- numerical_cdf_deriv(d, q, th, order = 1)
#' max(abs(fd$mu / distrib_grad_cdf(d, q, th, log = FALSE)$mu - 1))
S7::method(distrib_grad_cdf, PseudoHuberDistrib) <- partial_loc_scale_grad_cdf
