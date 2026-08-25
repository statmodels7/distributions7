#' @include distrib.R generics.R
NULL

#' @title Beta-Binomial Distribution Class, Two Shapes
#' @name BetaBinom2Distrib
#'
#' @description
#' The S7 class of the beta-binomial family in its canonical parametrization,
#' the two beta shapes \eqn{\alpha > 0} and \eqn{\beta > 0}, on the finite
#' support \eqn{\{0, 1, \dots, n\}}. It inherits from `discrete_distrib`, so it
#' answers every generic of the `distrib` contract; the seven methods listed
#' below are registered on it directly and everything else comes from the
#' parent.
#'
#' The class carries an extra property beyond the parent's, `size`: the number
#' of trials \eqn{n}, fixed at construction as it is for
#' [BinomialDistrib()]. Build one with [betabinom2_distrib()], which
#' validates `size`, supplies the two link functions and fills the properties
#' in. This page documents the raw S7 constructor, which validates none of the
#' relationships between them.
#'
#' @inheritParams distrib
#' @param size The number of trials \eqn{n}, a single positive integer. It
#'   belongs to the object, so an object cannot be reused across data sets
#'   whose group sizes differ.
#'
#' @return An S7 object of class `BetaBinom2Distrib`, inheriting from
#'   `discrete_distrib` and from `distrib`. Beyond `size` its properties are
#'   the parent's: `distrib_name`, `dimension`, `bounds`, `params`,
#'   `params_interpretation`, `n_params`, `params_bounds`, `link_params` and
#'   `params_smooth`. For an object built by [betabinom2_distrib()] they hold
#'   `"betabinom2 [size=n]"`, `"univariate"`, `c(0, size)`,
#'   `c("alpha", "beta")`, the interpretations
#'   `c(alpha = "shape", beta = "shape")`, `2`, the domain
#'   \eqn{(0, \infty)} for both, and the two links.
#'
#' @seealso [betabinom2_distrib()] to build one;
#'   [betabinom1_distrib()] for the same law in a mean proportion and a
#'   dispersion; [beta1_distrib()] for the mixing law;
#'   [distrib_pdf.BetaBinom2Distrib()] for the mass function.
#'
#' @section Methods:
#' Registered on this class:
#'   [`distrib_deriv3()`][distrib_deriv3.BetaBinom2Distrib],
#'   [`distrib_deriv4()`][distrib_deriv4.BetaBinom2Distrib],
#'   [`distrib_expected_hessian()`][distrib_expected_hessian.BetaBinom2Distrib],
#'   [`distrib_gradient()`][distrib_gradient.BetaBinom2Distrib],
#'   [`distrib_hessian()`][distrib_hessian.BetaBinom2Distrib],
#'   [`distrib_pdf()`][distrib_pdf.BetaBinom2Distrib],
#'   [`distrib_rng()`][distrib_rng.BetaBinom2Distrib].
#'
#' The four moments [`mean()`][mean.BetaBinom2Distrib],
#' [`variance()`][variance.BetaBinom2Distrib],
#' [`skewness()`][skewness.BetaBinom2Distrib] and
#' [`kurtosis()`][kurtosis.BetaBinom2Distrib] are registered in `moments.R`.
#'
#' The distribution and quantile functions come from
#' [discrete_distrib()] rather than from this class, and on a finite support
#' the parent's cumulative sum of the mass is exact.
#'
#' @examples
#' d <- betabinom2_distrib(size = 10)
#' S7::S7_inherits(d, discrete_distrib)
#'
#' # The trial count is a property of the object, not an entry of theta.
#' d@size
#' d@bounds
#'
#' # Both shapes are positive, so both ride a log by default.
#' d@params
#' vapply(d@link_params, function(l) l@link_name, character(1))
BetaBinom2Distrib <- S7::new_class("BetaBinom2Distrib",
  parent = discrete_distrib,
  properties = list(size = S7::class_numeric)
)

#' One Component of a Beta-Binomial Derivative in the Shapes
#'
#' @description
#' Returns the derivative component of a given order and multi-index for the
#' beta-binomial log-mass in its two shapes. The log-mass is a sum of log-gamma
#' terms, so its derivative of order \eqn{k+1} is the same sum with
#' \eqn{\psi^{(k)}} in place of \eqn{\log\Gamma}. All four orders come from
#' this one routine.
#'
#' @details
#' With \eqn{a} and \eqn{b} the shapes and \eqn{n} the size, the log-mass is
#' \deqn{\log\Gamma(y+a) + \log\Gamma(n-y+b) - \log\Gamma(n+a+b)
#'       - \log\Gamma(a) - \log\Gamma(b) + \log\Gamma(a+b)}
#' up to a term free of the parameters. A derivative in \eqn{a} alone
#' differentiates the first, third, fourth and sixth terms; one in \eqn{b}
#' alone the second, third, fifth and sixth; a **mixed** one only the two
#' terms carrying \eqn{a+b}, which is why a mixed component of any order is
#' \eqn{-\psi^{(k)}(n+a+b) + \psi^{(k)}(a+b)} and free of the data.
#'
#' @param y A numeric vector of counts.
#' @param a,b The two shapes, each a numeric vector of length 1 or of the
#'   length of `y`, strictly positive.
#' @param n The size, a single positive integer.
#' @param k The polygamma order, one less than the derivative order: 0 for the
#'   score, 1 for the Hessian, and so on.
#' @param i The number of \eqn{a} indices in the component.
#' @param j The number of \eqn{b} indices. `i + j` is the derivative order.
#'
#' @return A numeric vector of the recycled length of the inputs. When both
#'   `i` and `j` are positive the value is constant along it, a mixed
#'   component not depending on `y`.
#'
#' @seealso [betabinom2_derivs()], which calls this for every multi-index of
#'   an order, and [betabinom2_distrib()] for the family.
#'
#' @keywords internal
betabinom2_component <- function(y, a, b, n, k, i, j) {
  mixed <- -psigamma(n + a + b, deriv = k) + psigamma(a + b, deriv = k)
  if (i > 0 && j > 0) return(rep(mixed, length.out = length(y)))
  if (j == 0) {
    return(psigamma(y + a, deriv = k) - psigamma(a, deriv = k) + mixed)
  }
  psigamma(n - y + b, deriv = k) - psigamma(b, deriv = k) + mixed
}

#' Every Component of a Beta-Binomial Derivative in the Shapes
#'
#' @description
#' Assembles the full set of derivative components of a given order by calling
#' [betabinom2_component()] once per distinct multi-index, and names them as
#' [deriv_names()] names them.
#'
#' @param y A numeric vector of counts.
#' @param a,b The two shapes, each a numeric vector of length 1 or of the
#'   length of `y`, strictly positive.
#' @param n The size, a single positive integer.
#' @param order The derivative order, an integer from 1 to 4.
#' @param params The parameter names, `c("alpha", "beta")` for this family,
#'   used to build the component names and the multi-indices in the same
#'   order.
#'
#' @return A named list of component vectors, one per distinct multi-index of
#'   the given order: two at order 1, three at order 2, four at order 3 and
#'   five at order 4. Each has the recycled length of the inputs.
#'
#' @seealso [betabinom2_component()] for one component,
#'   [deriv_names()] and [deriv_indices()] for the enumeration, and
#'   [betabinom2_distrib()] for the family.
#'
#' @keywords internal
betabinom2_derivs <- function(y, a, b, n, order, params) {
  idx <- deriv_indices(params, order)
  nm <- deriv_names(params, order)
  stats::setNames(lapply(seq_along(nm), function(m) {
    id <- idx[[m]]
    betabinom2_component(y, a, b, n, order - 1L, sum(id == 1L), sum(id == 2L))
  }), nm)
}

# --- S7 METHODS IMPLEMENTATION ---

#' @title Beta-Binomial Mass Function in Its Shapes
#' @name distrib_pdf.BetaBinom2Distrib
#' @description
#' Computes the beta-binomial mass in the two shapes,
#' \deqn{P(Y = y) = \binom{n}{y}\dfrac{B(y+\alpha,\; n-y+\beta)}{B(\alpha,\beta)},}
#' with \eqn{B} the beta function and \eqn{n} the object's `size`. The support
#' is tested first, so a count outside \eqn{\{0, \dots, n\}} or not an integer
#' returns a mass of 0 rather than reaching [base::lchoose()], which would warn
#' on a non-integer argument.
#'
#' The evaluation route is chosen by the concentration; see
#' [betabinom_log_mass()], which this calls.
#'
#' @param distrib A `BetaBinom2Distrib` object, from [betabinom2_distrib()].
#' @param y A numeric vector of counts. A value outside \eqn{\{0, \dots, n\}}
#'   or not an integer gives a mass of 0, or `-Inf` with `log = TRUE`.
#' @param theta A named list with components `alpha` and `beta`, each a numeric
#'   vector of length 1 or of the length of `y`. A component of length 1 is
#'   recycled. Both must be strictly positive.
#' @param log Logical of length 1. When `TRUE` the log-mass is returned.
#'   Defaults to `FALSE`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of probabilities in \eqn{[0, 1]}, of length
#'   `length(y)`, one value per observation.
#'
#' @seealso [betabinom_log_mass()] for the two evaluation routes,
#'   [distrib_gradient.BetaBinom2Distrib()] for the derivatives of the
#'   log-mass, [distrib_pdf.BetaBinom1Distrib()] for the same quantity in the
#'   mean and dispersion, and [distrib_pdf()] for the generic.
#'
#' @examples
#' d <- betabinom2_distrib(size = 10)
#' th <- list(alpha = 2, beta = 3)
#'
#' # The mass over the whole support sums to one.
#' m <- distrib_pdf(d, 0:10, th)
#' round(m, 4)
#' sum(m)
#'
#' # It is the beta-binomial mass written out.
#' all.equal(m, choose(10, 0:10) * beta(0:10 + 2, 10 - 0:10 + 3) / beta(2, 3))
#'
#' # Off the support, and at a non-integer count, the mass is zero.
#' distrib_pdf(d, c(-1, 2.5, 11), th)
#'
#' # The same law as betabinom1 at mu = alpha / (alpha + beta) and
#' # sigma = 1 / (alpha + beta).
#' all.equal(m, distrib_pdf(betabinom1_distrib(size = 10), 0:10,
#'                          list(mu = 0.4, sigma = 0.2)))
S7::method(distrib_pdf, BetaBinom2Distrib) <- function(distrib, y, theta, log = FALSE, ...) {
  a <- theta[[1]]
  b <- theta[[2]]
  n <- distrib@size
  # The support is tested BEFORE lchoose(), which warns on a non-integer count
  # rather than returning nothing: a mass function evaluated off its support is
  # zero, not a numerical complaint.
  ok <- y >= 0 & y <= n & y == round(y)
  out <- rep(-Inf, length(y))
  if (any(ok)) {
    out[ok] <- betabinom_log_mass(y[ok], a, b, n)
  }
  if (log) out else exp(out)
}

#' Log-Mass of the Beta-Binomial
#'
#' @description
#' Returns
#' \eqn{\log\binom{n}{y} + \log B(y+\alpha, n-y+\beta) - \log B(\alpha, \beta)}
#' by whichever of two routes is accurate at the shapes given. The choice is
#' made from the shapes alone, so it is deterministic and costs one comparison.
#'
#' @details
#' # Why there are two routes
#'
#' The two beta functions are each of magnitude
#' \eqn{(\alpha+\beta)\log(\alpha+\beta)} while their difference is of order
#' one, so forming the mass as that difference carries an absolute error of
#' \eqn{\varepsilon} times the larger magnitude. That route is used only while
#' the error stays below `1e-8`.
#'
#' Beyond it the shifts \eqn{y}, \eqn{n-y} and \eqn{n} are **integers**, so
#' each log-gamma difference is an exact sum of logarithms,
#' \deqn{\log\Gamma(\alpha+y) - \log\Gamma(\alpha) =
#'       \sum_{j=0}^{y-1}\log(\alpha+j),}
#' and the mass follows from three such sums, none of which forms a quantity
#' larger than \eqn{n\log(\alpha+\beta)}. The sums also give the binomial limit
#' correctly as the shapes grow at a fixed ratio: measured at \eqn{n = 10},
#' \eqn{y = 3} and \eqn{\alpha/(\alpha+\beta) = 0.4}, the log-mass agrees with
#' the binomial one to twelve figures at a concentration of \eqn{10^{14}},
#' where the beta-function route is wrong in the third decimal.
#'
#' The cost is \eqn{O(n)} in the size, which is why the route is taken only
#' where it is needed.
#'
#' @param y A numeric vector of counts, already known to lie on the support.
#'   The caller tests that; this function does not.
#' @param a,b The two shapes, each a numeric vector of length 1 or of the
#'   length of `y`, strictly positive.
#' @param n The size, a single positive integer.
#'
#' @return A numeric vector of log-probabilities, of the recycled length of the
#'   inputs.
#'
#' @seealso [distrib_pdf.BetaBinom2Distrib()], which calls this after testing
#'   the support, and [betabinom2_distrib()] for the family.
#' @keywords internal
betabinom_log_mass <- function(y, a, b, n) {
  s <- a + b
  # the two large terms are grouped so that their difference is taken before
  # the term of order one is added to it
  if (all(is.finite(s)) &&
      max(lgamma(s + n)) * .Machine$double.eps < 1e-8) {
    return(lchoose(n, y) + (lbeta(y + a, n - y + b) - lbeta(a, b)))
  }
  m <- length(y)
  a <- rep_len(a, m)
  b <- rep_len(b, m)
  s1 <- numeric(m)
  s2 <- numeric(m)
  s3 <- numeric(m)
  for (j in seq_len(n) - 1) {
    s1 <- s1 + (j < y) * log(a + j)
    s2 <- s2 + (j < n - y) * log(b + j)
    s3 <- s3 + log(a + b + j)
  }
  lchoose(n, y) + s1 + s2 - s3
}

#' @title Beta-Binomial Random Generation in Its Shapes
#' @name distrib_rng.BetaBinom2Distrib
#' @description
#' Draws `n` independent beta-binomial counts by the two-stage hierarchy the
#' family is defined by: a success probability from
#' \eqn{\mathrm{Beta}(\alpha, \beta)}, then a count from
#' \eqn{\mathrm{Binomial}(n_{\mathrm{trials}}, p)} at that probability, one
#' fresh probability per draw. The draws depend on `.Random.seed` in the usual
#' way and consume two of R's streams per variate.
#'
#' @param distrib A `BetaBinom2Distrib` object, from [betabinom2_distrib()].
#' @param n A single positive integer, the number of draws. Note that the
#'   number of **trials** is the object's `size` property, not this argument.
#' @param theta A named list with components `alpha` and `beta`, each a numeric
#'   vector of length 1 or of length `n`. A component of length 1 is recycled.
#'   Both must be strictly positive.
#'
#' @return A numeric vector of `n` counts in \eqn{\{0, \dots, size\}}.
#'
#' @seealso [distrib_rng.BetaBinom1Distrib()] for the same draw in the mean and
#'   dispersion, [fit_distrib()] to estimate the shapes back from a sample, and
#'   [distrib_rng()] for the generic.
#'
#' @examples
#' d <- betabinom2_distrib(size = 10)
#' th <- list(alpha = 2, beta = 3)
#'
#' # The sample moments recover the closed forms.
#' set.seed(2)
#' z <- distrib_rng(d, 2e5, th)
#' rbind(sample = c(mean = mean(z), var = var(z)),
#'       theoretical = c(mean(d, th), variance(d, th)))
#'
#' # The counts are bounded by the trial count, which is the object's size and
#' # not the argument n.
#' range(z)
S7::method(distrib_rng, BetaBinom2Distrib) <- function(distrib, n, theta) {
  p <- stats::rbeta(n, shape1 = theta[[1]], shape2 = theta[[2]])
  stats::rbinom(n, size = distrib@size, prob = p)
}

#' @title Beta-Binomial Score in Its Shapes
#' @name distrib_gradient.BetaBinom2Distrib
#' @description
#' Computes the first derivatives of the beta-binomial log-mass with respect to
#' the two shapes, one value per observation, in closed form:
#' \deqn{\dfrac{\partial\ell}{\partial\alpha}
#'         = \psi(y+\alpha) - \psi(\alpha)
#'           - \psi(n+\alpha+\beta) + \psi(\alpha+\beta),}
#' and the same with \eqn{n-y} and \eqn{\beta} in place of \eqn{y} and
#' \eqn{\alpha}. The two share the term in \eqn{\alpha+\beta}, which is the
#' only part a mixed second derivative keeps.
#'
#' With `scale = "link"` the generic applies the chain rule for the links the
#' family carries before returning. This method always returns the parameter
#' scale; the transformation happens in the generic.
#'
#' @param distrib A `BetaBinom2Distrib` object, from [betabinom2_distrib()].
#' @param y A numeric vector of counts in \eqn{\{0, \dots, n\}}.
#' @param theta A named list with components `alpha` and `beta`, each a numeric
#'   vector of length 1 or of the length of `y`. A component of length 1 is
#'   recycled. Both must be strictly positive.
#' @param scale One of `"parameter"` (the default) or `"link"`, matched by
#'   [base::match.arg()]. Read by the generic, not by this method.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A named list of two numeric vectors, `alpha` and `beta`, each of
#'   length `max(length(y), length(alpha), length(beta))`.
#'
#' @section Notation:
#' \eqn{\ell} is the log-mass of one observation, \eqn{\alpha, \beta > 0} the
#' two beta shapes, \eqn{n} the trial count and \eqn{\psi} the digamma
#' function, `digamma()` in R.
#'
#' @section Large concentration:
#' Each digamma difference cancels to leading order as the shapes grow at a
#' fixed ratio, and this method writes it out directly. Measured at
#' \eqn{n = 10}, \eqn{y = 3} and a mean proportion of 0.4, the relative error
#' against the limiting \eqn{-5/(2S)} with \eqn{S = \alpha+\beta} is
#' \eqn{6\times10^{-8}} at \eqn{S = 10^8}, \eqn{4\times10^{-4}} at
#' \eqn{10^{12}} and 1.8 at \eqn{10^{15}}. Where a concentration of that size
#' is reachable, use [distrib_gradient.BetaBinom1Distrib()], whose compiled
#' kernel forms the same difference as a sum of reciprocals and holds nine
#' figures there.
#'
#' @seealso [distrib_hessian.BetaBinom2Distrib()] for the second derivatives,
#'   [distrib_expected_hessian.BetaBinom2Distrib()] for their expectation,
#'   [distrib_gradient.BetaBinom1Distrib()] for the same quantity in the mean
#'   and dispersion, and [distrib_gradient()] for the generic.
#'
#' @examples
#' d <- betabinom2_distrib(size = 10)
#' th <- list(alpha = 2, beta = 3)
#' g <- distrib_gradient(d, 0:10, th)
#'
#' # The alpha component, written out.
#' ds <- digamma(5) - digamma(15)
#' all.equal(g$alpha, digamma(0:10 + 2) - digamma(2) + ds)
#'
#' # It is the derivative of the log-mass, so a central difference reproduces
#' # it.
#' eps <- 1e-6
#' all.equal((distrib_pdf(d, 0:10, list(alpha = 2 + eps, beta = 3), log = TRUE) -
#'            distrib_pdf(d, 0:10, list(alpha = 2 - eps, beta = 3), log = TRUE)) /
#'             (2 * eps), g$alpha, tolerance = 1e-6)
#'
#' # The score has mean zero over the support: the first Bartlett identity.
#' w <- distrib_pdf(d, 0:10, th)
#' vapply(g, function(v) sum(w * v), numeric(1))
S7::method(distrib_gradient, BetaBinom2Distrib) <- function(distrib, y, theta,
                                                             scale = c("parameter", "link"), ...) {
  a <- theta[[1]]
  b <- theta[[2]]
  n <- distrib@size
  ds <- digamma(a + b) - digamma(n + a + b)
  list(alpha = digamma(y + a) - digamma(a) + ds,
       beta = digamma(n - y + b) - digamma(b) + ds)
}

#' @title Beta-Binomial Observed Hessian in Its Shapes
#' @name distrib_hessian.BetaBinom2Distrib
#' @description
#' Computes the three distinct second derivatives of the beta-binomial log-mass
#' with respect to the two shapes, one value per observation, in closed form.
#' They are the score's digamma differences with the trigamma function in its
#' place:
#' \deqn{\dfrac{\partial^2\ell}{\partial\alpha^2}
#'         = \psi_1(y+\alpha) - \psi_1(\alpha)
#'           - \psi_1(n+S) + \psi_1(S), \qquad S = \alpha+\beta,}
#' and the same with \eqn{n-y} and \eqn{\beta}. The **mixed** component keeps
#' only the shared term,
#' \eqn{-\psi_1(n+S) + \psi_1(S)}, the two shapes entering the log-mass
#' separately otherwise, so it does not depend on the data at all and equals
#' its own expectation at every observation.
#'
#' @param distrib A `BetaBinom2Distrib` object, from [betabinom2_distrib()].
#' @param y A numeric vector of counts in \eqn{\{0, \dots, n\}}.
#' @param theta A named list with components `alpha` and `beta`, each a numeric
#'   vector of length 1 or of the length of `y`. A component of length 1 is
#'   recycled. Both must be strictly positive.
#' @param scale One of `"parameter"` (the default) or `"link"`, matched by
#'   [base::match.arg()]. Read by the generic, not by this method.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A named list of three numeric vectors, `alpha_alpha`, `beta_beta`
#'   and `alpha_beta`, each of length
#'   `max(length(y), length(alpha), length(beta))`. The three name the distinct
#'   entries of a symmetric \eqn{2 \times 2} matrix per observation.
#'
#' @section Notation:
#' \eqn{\ell} is the log-mass of one observation, \eqn{\alpha, \beta > 0} the
#' two beta shapes, \eqn{n} the trial count and \eqn{\psi_1} the trigamma
#' function, `trigamma()` in R.
#'
#' @seealso [distrib_gradient.BetaBinom2Distrib()] for the score,
#'   [distrib_expected_hessian.BetaBinom2Distrib()] for the expectation of this
#'   quantity, [distrib_deriv3.BetaBinom2Distrib()] for the order above, and
#'   [distrib_hessian()] for the generic.
#'
#' @examples
#' d <- betabinom2_distrib(size = 10)
#' th <- list(alpha = 2, beta = 3)
#' h <- distrib_hessian(d, 0:10, th)
#' names(h)
#'
#' # The mixed component is free of the data: one number, repeated.
#' unique(h$alpha_beta)
#' -trigamma(10 + 5) + trigamma(5)
#'
#' # A central difference of the score reproduces the pure-alpha component.
#' eps <- 1e-5
#' up <- distrib_gradient(d, 0:10, list(alpha = 2 + eps, beta = 3))$alpha
#' dn <- distrib_gradient(d, 0:10, list(alpha = 2 - eps, beta = 3))$alpha
#' all.equal((up - dn) / (2 * eps), h$alpha_alpha, tolerance = 1e-6)
S7::method(distrib_hessian, BetaBinom2Distrib) <- function(distrib, y, theta,
                                                            scale = c("parameter", "link"), ...) {
  d <- betabinom2_derivs(y, theta[[1]], theta[[2]], distrib@size, 2L,
                         distrib@params)
  d[hess_names(distrib@params)]
}

#' @title Beta-Binomial Expected Hessian in Its Shapes
#' @name distrib_expected_hessian.BetaBinom2Distrib
#' @description
#' Returns the expectation of the observed Hessian under the model,
#' \eqn{\sum_{k=0}^{n} P(Y=k)\,\partial^2\ell/\partial\theta_i\partial\theta_j}
#' evaluated at \eqn{y = k}. The family being discrete on
#' \eqn{\{0, \dots, n\}}, that is an **exact finite sum** of at most \eqn{n+1}
#' terms, so the answer is the expectation to machine precision. `approx` and
#' `nsim` are therefore ignored: every strategy returns the same three
#' numbers.
#'
#' The mixed entry does not vanish, so the two shapes are not orthogonal and
#' their estimates are asymptotically correlated.
#'
#' @param distrib A `BetaBinom2Distrib` object, from [betabinom2_distrib()].
#' @param y A numeric vector of counts. Only its length is read, the
#'   expectation not depending on the data; the values are ignored.
#' @param theta A named list with components `alpha` and `beta`, each a numeric
#'   vector of length 1. Both must be strictly positive. One weighted sum is
#'   built for the whole call, so a parameter varying by observation is not
#'   supported here.
#' @param scale One of `"parameter"` (the default) or `"link"`, matched by
#'   [base::match.arg()]. Read by the generic, not by this method.
#' @param approx Ignored here, the expectation being an exact sum. Accepted so
#'   that the signature matches the generic's, where it selects between
#'   `"bartlett"`, `"integrate"`, `"mc"` and `"opg"`.
#' @param nsim Ignored here, for the same reason. Defaults to `10000`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A named list of three numeric vectors, `alpha_alpha`, `beta_beta`
#'   and `alpha_beta`, each of length `length(y)` and each constant along it.
#'
#' @seealso [betabinom2_expected()] for the summation,
#'   [distrib_hessian.BetaBinom2Distrib()] for the quantity this is the
#'   expectation of, and [distrib_expected_hessian()] for the generic.
#'
#' @examples
#' d <- betabinom2_distrib(size = 10)
#' th <- list(alpha = 2, beta = 3)
#' eh <- distrib_expected_hessian(d, 0:10, th)
#' vapply(eh, function(v) v[1], numeric(1))
#'
#' # It is the mass-weighted sum of the observed Hessian over the support,
#' # written out here by hand and agreeing exactly.
#' w <- distrib_pdf(d, 0:10, th)
#' vapply(distrib_hessian(d, 0:10, th), function(v) sum(w * v), numeric(1))
#'
#' # The strategy argument is inert, the expectation being an exact sum.
#' identical(eh, distrib_expected_hessian(d, 0:10, th, approx = "mc",
#'                                        nsim = 50))
S7::method(distrib_expected_hessian, BetaBinom2Distrib) <- function(distrib, y, theta,
                                                                     scale = c("parameter", "link"),
                                                                     approx = c("bartlett", "integrate", "mc", "opg"),
                                                                     nsim = 10000, ...) {
  betabinom2_expected(distrib, y, theta, 2L)
}

#' Expected Derivatives of the Beta-Binomial by Exact Summation
#'
#' @description
#' Averages every component of a derivative of the given order over the support
#' \eqn{\{0, \dots, n\}}, weighted by the mass function. The support being
#' finite, the average is the expectation exactly, with no quadrature error and
#' no sampling error.
#'
#' @param distrib A [BetaBinom2Distrib()] object, read for its `size` and
#'   `params`.
#' @param y A numeric vector, used only for its length: the result is one
#'   constant per component, repeated to that length.
#' @param theta A named list with components `alpha` and `beta`, each a numeric
#'   vector of length 1, strictly positive.
#' @param order The derivative order, an integer from 2 to 4.
#'
#' @return A named list of component vectors, each of length `length(y)` and
#'   constant along it. At order 2 the names and their order are
#'   [hess_names()]'s; above it they are [deriv_names()]'s.
#'
#' @seealso [distrib_expected_hessian.BetaBinom2Distrib()],
#'   [distrib_deriv3.BetaBinom2Distrib()] and
#'   [distrib_deriv4.BetaBinom2Distrib()], which call this;
#'   [betabinom2_derivs()] for the components being averaged.
#'
#' @keywords internal
betabinom2_expected <- function(distrib, y, theta, order) {
  a <- theta[[1]]
  b <- theta[[2]]
  n <- distrib@size
  supp <- 0:n
  w <- distrib_pdf(distrib, supp, theta)
  d <- betabinom2_derivs(supp, a, b, n, order, distrib@params)
  nm <- if (order == 2L) hess_names(distrib@params) else names(d)
  stats::setNames(lapply(nm, function(k) {
    rep(sum(w * d[[k]]), length.out = length(y))
  }), nm)
}

#' @title Beta-Binomial Third-Order Derivatives in Its Shapes
#' @name distrib_deriv3.BetaBinom2Distrib
#' @description
#' Computes the four distinct third derivatives of the beta-binomial log-mass
#' in the two shapes, in closed form. The log-mass is a sum of log-gamma terms,
#' so the third derivative is the same sum with \eqn{\psi_2}, the second
#' polygamma, in place of \eqn{\log\Gamma}. A component naming both shapes
#' keeps only the terms in \eqn{\alpha+\beta} and is free of the data.
#'
#' With `expected = TRUE` the expectation is an **exact finite sum** over the
#' support, so `approx` and `nsim` are ignored on both branches.
#'
#' @param distrib A `BetaBinom2Distrib` object, from [betabinom2_distrib()].
#' @param y A numeric vector of counts in \eqn{\{0, \dots, n\}}. With
#'   `expected = TRUE` only its length is read.
#' @param theta A named list with components `alpha` and `beta`, each a numeric
#'   vector of length 1 or of the length of `y`. A component of length 1 is
#'   recycled. Both must be strictly positive.
#' @param expected Logical of length 1. When `TRUE` the exact expectation under
#'   the model is returned in place of the value at the data. Defaults to
#'   `FALSE`.
#' @param scale One of `"parameter"` (the default) or `"link"`, matched by
#'   [base::match.arg()]. Read by the generic, not by this method.
#' @param approx Ignored here, both branches being exact. Accepted so that the
#'   signature matches the generic's.
#' @param nsim Ignored here, for the same reason. Defaults to `10000`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A named list of four numeric vectors, `alpha_alpha_alpha`,
#'   `alpha_alpha_beta`, `alpha_beta_beta` and `beta_beta_beta`, each of length
#'   `max(length(y), length(alpha), length(beta))`.
#'
#' @seealso [distrib_hessian.BetaBinom2Distrib()] for the order below,
#'   [distrib_deriv4.BetaBinom2Distrib()] for the order above,
#'   [betabinom2_derivs()] for the assembly, and [distrib_deriv3()] for the
#'   generic.
#'
#' @examples
#' d <- betabinom2_distrib(size = 10)
#' th <- list(alpha = 2, beta = 3)
#' d3 <- distrib_deriv3(d, 0:10, th)
#' names(d3)
#'
#' # A central difference of the Hessian reproduces the pure-alpha component.
#' eps <- 1e-5
#' up <- distrib_hessian(d, 0:10, list(alpha = 2 + eps, beta = 3))$alpha_alpha
#' dn <- distrib_hessian(d, 0:10, list(alpha = 2 - eps, beta = 3))$alpha_alpha
#' all.equal((up - dn) / (2 * eps), d3$alpha_alpha_alpha, tolerance = 1e-6)
#'
#' # The expected branch is the mass-weighted sum, exactly.
#' w <- distrib_pdf(d, 0:10, th)
#' rbind(expected = vapply(distrib_deriv3(d, 0, th, expected = TRUE),
#'                         function(v) v[1], numeric(1)),
#'       summed = vapply(d3, function(v) sum(w * v), numeric(1)))
S7::method(distrib_deriv3, BetaBinom2Distrib) <- function(distrib, y, theta, expected = FALSE,
                                                           scale = c("parameter", "link"),
                                                           approx = c("integrate", "bartlett", "mc", "opg"),
                                                           nsim = 10000, ...) {
  if (expected) {
    betabinom2_expected(distrib, y, theta, 3L)
  } else {
    betabinom2_derivs(y, theta[[1]], theta[[2]], distrib@size, 3L, distrib@params)
  }
}

#' @title Beta-Binomial Fourth-Order Derivatives in Its Shapes
#' @name distrib_deriv4.BetaBinom2Distrib
#' @description
#' Computes the five distinct fourth derivatives of the beta-binomial log-mass
#' in the two shapes, in closed form, by the construction
#' [distrib_deriv3.BetaBinom2Distrib()] describes carried one order further:
#' the same sum of terms with \eqn{\psi_3}, the third polygamma, in place of
#' \eqn{\log\Gamma}.
#'
#' With `expected = TRUE` the expectation is an **exact finite sum** over the
#' support, so `approx` and `nsim` are ignored on both branches.
#'
#' @param distrib A `BetaBinom2Distrib` object, from [betabinom2_distrib()].
#' @param y A numeric vector of counts in \eqn{\{0, \dots, n\}}. With
#'   `expected = TRUE` only its length is read.
#' @param theta A named list with components `alpha` and `beta`, each a numeric
#'   vector of length 1 or of the length of `y`. A component of length 1 is
#'   recycled. Both must be strictly positive.
#' @param expected Logical of length 1. When `TRUE` the exact expectation under
#'   the model is returned in place of the value at the data. Defaults to
#'   `FALSE`.
#' @param scale One of `"parameter"` (the default) or `"link"`, matched by
#'   [base::match.arg()]. Read by the generic, not by this method.
#' @param approx Ignored here, both branches being exact. Accepted so that the
#'   signature matches the generic's.
#' @param nsim Ignored here, for the same reason. Defaults to `10000`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A named list of five numeric vectors, `alpha_alpha_alpha_alpha`,
#'   `alpha_alpha_alpha_beta`, `alpha_alpha_beta_beta`,
#'   `alpha_beta_beta_beta` and `beta_beta_beta_beta`, each of length
#'   `max(length(y), length(alpha), length(beta))`.
#'
#' @seealso [distrib_deriv3.BetaBinom2Distrib()] for the order below,
#'   [betabinom2_derivs()] for the assembly, and [distrib_deriv4()] for the
#'   generic.
#'
#' @examples
#' d <- betabinom2_distrib(size = 10)
#' th <- list(alpha = 2, beta = 3)
#' d4 <- distrib_deriv4(d, 0:10, th)
#' names(d4)
#'
#' # A central difference of the third order reproduces a mixed component.
#' eps <- 1e-5
#' up <- distrib_deriv3(d, 0:10, list(alpha = 2, beta = 3 + eps))$alpha_alpha_beta
#' dn <- distrib_deriv3(d, 0:10, list(alpha = 2, beta = 3 - eps))$alpha_alpha_beta
#' all.equal((up - dn) / (2 * eps), d4$alpha_alpha_beta_beta, tolerance = 1e-5)
#'
#' # The expected branch is the mass-weighted sum, exactly.
#' w <- distrib_pdf(d, 0:10, th)
#' rbind(expected = vapply(distrib_deriv4(d, 0, th, expected = TRUE),
#'                         function(v) v[1], numeric(1)),
#'       summed = vapply(d4, function(v) sum(w * v), numeric(1)))
S7::method(distrib_deriv4, BetaBinom2Distrib) <- function(distrib, y, theta, expected = FALSE,
                                                           scale = c("parameter", "link"),
                                                           approx = c("integrate", "bartlett", "mc", "opg"),
                                                           nsim = 10000, ...) {
  if (expected) {
    betabinom2_expected(distrib, y, theta, 4L)
  } else {
    betabinom2_derivs(y, theta[[1]], theta[[2]], distrib@size, 4L, distrib@params)
  }
}


#' Beta-Binomial Distribution, Two Shapes
#'
#' @description
#' Builds the distribution object for the beta-binomial family in its canonical
#' parametrization, the two beta shapes \eqn{\alpha > 0} and \eqn{\beta > 0},
#' on the finite support \eqn{\{0, 1, \dots, n\}}. The returned object carries
#' closed-form derivatives of the log-mass to fourth order and expectations
#' that are exact finite sums over the support, so every generic of the toolkit
#' answers without a quadrature.
#'
#' This is the same law as [betabinom1_distrib()], which is written in a mean
#' proportion and a dispersion. Use this one when the shapes themselves are the
#' quantities of interest, and that one when a regression on the mean is.
#'
#' @param size The number of trials \eqn{n}, a single positive integer. It is a
#'   constant of the distribution and not a parameter, as for
#'   [binomial_distrib()], so an object cannot be reused across data sets whose
#'   group sizes differ. Anything else signals an error naming the argument.
#' @param link_alpha A `link` object from `linkfunctions7` for the shape
#'   \eqn{\alpha}. Defaults to [linkfunctions7::log_link()], which maps
#'   \eqn{(0, \infty)} onto the line and so keeps every fitted value positive.
#' @param link_beta A `link` object from `linkfunctions7` for the shape
#'   \eqn{\beta}. Defaults to [linkfunctions7::log_link()], for the same
#'   reason.
#'
#' @details
#' # The parametrization
#'
#' The mass on \eqn{y \in \{0, \dots, n\}} is
#' \deqn{P(Y=y) = \binom{n}{y}\frac{B(y+\alpha,\; n-y+\beta)}{B(\alpha, \beta)},}
#' with \eqn{B} the beta function, and the family is the binomial with its
#' success probability drawn from \eqn{\mathrm{Beta}(\alpha, \beta)}. Writing
#' \eqn{S = \alpha + \beta} for the concentration, the moments are
#' \deqn{\mathbb{E}[Y] = \frac{n\alpha}{S}, \qquad
#'       \operatorname{Var}(Y) = \frac{n\alpha\beta\,(S+n)}{S^{2}(S+1)}.}
#' The correspondence with [betabinom1_distrib()] is \eqn{\mu = \alpha/S} and
#' \eqn{\sigma = 1/S}, so a large concentration is a small dispersion and the
#' binomial limit.
#'
#' # Derivatives
#'
#' The log-mass is a sum of log-gamma terms, so a derivative of order \eqn{k}
#' replaces each by \eqn{\psi^{(k-1)}} and all four orders come from one
#' routine, [betabinom2_derivs()]. A component naming both shapes keeps only
#' the terms in \eqn{S} and is therefore free of the data at every order.
#'
#' Every expectation is an **exact finite sum** over \eqn{\{0, \dots, n\}}:
#' the support is bounded, so there is nothing to integrate over.
#'
#' # The cancellation at a large concentration
#'
#' The two beta functions of the mass are each of magnitude \eqn{S\log S} while
#' their difference is of order one, so writing the mass as that difference
#' loses one digit per factor of ten in \eqn{S}. [betabinom_log_mass()]
#' switches to a sum of logarithms past a measured threshold and stays exact:
#' at \eqn{n = 10}, \eqn{y = 3} and \eqn{\alpha/S = 0.4}, the log-mass agrees
#' with the binomial one to twelve figures at \eqn{S = 10^{14}}, where the
#' direct route is wrong in the third decimal.
#'
#' **The derivatives are not rewritten that way here**, and cede earlier than
#' the mass. Measured at the same setting, the relative error of the
#' \eqn{\alpha} score against its limit is \eqn{6\times10^{-8}} at
#' \eqn{S = 10^8}, \eqn{4\times10^{-4}} at \eqn{10^{12}} and 1.8 at
#' \eqn{10^{15}}. [betabinom1_distrib()]'s compiled kernel forms the same
#' differences as sums of reciprocals and holds nine figures at
#' \eqn{S = 10^{15}}, so that is the parametrization to use where such a
#' concentration is reachable.
#'
#' # Estimation
#'
#' [fit_distrib()] maximizes the log-likelihood on the link scale. Neither
#' shape is closed form, and the two are strongly correlated at a large
#' concentration, the data then determining their ratio far better than their
#' size.
#'
#' @section Notation:
#' \eqn{\ell} is the log-mass of one observation, \eqn{\alpha, \beta > 0} the
#' two beta shapes, \eqn{S = \alpha+\beta} the concentration, \eqn{n} the trial
#' count, \eqn{B} the beta function and \eqn{\psi^{(k)}} the polygamma
#' functions. \eqn{\eta} is a parameter on the unconstrained scale of its link,
#' with \eqn{\theta = g^{-1}(\eta)}.
#'
#' @return An S7 object of class `BetaBinom2Distrib`, inheriting from
#'   `discrete_distrib`, with `size` the trial count, `distrib_name`
#'   `"betabinom2 [size=n]"`, `dimension` `"univariate"`, `bounds`
#'   `c(0, size)`, `params` `c("alpha", "beta")`, `n_params` `2`,
#'   `params_bounds` the domain \eqn{(0, \infty)} for both, and `link_params`
#'   the two links given here.
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
#' @importFrom linkfunctions7 log_link
#' @importFrom stats rbeta rbinom
#'
#' @examples
#' d <- betabinom2_distrib(size = 10)
#' d
#'
#' # The mass over the support sums to one.
#' th <- list(alpha = 2, beta = 3)
#' sum(distrib_pdf(d, 0:10, th))
#'
#' # The closed-form moments, and the same numbers from the shapes.
#' c(mean = mean(d, th), var = variance(d, th),
#'   closed_mean = 10 * 2 / 5, closed_var = 10 * 2 * 3 * 15 / (25 * 6))
#'
#' # The same law as betabinom1 at mu = alpha / S and sigma = 1 / S.
#' all.equal(distrib_pdf(d, 0:10, th),
#'           distrib_pdf(betabinom1_distrib(size = 10), 0:10,
#'                       list(mu = 0.4, sigma = 0.2)))
#'
#' # At a large concentration the mass is still exact against the binomial
#' # limit, where forming it from two beta functions is not.
#' S <- 1e14
#' c(shipped = distrib_pdf(d, 3, list(alpha = 0.4 * S, beta = 0.6 * S),
#'                         log = TRUE),
#'   two_betas = lchoose(10, 3) +
#'     (lbeta(3 + 0.4 * S, 7 + 0.6 * S) - lbeta(0.4 * S, 0.6 * S)),
#'   binomial = dbinom(3, 10, 0.4, log = TRUE))
#'
#' # Fitting recovers both shapes.
#' set.seed(5)
#' z <- distrib_rng(d, 4000, th)
#' coef(fit_distrib(d, z))
#'
#' @seealso
#' [betabinom1_distrib()] for the same law in a mean proportion and a
#' dispersion, which is the one to model a mean with;
#' [binomial_distrib()] for the limit at a large concentration and
#' [beta1_distrib()] for the mixing law;
#' [fit_distrib()] to estimate the shapes; [check_distrib()] to validate a
#' family of your own against the same battery this one passes;
#' [BetaBinom2Distrib] for the class.
#' @export
betabinom2_distrib <- function(size, link_alpha = log_link(),
                               link_beta = log_link()) {
  if (length(size) != 1L || !is.finite(size) || size < 1 || size != round(size)) {
    stop("'size' must be a single positive integer.", call. = FALSE)
  }
  BetaBinom2Distrib(
    distrib_name = paste0("betabinom2 [size=", size, "]"),
    dimension = "univariate",
    bounds = c(0, size),
    params = c("alpha", "beta"),
    params_interpretation = c(alpha = "shape", beta = "shape"),
    n_params = 2,
    params_bounds = list(alpha = c(0, Inf), beta = c(0, Inf)),
    link_params = list(alpha = link_alpha, beta = link_beta),
    size = size
  )
}
