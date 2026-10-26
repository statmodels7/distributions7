#' @include distrib.R generics.R multivariate.R moments.R
NULL

#' @title Multinomial Distribution Class
#' @name MultinomialDistrib
#'
#' @description
#' This is the S7 class of the multinomial family on \eqn{p} categories and
#' \eqn{n} trials, parametrized by a probability vector carried on a
#' `parameters7` simplex. It inherits from `multivariate_distrib`. Its support
#' is a finite set of points: every vector of non-negative integers summing to
#' \eqn{n}.
#'
#' The class carries two properties beyond the parent's: `size`, the number of
#' trials, fixed at construction; and `param`, the simplex on which the
#' probabilities lie, whose free names become the family's own parameter names
#' prefixed by `probs_`. Build one with [multinomial_distrib()], which validates
#' `n_dim`, `size` and the simplex and fills the properties in. This page
#' documents the raw S7 constructor, which runs only the validators of
#' `multivariate_distrib` and `distrib` and does not check `size` or `param`.
#'
#' @inheritParams multivariate_distrib
#' @param size The number of trials \eqn{n}, a single positive integer stored
#'   as a numeric. It belongs to the object, so an object cannot be reused
#'   across data sets whose trial counts differ.
#' @param param The `parameters7` [parameters7::simplex()] carrying the
#'   probabilities. Its `n_free` is \eqn{p-1}, one fewer than the number of
#'   categories.
#'
#' @return An S7 object of class `MultinomialDistrib`, inheriting from
#'   `multivariate_distrib` and from `distrib`. Beyond `size` and `param` its
#'   properties are the parent's. For an object built by
#'   [multinomial_distrib()], `params` is the simplex's free names prefixed by
#'   `probs_`, `n_params` is \eqn{p-1}, `bounds` is `c(0, size)`, and every
#'   parameter carries an identity link, because its free value is already
#'   unconstrained. The family has no dispersion parameter.
#'
#' @seealso [multinomial_distrib()] to build one;
#'   [binomial_distrib()], which is a coordinate's marginal and the
#'   two-category case; [dirichlet_distrib()] for the continuous family on the
#'   same simplex, which is conjugate to this one;
#'   [parameters7::simplex()] for the chart.
#'
#' @section Methods:
#' The following methods are registered on this class, and each is linked to
#' the page that documents it:
#' - the distribution: [`distrib_pdf()`][distrib_pdf.MultinomialDistrib] and
#'   [`distrib_rng()`][distrib_rng.MultinomialDistrib];
#' - the derivatives in the parameters:
#'   [`distrib_gradient()`][distrib_gradient.MultinomialDistrib],
#'   [`distrib_hessian()`][distrib_hessian.MultinomialDistrib],
#'   [`distrib_expected_hessian()`][distrib_expected_hessian.MultinomialDistrib],
#'   [`distrib_deriv3()`][distrib_deriv3.MultinomialDistrib] and
#'   [`distrib_deriv4()`][distrib_deriv4.MultinomialDistrib];
#' - the moments: [`mean()`][mean.MultinomialDistrib] and
#'   [`variance()`][variance.MultinomialDistrib];
#' - the multivariate summaries:
#'   [`mv_location()`][mv_location.MultinomialDistrib],
#'   [`mv_sigma()`][mv_sigma.MultinomialDistrib],
#'   [`mv_marginal()`][mv_marginal.MultinomialDistrib] and
#'   [`mv_support()`][mv_support.MultinomialDistrib].
#'
#' Every other generic is inherited from [multivariate_distrib()], whose methods
#' signal an error for the distribution function, the quantile function, the
#' response gradient and Hessian, the mixed response derivative, the skewness
#' and the kurtosis.
#'
#' @examples
#' d <- multinomial_distrib(3, size = 5)
#' S7::S7_inherits(d, multivariate_distrib)
#'
#' # Three categories and five trials give two free probabilities.
#' c(categories = d@n_dim, trials = d@size, parameters = d@n_params)
#' d@params
#'
#' # The support is finite and enumerable.
#' nrow(mv_support(d, NULL))
MultinomialDistrib <- S7::new_class("MultinomialDistrib",
  parent = multivariate_distrib,
  properties = list(size = S7::class_numeric, param = S7::class_any)
)

#' The Pieces a Multinomial Derivative Needs
#'
#' @description
#' Evaluates the probability vector and the simplex's first two derivative
#' arrays once, so that a mass function or a derivative method computes them a
#' single time and shares them.
#'
#' @details
#' The probabilities sum to one at every free vector, so differentiating that
#' identity shows that the columns of \eqn{A = \partial p/\partial\eta} sum to
#' zero and so does every second-derivative vector \eqn{B_{\cdot,kl}}. The
#' second of those makes the expected information closed form: under expectation
#' the term carrying \eqn{B} becomes \eqn{n\sum_j B_{j,kl}}, which is zero. On
#' the additive log-ratio chart of [parameters7::simplex()] the
#' log-probabilities are formed directly by [mn_log_probs()], so a probability
#' that underflows keeps its logarithm.
#'
#' @param distrib A [MultinomialDistrib()] object, read for its `param`.
#' @param theta A named list of parameters on the parameter scale, the
#'   simplex's free values in its own order.
#'
#' @return A named list with `prob`, the probability vector of length \eqn{p}
#'   summing to one; `lp`, its logarithm; `alr`, whether the parametrization is
#'   the additive log-ratio chart; `A`, the \eqn{p \times (p-1)} matrix of first
#'   derivatives in the free values; `B`, the list of second-derivative vectors
#'   keyed by unordered tuple; `idx`, those tuples; and `k`, the number of free
#'   values.
#'
#' @seealso [distrib_gradient.MultinomialDistrib()] for the first consumer,
#'   [dir_parts()] for the Dirichlet's counterpart, and
#'   [multinomial_distrib()] for the family.
#'
#' @keywords internal
mn_parts <- function(distrib, theta) {
  s <- distrib@param
  eta <- mv_flat_theta(distrib, theta)
  # on the additive log-ratio chart the log-probabilities are formed directly,
  # so a probability that underflows keeps its logarithm
  alr <- mn_is_alr(s)
  lp <- if (alr) mn_log_probs(s, eta) else NULL
  prob <- if (alr) exp(lp) else as.numeric(parameters7::param_value(s, eta))
  list(prob = prob, lp = if (alr) lp else log(prob), alr = alr,
       A = do.call(cbind, parameters7::param_d1(s, eta)),
       B = parameters7::param_d2(s, eta),
       idx = parameters7::param_tuple_indices(s, 2L),
       k = s@n_free)
}

# --- S7 METHODS IMPLEMENTATION ---

#' @title Multinomial Probability Mass Function
#' @name distrib_pdf.MultinomialDistrib
#' @description
#' Computes the multinomial mass
#' \deqn{P(Y = y) = \dfrac{n!}{\prod_{j=1}^{p} y_j!}\prod_{j=1}^{p} p_j^{y_j},}
#' one value per row of `y`, on the weak compositions of \eqn{n} into \eqn{p}
#' parts. A row with a negative entry, a non-integer entry, or a sum other than
#' \eqn{n} is off the support and returns 0.
#'
#' The factorials are formed through [base::lgamma()] and the product through a
#' matrix multiplication by \eqn{\log p}, so the whole calculation runs on the
#' log scale and a mass far below the smallest representable double still has a
#' finite logarithm. On the additive log-ratio chart \eqn{\log p} is computed as
#' \eqn{\eta_j - \mathrm{LSE}(\eta)} (see [mn_log_probs()]), so a category whose
#' probability underflows and whose count is zero contributes zero, as it
#' should.
#'
#' @param distrib A `MultinomialDistrib` object, from [multinomial_distrib()].
#' @param y A numeric matrix with one row per observation and \eqn{p} columns,
#'   and each row is a vector of non-negative integers summing to the object's
#'   `size`. A single observation may be given as a plain vector and is read as
#'   one row. A row off the support gives a mass of 0, or `-Inf` with
#'   `log = TRUE`.
#' @param theta A named list of the simplex's free values on the parameter
#'   scale, and each value has length 1. A parameter may not vary by
#'   observation here, because the probability vector is one point of the
#'   simplex for the whole sample.
#' @param log Logical of length 1. When `TRUE` the log-mass is returned.
#'   Defaults to `FALSE`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of probabilities in \eqn{[0, 1]} of length
#'   `nrow(y)`, one per observation.
#'
#' @seealso [mv_support.MultinomialDistrib()] for the points on which the mass
#'   sits, [distrib_gradient.MultinomialDistrib()] for the derivatives of the
#'   log-mass, [mv_marginal.MultinomialDistrib()] for a coordinate's binomial
#'   marginal, and [distrib_pdf()] for the generic.
#'
#' @examples
#' d <- multinomial_distrib(3, size = 5)
#' th <- list(probs_alr1 = 0.3, probs_alr2 = -0.2)
#'
#' # The mass over the whole support sums to one, exactly.
#' supp <- mv_support(d, th)
#' c(points = nrow(supp), mass = sum(distrib_pdf(d, supp, th)))
#'
#' # It is stats::dmultinom at the implied probabilities.
#' pr <- mv_location(d, th) / 5
#' all.equal(distrib_pdf(d, supp, th),
#'           apply(supp, 1, function(r) dmultinom(r, prob = pr)))
#'
#' # A row that does not sum to the size, or is not integral, is off the
#' # support.
#' distrib_pdf(d, rbind(c(1, 1, 1), c(-1, 3, 3), c(1.5, 1.5, 2)), th)
S7::method(distrib_pdf, MultinomialDistrib) <- function(distrib, y, theta, log = FALSE, ...) {
  p <- mn_parts(distrib, theta)
  y <- if (is.matrix(y)) y else matrix(y, nrow = 1L)
  n <- distrib@size
  out <- lgamma(n + 1) - rowSums(lgamma(y + 1)) +
    as.numeric(y %*% p$lp)
  bad <- apply(y, 1L, function(r) {
    any(r < 0) || any(r != round(r)) || abs(sum(r) - n) > 1e-8
  })
  out[bad] <- -Inf
  if (log) out else exp(out)
}

#' @title Multinomial Random Generation
#' @name distrib_rng.MultinomialDistrib
#' @description
#' Draws `n` independent multinomial counts by calling [stats::rmultinom()] and
#' transposing, so that one row of the result is one observation: R returns one
#' column per draw, and the multivariate families of the package use one row per
#' observation. The draws depend on `.Random.seed` in the usual way.
#'
#' @param distrib A `MultinomialDistrib` object, from [multinomial_distrib()].
#' @param n A single positive integer, the number of draws. The number of trials
#'   is the object's `size` property.
#' @param theta A named list of the simplex's free values on the parameter
#'   scale, and each value has length 1.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric matrix with `n` rows and \eqn{p} columns, and each row is
#'   a vector of non-negative integers summing to the object's `size`.
#'
#' @seealso [distrib_pdf.MultinomialDistrib()] for the mass,
#'   [mv_support.MultinomialDistrib()] for the points on which it can land,
#'   [fit_distrib()] to estimate the probabilities back from a sample, and
#'   [distrib_rng()] for the generic.
#'
#' @examples
#' d <- multinomial_distrib(3, size = 5)
#' th <- list(probs_alr1 = 0.3, probs_alr2 = -0.2)
#'
#' # Each row is one observation and sums to the trial count.
#' set.seed(1)
#' Y <- distrib_rng(d, 4, th)
#' Y
#' rowSums(Y)
#'
#' # The sample recovers the mean vector and the coordinate variances.
#' set.seed(4)
#' Z <- distrib_rng(d, 2e4, th)
#' rbind(sample = c(colMeans(Z), diag(var(Z))),
#'       theoretical = c(mv_location(d, th), diag(mv_sigma(d, th))))
S7::method(distrib_rng, MultinomialDistrib) <- function(distrib, n, theta, ...) {
  p <- mn_parts(distrib, theta)
  t(stats::rmultinom(n, size = distrib@size, prob = p$prob))
}

#' @title The Support Points of a Multinomial
#' @name mv_support.MultinomialDistrib
#' @description
#' Enumerates the whole support: every vector of \eqn{p} non-negative integers
#' summing to the trial count \eqn{n}, one row per point. These are the weak
#' compositions of \eqn{n} into \eqn{p} parts, supplied by
#' [numericals7::compositions()].
#'
#' Because the support is enumerated, every expectation is an exact sum.
#' [check_distrib()] uses it to test that the mass adds to one, an equality with
#' tolerance `1e-10`; an importance-sampling check compares only against Monte
#' Carlo error, and a normalization wrong by a thousandth would pass it.
#'
#' @param distrib A `MultinomialDistrib` object, from [multinomial_distrib()],
#'   read for its `size` and `n_dim`.
#' @param theta Ignored. The support depends only on \eqn{n} and \eqn{p},
#'   so any value, including `NULL`, gives the same answer.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric matrix with \eqn{\binom{n+p-1}{p-1}} rows and \eqn{p}
#'   columns. The count grows quickly in both arguments: 21 rows at
#'   \eqn{n = 5, p = 3}, and 1001 at \eqn{n = 10, p = 5}.
#'
#' @seealso [numericals7::compositions()] for the enumeration,
#'   [distrib_pdf.MultinomialDistrib()] for the mass on these points,
#'   [check_distrib()], which consumes this, and [mv_support()] for the
#'   generic.
#'
#' @examples
#' d <- multinomial_distrib(3, size = 5)
#' supp <- mv_support(d, NULL)
#' dim(supp)
#' head(supp)
#'
#' # Every row sums to the trial count, and the count of rows is the number of
#' # weak compositions.
#' c(all_sum_to_5 = all(rowSums(supp) == 5),
#'   rows = nrow(supp), formula = choose(5 + 3 - 1, 3 - 1))
#'
#' # The mass over the support adds to one exactly.
#' sum(distrib_pdf(d, supp, list(probs_alr1 = 0.3, probs_alr2 = -0.2))) - 1
S7::method(mv_support, MultinomialDistrib) <- function(distrib, theta, ...) {
  numericals7::compositions(distrib@size, distrib@n_dim)
}

#' @title Multinomial Score
#' @name distrib_gradient.MultinomialDistrib
#' @description
#' Computes the first derivatives of the multinomial log-mass with respect to
#' the free values of the simplex, one value per observation, in closed form.
#' With \eqn{A = \partial p/\partial\eta} the Jacobian of the chart,
#' \deqn{\dfrac{\partial\ell}{\partial\eta_k} = \sum_{j=1}^{p}
#' \dfrac{y_j}{p_j}A_{jk}.} The factorial term of the mass carries no parameter,
#' and the data enter only through the counts. On the additive log-ratio chart,
#' where \eqn{A_{jk} = p_j(\delta_{jk} - p_k)}, the sum is computed as \eqn{y_k
#' - N p_k} with \eqn{N} the total of the row, without dividing by a
#' probability.
#'
#' With `scale = "link"` the generic applies the chain rule for the links of the
#' family to the result of the method. Every parameter takes the identity link
#' here, because its free value is already unconstrained, so the two scales
#' coincide.
#'
#' @param distrib A `MultinomialDistrib` object, from [multinomial_distrib()].
#' @param y A numeric matrix with one row per observation and \eqn{p} columns,
#'   and each row is on the support. A single observation may be given as a
#'   plain vector.
#' @param theta A named list of the simplex's free values on the parameter
#'   scale, and each value has length 1.
#' @param scale One of `"parameter"` (the default) or `"link"`, matched by
#'   [base::match.arg()]. The generic reads it, and this method does not.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A named list with one numeric vector per free value, in the order
#'   that `distrib@params` gives, and each vector has length `nrow(y)`.
#'
#' @section Notation:
#' \eqn{\ell} is the log-mass of one observation, \eqn{p} the probability
#' vector, \eqn{n} the trial count, \eqn{\eta} the free vector of the simplex
#' chart and \eqn{A = \partial p/\partial\eta} its Jacobian.
#'
#' @seealso [distrib_hessian.MultinomialDistrib()] for the second derivatives,
#'   [distrib_expected_hessian.MultinomialDistrib()] for their expectation,
#'   [mn_parts()] for the shared pieces, and [distrib_gradient()] for the
#'   generic.
#'
#' @examples
#' d <- multinomial_distrib(3, size = 5)
#' th <- list(probs_alr1 = 0.3, probs_alr2 = -0.2)
#' set.seed(1)
#' Y <- distrib_rng(d, 4, th)
#' g <- distrib_gradient(d, Y, th)
#' names(g)
#'
#' # It is the derivative of the log-mass, so numDeriv on the summed log-mass
#' # reproduces the summed score.
#' fn <- function(v)
#'   sum(distrib_pdf(d, Y, as.list(stats::setNames(v, d@params)), log = TRUE))
#' rbind(numeric = numDeriv::grad(fn, unlist(th)),
#'       closed = vapply(g, sum, numeric(1)))
#'
#' # The score has mean zero over the support: the first Bartlett identity,
#' # here an exact sum.
#' supp <- mv_support(d, th)
#' w <- distrib_pdf(d, supp, th)
#' vapply(distrib_gradient(d, supp, th), function(v) sum(w * v), numeric(1))
S7::method(distrib_gradient, MultinomialDistrib) <- function(distrib, y, theta,
                                                              scale = c("parameter", "link"), ...) {
  p <- mn_parts(distrib, theta)
  y <- if (is.matrix(y)) y else matrix(y, nrow = 1L)
  if (p$alr) {
    # y_k - N p_k, with N the row's total: no division by a probability
    tot <- rowSums(y)
    return(stats::setNames(
      lapply(seq_len(p$k), function(k) y[, k] - tot * p$prob[k]),
      distrib@params
    ))
  }
  yp <- sweep(y, 2L, p$prob, "/")
  stats::setNames(
    lapply(seq_len(p$k), function(k) as.numeric(yp %*% p$A[, k])),
    distrib@params
  )
}

#' @title Multinomial Observed Hessian
#' @name distrib_hessian.MultinomialDistrib
#' @description
#' Computes the distinct second derivatives of the multinomial log-mass with
#' respect to the simplex's free values, one value per observation, in closed
#' form. With \eqn{A = \partial p/\partial\eta} and
#' \eqn{B_{\cdot,kl} = \partial^2 p/\partial\eta_k\partial\eta_l},
#' \deqn{\ell^{(\eta_k\eta_l)} = \sum_{j=1}^{p}\left(\dfrac{y_j}{p_j}B_{j,kl}
#'       - \dfrac{y_j}{p_j^2}A_{jk}A_{jl}\right).}
#' On the additive log-ratio chart the two terms combine into
#' \eqn{-N(\delta_{kl}p_k - p_kp_l)}, with \eqn{N} the total of the row, the
#' second cumulant of the categorical law, which the method computes directly
#' and without dividing by a probability; for a row on the support \eqn{N = n},
#' so the observed matrix equals its expectation. On another chart the two terms
#' need not combine in this way.
#'
#' @param distrib A `MultinomialDistrib` object, from [multinomial_distrib()].
#' @param y A numeric matrix with one row per observation and \eqn{p} columns,
#'   and each row is on the support. A single observation may be given as a
#'   plain vector.
#' @param theta A named list of the simplex's free values on the parameter
#'   scale, and each value has length 1.
#' @param scale One of `"parameter"` (the default) or `"link"`, matched by
#'   [base::match.arg()]. The generic reads it, and this method does not.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A named list of numeric vectors, one per distinct entry of the
#'   symmetric matrix and named as [hess_names()] names them, and each vector
#'   has length `nrow(y)`. For \eqn{p} categories there are \eqn{p(p-1)/2} of
#'   them.
#'
#' @section Notation:
#' \eqn{\ell} is the log-mass of one observation, \eqn{p} the probability
#' vector, \eqn{\eta} the free vector of the simplex chart, \eqn{A} its
#' Jacobian and \eqn{B} its second-derivative arrays.
#'
#' @seealso [distrib_gradient.MultinomialDistrib()] for the score,
#'   [distrib_expected_hessian.MultinomialDistrib()] for the expectation of
#'   this quantity, [distrib_deriv3.MultinomialDistrib()] for the order above,
#'   and [distrib_hessian()] for the generic.
#'
#' @examples
#' d <- multinomial_distrib(3, size = 5)
#' th <- list(probs_alr1 = 0.3, probs_alr2 = -0.2)
#' set.seed(1)
#' Y <- distrib_rng(d, 4, th)
#' h <- distrib_hessian(d, Y, th)
#' names(h)
#'
#' # numDeriv on the summed log-mass reproduces the summed matrix.
#' fn <- function(v)
#'   sum(distrib_pdf(d, Y, as.list(stats::setNames(v, d@params)), log = TRUE))
#' H <- numDeriv::hessian(fn, unlist(th))
#' rbind(numeric = c(H[1, 1], H[2, 2], H[1, 2]),
#'       closed = vapply(h, sum, numeric(1)))
#'
#' # The mass-weighted sum over the support is the expected Hessian, exactly.
#' supp <- mv_support(d, th)
#' w <- distrib_pdf(d, supp, th)
#' rbind(summed = vapply(distrib_hessian(d, supp, th),
#'                       function(v) sum(w * v), numeric(1)),
#'       expected = vapply(distrib_expected_hessian(d, Y, th),
#'                         function(v) v[1], numeric(1)))
S7::method(distrib_hessian, MultinomialDistrib) <- function(distrib, y, theta,
                                                             scale = c("parameter", "link"), ...) {
  p <- mn_parts(distrib, theta)
  y <- if (is.matrix(y)) y else matrix(y, nrow = 1L)
  nm <- hess_names(distrib@params)
  pr <- hess_pairs(distrib@params)
  pos <- stats::setNames(seq_along(distrib@params), distrib@params)
  if (p$alr) {
    # -N (p_a delta_ab - p_a p_b), the second cumulant of the categorical law
    tot <- rowSums(y)
    return(stats::setNames(lapply(seq_along(nm), function(m) {
      a <- pos[[pr[[m]][1]]]
      b <- pos[[pr[[m]][2]]]
      -tot * mn_cumulant(p$prob, c(a, b))
    }), nm))
  }
  yp <- sweep(y, 2L, p$prob, "/")
  yp2 <- sweep(y, 2L, p$prob^2, "/")

  stats::setNames(lapply(seq_along(nm), function(m) {
    a <- pos[[pr[[m]][1]]]
    b <- pos[[pr[[m]][2]]]
    bi <- which(vapply(p$idx, function(t) {
      length(t) == 2L && identical(sort(as.integer(t)),
                                   sort(as.integer(c(a, b))))
    }, logical(1)))[1L]
    as.numeric(yp %*% p$B[[bi]]) -
      as.numeric(yp2 %*% (p$A[, a] * p$A[, b]))
  }), nm)
}

#' @title Multinomial Expected Hessian
#' @name distrib_expected_hessian.MultinomialDistrib
#' @description
#' Returns the expectation of the observed Hessian under the model, in closed
#' form and with no quadrature, simulation or summation over the support:
#' \deqn{\mathbb{E}[\ell^{(\eta_k\eta_l)}] = -n\sum_{j=1}^{p}
#' \dfrac{A_{jk}A_{jl}}{p_j},} with \eqn{A = \partial p/\partial\eta}; on the
#' additive log-ratio chart it is computed as \eqn{-n(\delta_{kl}p_k -
#' p_kp_l)}, as shown below.
#'
#' @details
#' # The second-derivative term
#'
#' Under the model \eqn{\mathbb{E}[y_j] = n p_j}, so the term of the observed
#' Hessian carrying \eqn{B} becomes \eqn{n\sum_j B_{j,kl}}. The probabilities
#' sum to one at every free vector, so every derivative of that sum is zero and
#' the whole term drops out, leaving the display above. The same argument, one
#' order lower, gives the score mean zero.
#'
#' # The information on the additive log-ratio chart
#'
#' On the default chart, [parameters7::simplex()]'s additive log-ratio
#' coordinates, the Jacobian is \eqn{A_{jk} = p_j(\delta_{jk} - p_k)}, so
#' \deqn{\sum_j \dfrac{A_{jk}A_{jl}}{p_j}
#'   = \sum_j p_j(\delta_{jk}-p_k)(\delta_{jl}-p_l)
#'   = \delta_{kl}p_k - p_k p_l,}
#' and the expected information \eqn{-\mathbb{E}[\ell^{(\eta_k\eta_l)}]} is
#' exactly \eqn{\operatorname{Cov}(Y_k, Y_l)}, the leading
#' \eqn{(p-1)\times(p-1)} block of [mv_sigma.MultinomialDistrib()]. The identity
#' belongs to that chart; a different simplex parametrization has a different
#' Jacobian.
#'
#' @param distrib A `MultinomialDistrib` object, from [multinomial_distrib()].
#' @param y A numeric matrix of observations. Only its row count is read
#'   through [n_obs()], because the expectation does not depend on the data; the
#'   values are ignored.
#' @param theta A named list of the simplex's free values on the parameter
#'   scale, and each value has length 1.
#' @param scale One of `"parameter"` (the default) or `"link"`, matched by
#'   [base::match.arg()]. The generic reads it, and this method does not.
#' @param method `NULL` or an object inheriting from [ExpectMethod]. The generic
#'   applies an object before dispatching on the family, so this method
#'   receives `NULL`, and it returns the closed form of the family. See
#'   [expected_derivative_methods].
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A named list of numeric vectors, one per distinct entry of the
#'   symmetric matrix and named as [hess_names()] names them, and each vector
#'   has length `n_obs(distrib, y)` and is constant along it.
#'
#' @section Notation:
#' \eqn{\ell} is the log-mass of one observation, \eqn{p} the probability
#' vector, \eqn{n} the trial count, \eqn{\eta} the free vector of the simplex
#' chart and \eqn{A = \partial p/\partial\eta} its Jacobian.
#'
#' @seealso [distrib_hessian.MultinomialDistrib()] for the quantity whose
#'   expectation this is, [mv_sigma.MultinomialDistrib()] for the covariance
#'   with which it coincides, [mv_support.MultinomialDistrib()] for the exact
#'   sum against which it can be checked, and [distrib_expected_hessian()] for
#'   the generic.
#'
#' @examples
#' d <- multinomial_distrib(3, size = 5)
#' th <- list(probs_alr1 = 0.3, probs_alr2 = -0.2)
#' eh <- distrib_expected_hessian(d, matrix(0, 1, 3), th)
#' vapply(eh, function(v) v[1], numeric(1))
#'
#' # The exact sum over the support gives the same three numbers.
#' supp <- mv_support(d, th)
#' w <- distrib_pdf(d, supp, th)
#' vapply(distrib_hessian(d, supp, th), function(v) sum(w * v), numeric(1))
#'
#' # On the default chart the information is the covariance of the first
#' # p - 1 counts.
#' S <- mv_sigma(d, th)
#' rbind(information = -vapply(eh, function(v) v[1], numeric(1)),
#'       covariance = c(S[1, 1], S[2, 2], S[1, 2]))
#'
#' # The method argument is inert, because the expectation is exact.
#' identical(eh, distrib_expected_hessian(d, matrix(0, 1, 3), th,
#'                                        method = expect_mc(nsim = 50)))
S7::method(distrib_expected_hessian, MultinomialDistrib) <- function(distrib, y, theta,
                                                                      scale = c("parameter", "link"),
                                                                      method = NULL, ...) {
  p <- mn_parts(distrib, theta)
  n <- n_obs(distrib, y)
  nm <- hess_names(distrib@params)
  pr <- hess_pairs(distrib@params)
  pos <- stats::setNames(seq_along(distrib@params), distrib@params)

  stats::setNames(lapply(seq_along(nm), function(m) {
    a <- pos[[pr[[m]][1]]]
    b <- pos[[pr[[m]][2]]]
    v <- if (p$alr) {
      -distrib@size * mn_cumulant(p$prob, c(a, b))
    } else {
      -distrib@size * sum(p$A[, a] * p$A[, b] / p$prob)
    }
    rep(v, n)
  }), nm)
}

#' @title Multinomial Mean Vector
#' @name mv_location.MultinomialDistrib
#' @description
#' Returns \eqn{\mathbb{E}[Y] = n p}, the trial count times the probability
#' vector, so the result sums to \eqn{n} and not to one. Divide by the object's
#' `size` to recover the probabilities themselves.
#' [mean.MultinomialDistrib()] delegates here.
#'
#' @param distrib A `MultinomialDistrib` object, from [multinomial_distrib()].
#' @param theta A named list of the simplex's free values on the parameter
#'   scale, and each value has length 1.
#'
#' @return A numeric vector of length \eqn{p}, without names, summing to the
#'   object's `size`.
#'
#' @seealso [mv_sigma.MultinomialDistrib()] for the covariance,
#'   [mean.MultinomialDistrib()], which calls this, and [mv_location()] for the
#'   generic.
#'
#' @examples
#' d <- multinomial_distrib(3, size = 5)
#' th <- list(probs_alr1 = 0.3, probs_alr2 = -0.2)
#' m <- mv_location(d, th)
#' m
#' sum(m)
#'
#' # Dividing by the trial count gives the probabilities, which sum to one.
#' m / d@size
S7::method(mv_location, MultinomialDistrib) <- function(distrib, theta) {
  distrib@size * mn_parts(distrib, theta)$prob
}

#' @title Multinomial Covariance Matrix
#' @name mv_sigma.MultinomialDistrib
#' @description
#' Returns the covariance matrix of the counts,
#' \deqn{\operatorname{Cov}(Y_i, Y_j) = n(\delta_{ij}p_i - p_i p_j),}
#' so the coordinate variances are \eqn{np_j(1-p_j)}, the binomial ones, and
#' every covariance is negative, because the counts share a fixed total.
#'
#' The matrix is singular, because the coordinates sum to \eqn{n}, so the vector
#' of ones is in its null space and the rank is \eqn{p-1}. Its leading
#' \eqn{(p-1)\times(p-1)} block is the expected information on the default
#' chart; see [distrib_expected_hessian.MultinomialDistrib()].
#'
#' @param distrib A `MultinomialDistrib` object, from [multinomial_distrib()].
#' @param theta A named list of the simplex's free values on the parameter
#'   scale, and each value has length 1.
#'
#' @return A symmetric \eqn{p \times p} numeric matrix of rank \eqn{p-1},
#'   without dimnames.
#'
#' @seealso [mv_location.MultinomialDistrib()] for the mean,
#'   [variance.MultinomialDistrib()], which calls this,
#'   [mv_marginal.MultinomialDistrib()] for a coordinate's binomial law, and
#'   [mv_sigma()] for the generic.
#'
#' @examples
#' d <- multinomial_distrib(3, size = 5)
#' th <- list(probs_alr1 = 0.3, probs_alr2 = -0.2)
#' S <- mv_sigma(d, th)
#' round(S, 6)
#'
#' # The matrix is singular by construction: its rank is p - 1, and the ones
#' # vector is in its null space.
#' c(rank = qr(S)$rank, dim = ncol(S))
#' round(S %*% rep(1, 3), 12)
#'
#' # The diagonal is the binomial variance of each category.
#' pr <- mv_location(d, th) / 5
#' rbind(diagonal = diag(S), binomial = 5 * pr * (1 - pr))
S7::method(mv_sigma, MultinomialDistrib) <- function(distrib, theta) {
  p <- mn_parts(distrib, theta)$prob
  distrib@size * (diag(p, nrow = length(p)) - tcrossprod(p))
}

#' @title Mean of a Multinomial
#' @name mean.MultinomialDistrib
#' @description
#' Returns the mean count vector \eqn{\mathbb{E}[Y] = np} by delegating to
#' [mv_location.MultinomialDistrib()]. The result sums to the trial count, not
#' to one.
#'
#' @param x A `MultinomialDistrib` object, from [multinomial_distrib()]. The
#'   argument is named `x` because the generic is [base::mean()].
#' @param theta A named list of the simplex's free values on the parameter
#'   scale, and each value has length 1.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of length \eqn{p} summing to the object's `size`.
#'
#' @seealso [variance.MultinomialDistrib()] for the covariance,
#'   [mv_location.MultinomialDistrib()], which this calls, and
#'   [mean.distrib()] for the generic's default.
#' @keywords internal
S7::method(mean, MultinomialDistrib) <- function(x, theta, ...) {
  mv_location(x, theta)
}

#' @title Variance of a Multinomial
#' @name variance.MultinomialDistrib
#' @description
#' Returns the covariance matrix of the counts by delegating to
#' [mv_sigma.MultinomialDistrib()]. It is singular because the counts sum to the
#' trial count, so a computation that inverts a covariance works with the
#' marginals or with the free vector of the simplex.
#'
#' @param x A `MultinomialDistrib` object, from [multinomial_distrib()]. The
#'   argument is named `x` for consistency with [base::mean()], whose signature
#'   the moment generics follow.
#' @param theta A named list of the simplex's free values on the parameter
#'   scale, and each value has length 1.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A symmetric \eqn{p \times p} numeric matrix of rank \eqn{p-1}.
#'
#' @seealso [mean.MultinomialDistrib()] for the mean vector,
#'   [mv_sigma.MultinomialDistrib()], which this calls, and [variance()] for
#'   the generic.
#' @keywords internal
S7::method(variance, MultinomialDistrib) <- function(x, theta, ...) {
  mv_sigma(x, theta)
}

#' @title Multinomial Marginal
#' @name mv_marginal.MultinomialDistrib
#' @description
#' Returns the law of one coordinate, which is
#' \eqn{\mathrm{Binomial}(n, p_j)}: the other categories collapse into a single
#' failure, and the trial count is unchanged. The returned object carries the
#' same `size` as the multinomial and its `mu` is \eqn{p_j}.
#'
#' A vector of several coordinates signals an error: a sub-vector is again
#' multinomial only once the remaining outcomes are collapsed into a category of
#' their own, which is a different object from the sub-vector, and the error
#' message gives this reason. For the same reason [plot.multivariate_distrib()],
#' whose panels are marginals of two coordinates, signals that error for this
#' family.
#'
#' @param distrib A `MultinomialDistrib` object, from [multinomial_distrib()].
#' @param theta A named list of the simplex's free values on the parameter
#'   scale, and each value has length 1.
#' @param which A single integer in \eqn{1, \dots, p}, the category whose
#'   marginal is returned. A vector of two or more distinct categories signals
#'   an error explaining why; an empty, duplicated, non-integer or out-of-range
#'   `which` is rejected earlier by [mv_marginal()].
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A named list with `distrib`, a [binomial_distrib()] object at the
#'   same `size`, and `theta`, a list holding the marginal's `mu`.
#'
#' @seealso [binomial_distrib()] for the family returned,
#'   [mv_sigma.MultinomialDistrib()] for the coordinate variances, and
#'   [mv_marginal()] for the generic.
#'
#' @examples
#' d <- multinomial_distrib(3, size = 5)
#' th <- list(probs_alr1 = 0.3, probs_alr2 = -0.2)
#' m <- mv_marginal(d, th, which = 2)
#' m$distrib
#' m$theta
#'
#' # The marginal's mean and variance are the second coordinate's.
#' rbind(marginal = c(mean(m$distrib, m$theta), variance(m$distrib, m$theta)),
#'       joint = c(mv_location(d, th)[2], mv_sigma(d, th)[2, 2]))
#'
#' # Every marginal carries the same trial count.
#' vapply(1:3, function(j) mv_marginal(d, th, which = j)$distrib@size,
#'        numeric(1))
S7::method(mv_marginal, MultinomialDistrib) <- function(distrib, theta, which, ...) {
  if (length(which) != 1L) {
    stop(paste0(
      "A multinomial marginal is returned one coordinate at a time. Several\n",
      "  coordinates are again multinomial, but only after the remaining\n",
      "  outcomes are collapsed into a category of their own, which is a\n",
      "  different object from a sub-vector."
    ), call. = FALSE)
  }
  list(distrib = binomial_distrib(size = distrib@size),
       theta = list(mu = mn_parts(distrib, theta)$prob[which]))
}

# --- CONSTRUCTOR WRAPPER ---

#' Multinomial Distribution
#'
#' @description
#' Builds the distribution object for the multinomial family on \eqn{p}
#' categories and \eqn{n} trials, parametrized by a probability vector carried
#' on a `parameters7` simplex. The returned object carries closed-form
#' derivatives of the log-mass to fourth order, a closed-form expected
#' information and closed-form expectations of the third and fourth derivatives,
#' and its support is small enough to enumerate, so every expectation can be
#' checked as an exact sum.
#'
#' @param n_dim The number of categories \eqn{p}, a single whole number of at
#'   least 2. Any other value signals an error naming the argument.
#' @param size The number of trials \eqn{n}, a single positive integer. It is a
#'   constant of the distribution and not a parameter, as for
#'   [binomial_distrib()], so an object cannot be reused across data sets whose
#'   trial counts differ.
#' @param probs A `parameters7` parameter producing \eqn{p} coordinates that
#'   sum to one, normally a [parameters7::simplex()]. Defaults to
#'   `parameters7::simplex(n_dim)`. An object that is not a `parameters7`
#'   parameter, or that produces a different number of coordinates, signals an
#'   error. Its free values are unconstrained already and carry the identity
#'   link.
#'
#' @details
#' # The parametrization
#'
#' The mass on the weak compositions of \eqn{n} into \eqn{p} parts is
#' \deqn{P(Y = y) = \frac{n!}{\prod_{j=1}^{p} y_j!}\prod_{j=1}^{p} p_j^{y_j},}
#' with
#' \deqn{\mathbb{E}[Y_j] = n p_j, \qquad
#'       \operatorname{Var}(Y_j) = n p_j (1 - p_j), \qquad
#'       \operatorname{Cov}(Y_j, Y_k) = -n p_j p_k.}
#' The family has no dispersion parameter, and its parameters are the
#' probabilities, of which \eqn{p-1} are free.
#'
#' The probabilities are carried by a `parameters7` simplex and flattened into
#' scalars with identity links, as a covariance is for the multivariate
#' gaussian. The constraint that they be positive and sum to one is carried by
#' the parameter, which a scalar link could not express.
#'
#' # The support
#'
#' The family is multivariate and discrete, so its support is a finite set of
#' points. [mv_support.MultinomialDistrib()] enumerates them, which makes an
#' expectation an exact sum and lets [check_distrib()] test the total mass by
#' addition. An importance-sampling check compares only against Monte Carlo
#' error, and a normalization wrong by a thousandth would pass it.
#'
#' The count of support points is \eqn{\binom{n+p-1}{p-1}} and grows quickly,
#' so the enumeration is a validation tool rather than a fitting route.
#'
#' # Score and information
#'
#' With \eqn{A = \partial p/\partial\eta} the chart's Jacobian,
#' \deqn{\dfrac{\partial\ell}{\partial\eta_k} = \sum_j \dfrac{y_j}{p_j}A_{jk},
#'       \qquad
#'       \mathbb{E}[\ell^{(\eta_k\eta_l)}]
#'         = -n\sum_j \dfrac{A_{jk}A_{jl}}{p_j}.}
#' The expected form is closed because \eqn{\mathbb{E}[y_j] = np_j} turns the
#' second-derivative term into \eqn{n\sum_j B_{j,kl}}, which vanishes: the
#' probabilities sum to one, so every derivative of their sum is zero.
#'
#' On the default additive log-ratio chart the resulting information is exactly
#' the covariance of the first \eqn{p-1} counts, since \eqn{A_{jk} =
#' p_j(\delta_{jk}-p_k)} collapses the sum to \eqn{\delta_{kl}p_k - p_kp_l}. On
#' that chart every derivative of the log-mass of order two or more is \eqn{-N}
#' times a cumulant of the categorical law (see [mn_cumulant()]), a sum of
#' products of probabilities that the methods compute without dividing by one,
#' so a category whose probability is tiny or underflows gives finite
#' derivatives.
#'
#' # Marginals and conjugacy
#'
#' Coordinate \eqn{j} is \eqn{\mathrm{Binomial}(n, p_j)}, because the other
#' categories collapse into a single failure, so
#' [mv_marginal.MultinomialDistrib()] returns an object. The Dirichlet is the
#' conjugate prior for the probability vector, and
#' [dirichlet_distrib()] is written on the same simplex.
#'
#' # Estimation
#'
#' [fit_distrib()] maximizes the log-likelihood over the simplex's free values,
#' whose links are the identity, so the two scales coincide. The maximum
#' likelihood estimate of \eqn{p_j} is the pooled sample proportion.
#'
#' @section Notation:
#' \eqn{\ell} is the log-mass of one observation, \eqn{p} the probability
#' vector, \eqn{n} the trial count, \eqn{\eta} the free vector of the simplex
#' chart, \eqn{A = \partial p/\partial\eta} its Jacobian and \eqn{B} its
#' second-derivative arrays.
#'
#' @return An S7 object of class `MultinomialDistrib`, inheriting from
#'   `multivariate_distrib`, with `size` the trial count, `param` the simplex
#'   given here, `distrib_name` `"multinomial [<p>d, size=<n>, probs=<name>]"`
#'   with `<name>` the `param_name` of the probability parameter,
#'   `dimension` `"multivariate"`, `n_dim` \eqn{p}, `bounds` `c(0, size)`,
#'   `params` the simplex's free names prefixed by `probs_`, `n_params`
#'   \eqn{p-1}, and `link_params` the identity for each.
#'
#' @references
#' Johnson, N. L., Kotz, S. and Balakrishnan, N. (1997).
#' *Discrete Multivariate Distributions*, Chapter 35. Wiley, New York.
#'
#' Agresti, A. (2013). *Categorical Data Analysis*, 3rd edition, Section 1.2.
#' Wiley, Hoboken.
#'
#' @importFrom linkfunctions7 identity_link
#' @importFrom stats rmultinom
#'
#' @examples
#' d <- multinomial_distrib(3, size = 5)
#' d
#'
#' # Three categories give two free probabilities, and the family has no
#' # dispersion parameter.
#' th <- list(probs_alr1 = 0.3, probs_alr2 = -0.2)
#' mv_location(d, th)
#'
#' # The support is a finite set of points, so the mass sums exactly.
#' supp <- mv_support(d, th)
#' c(points = nrow(supp), mass = sum(distrib_pdf(d, supp, th)))
#'
#' # The covariance is singular, and its diagonal holds the binomial variances.
#' pr <- mv_location(d, th) / 5
#' S <- mv_sigma(d, th)
#' c(rank = qr(S)$rank, dim = 3)
#' rbind(diagonal = diag(S), binomial = 5 * pr * (1 - pr))
#'
#' # On the default chart the expected information is the covariance of the
#' # first p - 1 counts.
#' eh <- distrib_expected_hessian(d, matrix(0, 1, 3), th)
#' rbind(information = -vapply(eh, function(v) v[1], numeric(1)),
#'       covariance = c(S[1, 1], S[2, 2], S[1, 2]))
#'
#' # Fitting recovers the probabilities.
#' set.seed(4)
#' coef(fit_distrib(d, distrib_rng(d, 1500, th)))
#'
#' @seealso
#' [binomial_distrib()] for a category's marginal and the two-category case;
#' [dirichlet_distrib()] for the conjugate family on the same simplex;
#' [betabinom1_distrib()] for the overdispersed two-category count;
#' [parameters7::simplex()] for the chart that carries the probabilities;
#' [numericals7::compositions()] for the support enumeration; [fit_distrib()] to
#' estimate the probabilities; [check_distrib()] to validate a user-defined
#' family against the same battery that this one passes; [MultinomialDistrib]
#' for the class.
#' @export
multinomial_distrib <- function(n_dim, size,
                                probs = parameters7::simplex(n_dim)) {
  if (!is.numeric(n_dim) || length(n_dim) != 1L || !is.finite(n_dim) ||
      n_dim != round(n_dim) || n_dim < 2) {
    stop("'n_dim' must be a single integer of at least 2.", call. = FALSE)
  }
  p <- as.integer(n_dim)
  if (!is.numeric(size) || length(size) != 1L || !is.finite(size) ||
      size < 1 || size != round(size)) {
    stop("'size' must be a single positive integer.", call. = FALSE)
  }
  if (!S7::S7_inherits(probs, parameters7::parameter)) {
    stop("'probs' must be a parameters7 parameter.", call. = FALSE)
  }
  if (length(parameters7::param_value(probs, numeric(probs@n_free))) != p) {
    stop("The probability parameter does not have the stated dimension.",
         call. = FALSE)
  }

  params <- paste0("probs_", probs@free_names)
  MultinomialDistrib(
    distrib_name = sprintf("multinomial [%dd, size=%g, probs=%s]", p, size,
                           probs@param_name),
    dimension = "multivariate", n_dim = p, bounds = c(0, size),
    params = params,
    params_interpretation = stats::setNames(
      rep("probability", probs@n_free), params
    ),
    n_params = probs@n_free,
    params_bounds = stats::setNames(
      rep(list(c(-Inf, Inf)), probs@n_free), params
    ),
    link_params = stats::setNames(
      rep(list(linkfunctions7::identity_link()), probs@n_free), params
    ),
    size = as.numeric(size),
    param = probs
  )
}
