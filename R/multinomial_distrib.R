#' @include distrib.R generics.R multivariate.R moments.R
NULL

#' @title Multinomial Distribution Class
#' @name MultinomialDistrib
#'
#' @description
#' The S7 class of the multinomial family on \eqn{p} categories and \eqn{n}
#' trials, parametrized by a probability vector carried on a `parameters7`
#' simplex. It inherits from `multivariate_distrib`, and it is the first
#' multivariate family here whose support is a **finite set of points** rather
#' than a region: every vector of non-negative integers summing to \eqn{n}.
#' The nine methods listed below are registered on it in this file and two more
#' in `mv_higher.R`.
#'
#' The class carries two properties beyond the parent's: `size`, the number of
#' trials, fixed at construction; and `param`, the
#' simplex the probabilities lie on, whose free names become the family's own
#' parameter names prefixed by `probs_`. Build one with
#' [multinomial_distrib()], which validates all three and fills the properties
#' in. This page documents the raw S7 constructor, which validates none of
#' that.
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
#'   parameter carries an identity link, its free value being unconstrained
#'   already. There is **no dispersion parameter**: the probabilities are all
#'   there is.
#'
#' @seealso [multinomial_distrib()] to build one;
#'   [binomial_distrib()], which is a coordinate's marginal and the
#'   two-category case; [dirichlet_distrib()] for the continuous family on the
#'   same simplex, which is conjugate to this one;
#'   [parameters7::simplex()] for the chart.
#'
#' @section Methods:
#' Registered on this class in this file:
#'   [`distrib_expected_hessian()`][distrib_expected_hessian.MultinomialDistrib],
#'   [`distrib_gradient()`][distrib_gradient.MultinomialDistrib],
#'   [`distrib_hessian()`][distrib_hessian.MultinomialDistrib],
#'   [`distrib_pdf()`][distrib_pdf.MultinomialDistrib],
#'   [`distrib_rng()`][distrib_rng.MultinomialDistrib],
#'   [`mv_location()`][mv_location.MultinomialDistrib],
#'   [`mv_marginal()`][mv_marginal.MultinomialDistrib],
#'   [`mv_sigma()`][mv_sigma.MultinomialDistrib],
#'   [`mv_support()`][mv_support.MultinomialDistrib],
#'   [`mean()`][mean.MultinomialDistrib] and
#'   [`variance()`][variance.MultinomialDistrib].
#'
#' Two more are registered in `mv_higher.R`:
#'   [`distrib_deriv3()`][distrib_deriv3.MultinomialDistrib] and
#'   [`distrib_deriv4()`][distrib_deriv4.MultinomialDistrib].
#'
#' Everything else is inherited from [multivariate_distrib()], which refuses
#' the distribution function, the quantile and the response derivatives.
#'
#' @examples
#' d <- multinomial_distrib(3, size = 5)
#' S7::S7_inherits(d, multivariate_distrib)
#'
#' # Three categories and five trials, so two free probabilities.
#' c(categories = d@n_dim, trials = d@size, parameters = d@n_params)
#' d@params
#'
#' # The support is finite and enumerable, which is what the class is for.
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
#' second of those makes the expected information closed form: under
#' expectation the term carrying \eqn{B} becomes \eqn{n\sum_j B_{j,kl}}, which
#' is zero.
#'
#' @param distrib A [MultinomialDistrib()] object, read for its `param`.
#' @param theta A named list of parameters on the parameter scale, the
#'   simplex's free values in its own order.
#'
#' @return A named list with `prob`, the probability vector of length \eqn{p}
#'   summing to one; `A`, the \eqn{p \times (p-1)} matrix of first derivatives
#'   in the free values; `B`, the list of second-derivative vectors keyed by
#'   unordered tuple; `idx`, those tuples; and `k`, the number of free values.
#'
#' @seealso [distrib_gradient.MultinomialDistrib()] for the first consumer,
#'   [dir_parts()] for the Dirichlet's counterpart, and
#'   [multinomial_distrib()] for the family.
#'
#' @keywords internal
mn_parts <- function(distrib, theta) {
  s <- distrib@param
  eta <- mv_flat_theta(distrib, theta)
  list(prob = as.numeric(parameters7::param_value(s, eta)),
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
#' finite logarithm.
#'
#' @param distrib A `MultinomialDistrib` object, from [multinomial_distrib()].
#' @param y A numeric matrix with one row per observation and \eqn{p} columns,
#'   each row a vector of non-negative integers summing to the object's `size`.
#'   A single observation may be given as a plain vector and is read as one
#'   row. A row off the support gives a mass of 0, or `-Inf` with
#'   `log = TRUE`.
#' @param theta A named list of the simplex's free values on the parameter
#'   scale, each of length 1. A parameter may not vary by observation here, the
#'   probability vector being one point of the simplex for the whole sample.
#' @param log Logical of length 1. When `TRUE` the log-mass is returned.
#'   Defaults to `FALSE`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of probabilities in \eqn{[0, 1]} of length
#'   `nrow(y)`, one per observation.
#'
#' @seealso [mv_support.MultinomialDistrib()] for the points the mass sits on,
#'   [distrib_gradient.MultinomialDistrib()] for the derivatives of the
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
    as.numeric(y %*% log(p$prob))
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
#' transposing, so that one **row** of the result is one observation. R returns
#' one column per draw; the package's convention throughout the multivariate
#' families is one row per observation, and the transpose is the whole of the
#' difference. The draws depend on `.Random.seed` in the usual way.
#'
#' @param distrib A `MultinomialDistrib` object, from [multinomial_distrib()].
#' @param n A single positive integer, the number of draws. Note that the
#'   number of **trials** is the object's `size` property, not this argument.
#' @param theta A named list of the simplex's free values on the parameter
#'   scale, each of length 1.
#'
#' @return A numeric matrix with `n` rows and \eqn{p} columns, each row a
#'   vector of non-negative integers summing to the object's `size`.
#'
#' @seealso [distrib_pdf.MultinomialDistrib()] for the mass,
#'   [mv_support.MultinomialDistrib()] for the points it can land on,
#'   [fit_distrib()] to estimate the probabilities back from a sample, and
#'   [distrib_rng()] for the generic.
#'
#' @examples
#' d <- multinomial_distrib(3, size = 5)
#' th <- list(probs_alr1 = 0.3, probs_alr2 = -0.2)
#'
#' # One row per observation, each summing to the trial count.
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
S7::method(distrib_rng, MultinomialDistrib) <- function(distrib, n, theta) {
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
#' Having the support in hand is what turns every expectation into an **exact
#' sum**. [check_distrib()] uses it to test that the mass adds to one and that
#' the closed-form information matches the summed observed Hessian, both to
#' machine precision; an importance-sampling check could only ever compare
#' against Monte Carlo error, and a normalization wrong by a thousandth would
#' pass it.
#'
#' @param distrib A `MultinomialDistrib` object, from [multinomial_distrib()],
#'   read for its `size` and `n_dim`.
#' @param theta Ignored. The support is a property of \eqn{n} and \eqn{p}
#'   alone, so any value, `NULL` included, gives the same answer.
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
#' # The mass over it adds to one exactly, which no sampling check could say.
#' sum(distrib_pdf(d, supp, list(probs_alr1 = 0.3, probs_alr2 = -0.2))) - 1
S7::method(mv_support, MultinomialDistrib) <- function(distrib, theta, ...) {
  numericals7::compositions(distrib@size, distrib@n_dim)
}

#' @title Multinomial Score
#' @name distrib_gradient.MultinomialDistrib
#' @description
#' Computes the first derivatives of the multinomial log-mass with respect to
#' the simplex's free values, one value per observation, in closed form. With
#' \eqn{A = \partial p/\partial\eta} the chart's own Jacobian,
#' \deqn{\dfrac{\partial\ell}{\partial\eta_k}
#'       = \sum_{j=1}^{p} \dfrac{y_j}{p_j}A_{jk}.}
#' The factorial term of the mass carries no parameter, so it contributes
#' nothing, and the data enter only through the counts themselves.
#'
#' With `scale = "link"` the generic applies the chain rule for the links the
#' family carries. Every parameter rides the identity here, its free value
#' being unconstrained already, so the two scales coincide.
#'
#' @param distrib A `MultinomialDistrib` object, from [multinomial_distrib()].
#' @param y A numeric matrix with one row per observation and \eqn{p} columns,
#'   each row on the support. A single observation may be given as a plain
#'   vector.
#' @param theta A named list of the simplex's free values on the parameter
#'   scale, each of length 1.
#' @param scale One of `"parameter"` (the default) or `"link"`, matched by
#'   [base::match.arg()]. Read by the generic, not by this method.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A named list with one numeric vector per free value, in the order
#'   `distrib@params` gives, each of length `nrow(y)`.
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
#' Both terms carry the data, so the observed matrix differs from its
#' expectation at every observation.
#'
#' @param distrib A `MultinomialDistrib` object, from [multinomial_distrib()].
#' @param y A numeric matrix with one row per observation and \eqn{p} columns,
#'   each row on the support. A single observation may be given as a plain
#'   vector.
#' @param theta A named list of the simplex's free values on the parameter
#'   scale, each of length 1.
#' @param scale One of `"parameter"` (the default) or `"link"`, matched by
#'   [base::match.arg()]. Read by the generic, not by this method.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A named list of numeric vectors, one per distinct entry of the
#'   symmetric matrix and named as [hess_names()] names them, each of length
#'   `nrow(y)`. For \eqn{p} categories there are \eqn{p(p-1)/2} of them.
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
  yp <- sweep(y, 2L, p$prob, "/")
  yp2 <- sweep(y, 2L, p$prob^2, "/")
  nm <- hess_names(distrib@params)
  pr <- hess_pairs(distrib@params)
  pos <- stats::setNames(seq_along(distrib@params), distrib@params)

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
#' \deqn{\mathbb{E}[\ell^{(\eta_k\eta_l)}]
#'       = -n\sum_{j=1}^{p} \dfrac{A_{jk}A_{jl}}{p_j},}
#' with \eqn{A = \partial p/\partial\eta}. `approx` and `nsim` are ignored:
#' every strategy returns the same matrix.
#'
#' @details
#' # Why the second-derivative term vanishes
#'
#' Under the model \eqn{\mathbb{E}[y_j] = n p_j}, so the term of the observed
#' Hessian carrying \eqn{B} becomes \eqn{n\sum_j B_{j,kl}}. The probabilities
#' sum to one at every free vector, so every derivative of that sum is zero and
#' the whole term drops out, leaving the display above. The same argument, one
#' order lower, gives the score mean zero.
#'
#' # The information is the covariance of the first p-1 counts
#'
#' On the default chart, [parameters7::simplex()]'s additive log-ratio
#' coordinates, the Jacobian is \eqn{A_{jk} = p_j(\delta_{jk} - p_k)}, so
#' \deqn{\sum_j \dfrac{A_{jk}A_{jl}}{p_j}
#'   = \sum_j p_j(\delta_{jk}-p_k)(\delta_{jl}-p_l)
#'   = \delta_{kl}p_k - p_k p_l,}
#' and the expected information \eqn{-\mathbb{E}[\ell^{(\eta_k\eta_l)}]} is
#' exactly \eqn{\operatorname{Cov}(Y_k, Y_l)}, the leading
#' \eqn{(p-1)\times(p-1)} block of [mv_sigma.MultinomialDistrib()]. Measured,
#' the two agree to `2e-16` at every dimension and trial count tried. The
#' identity belongs to that chart; a different simplex parametrization has a
#' different Jacobian and no such coincidence.
#'
#' @param distrib A `MultinomialDistrib` object, from [multinomial_distrib()].
#' @param y A numeric matrix of observations. Only its row count is read
#'   through [n_obs()], the expectation not depending on the data; the values
#'   are ignored.
#' @param theta A named list of the simplex's free values on the parameter
#'   scale, each of length 1.
#' @param scale One of `"parameter"` (the default) or `"link"`, matched by
#'   [base::match.arg()]. Read by the generic, not by this method.
#' @param approx Ignored here, the expectation being exact. Accepted so that
#'   the signature matches the generic's, where it selects between
#'   `"bartlett"`, `"integrate"`, `"mc"` and `"opg"`.
#' @param nsim Ignored here, for the same reason. Defaults to `10000`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A named list of numeric vectors, one per distinct entry of the
#'   symmetric matrix and named as [hess_names()] names them, each of length
#'   `n_obs(distrib, y)` and constant along it.
#'
#' @section Notation:
#' \eqn{\ell} is the log-mass of one observation, \eqn{p} the probability
#' vector, \eqn{n} the trial count, \eqn{\eta} the free vector of the simplex
#' chart and \eqn{A = \partial p/\partial\eta} its Jacobian.
#'
#' @seealso [distrib_hessian.MultinomialDistrib()] for the quantity this is the
#'   expectation of, [mv_sigma.MultinomialDistrib()] for the covariance it
#'   coincides with, [mv_support.MultinomialDistrib()] for the exact sum it can
#'   be checked against, and [distrib_expected_hessian()] for the generic.
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
#' # The strategy argument is inert, the expectation being exact.
#' identical(eh, distrib_expected_hessian(d, matrix(0, 1, 3), th,
#'                                        approx = "mc", nsim = 50))
S7::method(distrib_expected_hessian, MultinomialDistrib) <- function(distrib, y, theta,
                                                                      scale = c("parameter", "link"),
                                                                      approx = c("bartlett", "integrate", "mc", "opg"),
                                                                      nsim = 10000, ...) {
  p <- mn_parts(distrib, theta)
  n <- n_obs(distrib, y)
  nm <- hess_names(distrib@params)
  pr <- hess_pairs(distrib@params)
  pos <- stats::setNames(seq_along(distrib@params), distrib@params)

  stats::setNames(lapply(seq_along(nm), function(m) {
    a <- pos[[pr[[m]][1]]]
    b <- pos[[pr[[m]][2]]]
    rep(-distrib@size * sum(p$A[, a] * p$A[, b] / p$prob), n)
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
#'   scale, each of length 1.
#'
#' @return A numeric vector of length \eqn{p}, strictly positive and summing to
#'   the object's `size`.
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
#' every covariance is negative: the counts share a fixed total, so one
#' category can only grow at the others' expense.
#'
#' The matrix is **singular by construction**, the coordinates summing to
#' \eqn{n}, so the vector of ones is in its null space and the rank is
#' \eqn{p-1}. Its leading \eqn{(p-1)\times(p-1)} block is the expected
#' information on the default chart; see
#' [distrib_expected_hessian.MultinomialDistrib()].
#'
#' @param distrib A `MultinomialDistrib` object, from [multinomial_distrib()].
#' @param theta A named list of the simplex's free values on the parameter
#'   scale, each of length 1.
#'
#' @return A symmetric \eqn{p \times p} numeric matrix of rank \eqn{p-1}.
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
#' # Singular by construction: rank p - 1, with the ones vector in the null
#' # space.
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
#'   scale, each of length 1.
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
#' [mv_sigma.MultinomialDistrib()]. It is singular, the counts summing to the
#' trial count, so a caller wanting something to invert should take the
#' marginals or work on the simplex's free vector.
#'
#' @param x A `MultinomialDistrib` object, from [multinomial_distrib()]. The
#'   argument is named `x` for consistency with [base::mean()], whose signature
#'   the moment generics follow.
#' @param theta A named list of the simplex's free values on the parameter
#'   scale, each of length 1.
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
#' Several coordinates at once are refused. A sub-vector is again multinomial,
#' but only after the remaining outcomes are collapsed into a category of their
#' own, and returning that object under this name would mislead. The error
#' says so.
#'
#' @param distrib A `MultinomialDistrib` object, from [multinomial_distrib()].
#' @param theta A named list of the simplex's free values on the parameter
#'   scale, each of length 1.
#' @param which A single integer in \eqn{1, \dots, p}, the category wanted. A
#'   vector of length other than 1 signals an error explaining why.
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
#' derivatives of the log-mass to fourth order and a closed-form expected
#' information, and its support is small enough to enumerate, so every
#' expectation can be checked as an exact sum.
#'
#' @param n_dim The number of categories \eqn{p}, a single integer of at least
#'   2. Anything else signals an error naming the argument.
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
#' There is **no dispersion parameter**: the probabilities are all the family
#' has, and \eqn{p-1} of them are free.
#'
#' The probabilities are carried by a `parameters7` simplex and flattened into
#' scalars with identity links, exactly as a covariance is for the multivariate
#' gaussian. The constraint that they be positive and sum to one lives in the
#' parameter, where a scalar link could not express it.
#'
#' # A multivariate family that is discrete
#'
#' This is the first family of the package that is multivariate and discrete,
#' so its support is a finite set of points rather than a region.
#' [mv_support.MultinomialDistrib()] enumerates them, which lets an expectation
#' be an exact sum and lets [check_distrib()] test the total mass by addition:
#' measured, the mass over the support comes back at 1 to `8e-16` and the
#' closed-form information agrees with the summed observed Hessian to the same
#' order. An importance-sampling check would only ever compare against Monte
#' Carlo error, and a normalization wrong by a thousandth would pass it.
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
#' the covariance of the first \eqn{p-1} counts, since
#' \eqn{A_{jk} = p_j(\delta_{jk}-p_k)} collapses the sum to
#' \eqn{\delta_{kl}p_k - p_kp_l}. Measured, the two agree to `2e-16` at every
#' dimension and trial count tried.
#'
#' # Marginals and conjugacy
#'
#' Coordinate \eqn{j} is \eqn{\mathrm{Binomial}(n, p_j)}, the other categories
#' collapsing into a single failure, so
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
#'   given here, `distrib_name` `"multinomial [pd, size=n, probs=<chart>]"`,
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
#' # Two free probabilities for three categories; no dispersion parameter.
#' th <- list(probs_alr1 = 0.3, probs_alr2 = -0.2)
#' mv_location(d, th)
#'
#' # The support is a finite set of points, so the mass sums exactly.
#' supp <- mv_support(d, th)
#' c(points = nrow(supp), mass = sum(distrib_pdf(d, supp, th)))
#'
#' # The covariance is singular, its diagonal the binomial variances.
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
#' [parameters7::simplex()] for the chart the probabilities ride;
#' [numericals7::compositions()] for the support enumeration;
#' [fit_distrib()] to estimate the probabilities; [check_distrib()] to validate
#' a family of your own against the same battery this one passes;
#' [MultinomialDistrib] for the class.
#' @export
multinomial_distrib <- function(n_dim, size,
                                probs = parameters7::simplex(n_dim)) {
  p <- as.integer(n_dim)
  if (length(p) != 1L || is.na(p) || p < 2L) {
    stop("'n_dim' must be a single integer of at least 2.", call. = FALSE)
  }
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
