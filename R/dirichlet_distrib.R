#' @include distrib.R generics.R multivariate.R moments.R
NULL

#' @title Dirichlet Distribution Class, Mean Vector and Concentration
#' @name DirichletDistrib
#'
#' @description
#' The S7 class of the Dirichlet family on the open simplex of \eqn{p}
#' coordinates, parametrized by a mean vector \eqn{\mu} carried on a
#' `parameters7` simplex and a concentration \eqn{\phi > 0}. It inherits from
#' `multivariate_distrib`, so the distribution function and the quantile are
#' refused rather than approximated; the eleven methods listed below are
#' registered on it in this file and two more in `mv_higher.R`.
#'
#' The class carries an extra property beyond the parent's, `param`: the
#' `parameters7` simplex the mean lies on, whose free names become the
#' family's own parameter names prefixed by `mean_`. Build one with
#' [dirichlet_distrib()], which validates the simplex against `n_dim`,
#' supplies the concentration's link and fills the properties in. This page
#' documents the raw S7 constructor, which validates none of that.
#'
#' @inheritParams multivariate_distrib
#' @param param The `parameters7` [parameters7::simplex()] carrying the mean.
#'   Its `n_free` is \eqn{p-1}, one fewer than the number of coordinates, the
#'   simplex being a set of dimension \eqn{p-1} in \eqn{\mathbb{R}^p}.
#'
#' @return An S7 object of class `DirichletDistrib`, inheriting from
#'   `multivariate_distrib` and from `distrib`. Beyond `param` its properties
#'   are the parent's: `distrib_name`, `dimension`, `n_dim`, `bounds`,
#'   `params`, `params_interpretation`, `n_params`, `params_bounds`,
#'   `link_params` and `params_smooth`. For an object built by
#'   [dirichlet_distrib()] `params` is `c("mean_<free names>", "phi")`,
#'   `n_params` is \eqn{p}, and every mean coordinate carries an identity link,
#'   its free value being unconstrained already.
#'
#' @seealso [dirichlet_distrib()] to build one;
#'   [beta1_distrib()], which is a coordinate's marginal and the
#'   two-coordinate case seen on the line;
#'   [multinomial_distrib()] for the discrete family on the same simplex;
#'   [parameters7::simplex()] for the chart.
#'
#' @section Methods:
#' Registered on this class in this file:
#'   [`distrib_expected_hessian()`][distrib_expected_hessian.DirichletDistrib],
#'   [`distrib_gradient()`][distrib_gradient.DirichletDistrib],
#'   [`distrib_hessian()`][distrib_hessian.DirichletDistrib],
#'   [`distrib_pdf()`][distrib_pdf.DirichletDistrib],
#'   [`distrib_rng()`][distrib_rng.DirichletDistrib],
#'   [`mv_location()`][mv_location.DirichletDistrib],
#'   [`mv_marginal()`][mv_marginal.DirichletDistrib],
#'   [`mv_reference_draw()`][mv_reference_draw.DirichletDistrib],
#'   [`mv_sigma()`][mv_sigma.DirichletDistrib],
#'   [`mean()`][mean.DirichletDistrib] and
#'   [`variance()`][variance.DirichletDistrib].
#'
#' Two more are registered in `mv_higher.R`:
#'   [`distrib_deriv3()`][distrib_deriv3.DirichletDistrib] and
#'   [`distrib_deriv4()`][distrib_deriv4.DirichletDistrib], both closed form.
#'
#' Everything else is inherited from [multivariate_distrib()], which refuses
#' the distribution function, the quantile and the response derivatives.
#'
#' @examples
#' d <- dirichlet_distrib(3)
#' S7::S7_inherits(d, multivariate_distrib)
#'
#' # Three coordinates, so two free mean values plus the concentration.
#' d@n_dim
#' d@params
#' d@params_interpretation
#'
#' # The mean lives on a parameters7 simplex, which is the extra property.
#' d@param@param_name
#' d@param@n_free
DirichletDistrib <- S7::new_class("DirichletDistrib",
  parent = multivariate_distrib,
  properties = list(param = S7::class_any)
)

#' The Pieces a Dirichlet Derivative Needs
#'
#' @description
#' Evaluates the mean vector, the concentration, the shapes
#' \eqn{\alpha = \phi\mu} and the simplex's first two derivative arrays once,
#' so that a density or a derivative method computes them a single time and
#' shares them.
#'
#' @details
#' Two identities keep every formula on the family's pages short, and both
#' follow from differentiating \eqn{\sum_j \mu_j = 1}: once, the columns of
#' \eqn{A = \partial\mu/\partial\eta} sum to zero; twice, so does every
#' second-derivative vector \eqn{B_{\cdot,kl}}. They are what makes the
#' expected information closed form, since every term carrying the data is a
#' constant times one of those zero sums and drops out under expectation.
#'
#' @param distrib A [DirichletDistrib()] object, read for its `param` and its
#'   parameter names.
#' @param theta A named list of parameters on the parameter scale: the mean's
#'   free values, in the simplex's own order, followed by the concentration
#'   `phi`.
#'
#' @return A named list with `mu`, the mean vector of length \eqn{p} summing to
#'   one; `phi`, the concentration; `alpha`, the shapes \eqn{\phi\mu}; `A`, the
#'   \eqn{p \times (p-1)} matrix of first derivatives of the mean in the free
#'   values; `B`, the list of second-derivative vectors keyed by unordered
#'   tuple; `idx`, those tuples; and `k`, the number of free mean values.
#'
#' @seealso [dir_b_index()] for locating an entry of `B`,
#'   [distrib_gradient.DirichletDistrib()] for the first consumer, and
#'   [dirichlet_distrib()] for the family.
#'
#' @keywords internal
dir_parts <- function(distrib, theta) {
  s <- distrib@param
  v <- mv_flat_theta(distrib, theta)
  k <- s@n_free
  eta <- v[seq_len(k)]
  phi <- v[[k + 1L]]
  mu <- as.numeric(parameters7::param_value(s, eta))
  list(mu = mu, phi = phi, alpha = phi * mu,
       A = do.call(cbind, parameters7::param_d1(s, eta)),
       B = parameters7::param_d2(s, eta),
       idx = parameters7::param_tuple_indices(s, 2L),
       k = k)
}

#' The Position of a Second-Derivative Block
#'
#' @description
#' Returns the position of the \eqn{(k, l)} entry in the list of second
#' derivatives a `parameters7` parameter supplies. That list is keyed by
#' unordered tuple rather than by position, so a caller holding two free-value
#' indices has to search for the entry naming them, in either order.
#'
#' @param idx The tuple index list, as returned by
#'   [parameters7::param_tuple_indices()] at order 2.
#' @param k,l The two free-value positions, in either order.
#'
#' @return A single integer, the position into the second-derivative list. `NA`
#'   if no tuple matches, which cannot happen for indices inside the
#'   parameter's own range.
#'
#' @seealso [dir_parts()], which supplies both the list and the tuples, and
#'   [dirichlet_distrib()] for the family.
#'
#' @keywords internal
dir_b_index <- function(idx, k, l) {
  which(vapply(idx, function(t) {
    length(t) == 2L && identical(sort(as.integer(t)), sort(as.integer(c(k, l))))
  }, logical(1)))[1L]
}

# --- S7 METHODS IMPLEMENTATION ---

#' @title Dirichlet Density
#' @name distrib_pdf.DirichletDistrib
#' @description
#' Computes the Dirichlet density
#' \deqn{f(y) = \dfrac{\Gamma(\phi)}{\prod_{j=1}^{p} \Gamma(\alpha_j)}
#'       \prod_{j=1}^{p} y_j^{\alpha_j - 1}, \qquad \alpha = \phi\mu,}
#' one value per row of `y`, with respect to the \eqn{(p-1)}-dimensional
#' measure the simplex carries.
#'
#' The support is tested **before** the logarithm is taken. A negative
#' coordinate would make [base::log()] warn, and a density evaluated off its
#' support is 0, so a row that is not on the simplex returns 0 without a
#' numerical complaint.
#'
#' @param distrib A `DirichletDistrib` object, from [dirichlet_distrib()].
#' @param y A numeric matrix with one row per observation and \eqn{p} columns,
#'   each row strictly positive and summing to one within `1e-8`. A single
#'   observation may be given as a plain vector and is read as one row. A row
#'   failing either test gives a density of 0, or `-Inf` with `log = TRUE`.
#' @param theta A named list of parameters on the parameter scale: the mean's
#'   free values followed by `phi`, each of length 1. A parameter may not vary
#'   by observation here, the mean being one point of the simplex for the whole
#'   sample.
#' @param log Logical of length 1. When `TRUE` the log-density is returned.
#'   Defaults to `FALSE`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of densities of length `nrow(y)`, one per
#'   observation.
#'
#' @section Notation:
#' \eqn{\mu} is the mean vector, a point of the simplex; \eqn{\phi > 0} the
#' concentration; \eqn{\alpha = \phi\mu} the shapes; and \eqn{p} the number of
#' coordinates.
#'
#' @seealso [distrib_gradient.DirichletDistrib()] for the derivatives of the
#'   log-density, [distrib_rng.DirichletDistrib()] for draws,
#'   [mv_marginal.DirichletDistrib()] for a coordinate's beta marginal, and
#'   [distrib_pdf()] for the generic.
#'
#' @examples
#' d <- dirichlet_distrib(3)
#' th <- list(mean_alr1 = 0.3, mean_alr2 = -0.2, phi = 12)
#' mu <- mv_location(d, th)
#'
#' # The density at four draws, against the formula written out.
#' set.seed(1)
#' Y <- distrib_rng(d, 4, th)
#' al <- 12 * mu
#' all.equal(distrib_pdf(d, Y, th),
#'           as.numeric(exp(lgamma(12) - sum(lgamma(al))) *
#'                        apply(Y, 1, function(r) prod(r^(al - 1)))))
#'
#' # A row that does not sum to one, or has a negative coordinate, is off the
#' # support and returns zero rather than warning.
#' distrib_pdf(d, rbind(c(0.5, 0.5, 0.5), c(-0.1, 0.6, 0.5)), th)
#'
#' # A single observation may be given as a vector.
#' distrib_pdf(d, c(0.4, 0.25, 0.35), th)
S7::method(distrib_pdf, DirichletDistrib) <- function(distrib, y, theta, log = FALSE, ...) {
  p <- dir_parts(distrib, theta)
  y <- if (is.matrix(y)) y else matrix(y, nrow = 1L)
  const <- lgamma(p$phi) - sum(lgamma(p$alpha))
  # The support is tested before the logarithm rather than after: log() of a
  # negative coordinate warns, and a density evaluated off its support is zero
  # rather than a numerical complaint.
  ok <- apply(y, 1L, function(r) all(r > 0) && abs(sum(r) - 1) <= 1e-8)
  out <- rep(-Inf, nrow(y))
  if (any(ok)) {
    out[ok] <- const + as.numeric(log(y[ok, , drop = FALSE]) %*% (p$alpha - 1))
  }
  if (log) out else exp(out)
}

#' @title Dirichlet Random Generation
#' @name distrib_rng.DirichletDistrib
#' @description
#' Draws `n` independent Dirichlet vectors by the representation the family is
#' defined by: independent gamma variates with shapes \eqn{\alpha_j = \phi\mu_j}
#' and a common rate, divided by their sum. The normalization removes the rate,
#' so [stats::rgamma()]'s default of 1 is used, and every row of the result
#' sums to one by construction. The draws depend on `.Random.seed` in the usual
#' way and consume \eqn{p} of R's streams per observation.
#'
#' @param distrib A `DirichletDistrib` object, from [dirichlet_distrib()].
#' @param n A single positive integer, the number of draws.
#' @param theta A named list of parameters on the parameter scale: the mean's
#'   free values followed by `phi`, each of length 1.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric matrix with `n` rows and \eqn{p} columns, each row
#'   strictly positive and summing to one.
#'
#' @seealso [distrib_pdf.DirichletDistrib()] for the density,
#'   [mv_reference_draw.DirichletDistrib()] for the uniform proposal
#'   [check_distrib()] integrates against, [fit_distrib()] to estimate the
#'   parameters back from a sample, and [distrib_rng()] for the generic.
#'
#' @examples
#' d <- dirichlet_distrib(3)
#' th <- list(mean_alr1 = 0.3, mean_alr2 = -0.2, phi = 12)
#'
#' # Every row is on the simplex.
#' set.seed(1)
#' Y <- distrib_rng(d, 4, th)
#' round(Y, 4)
#' rowSums(Y)
#'
#' # The sample recovers the mean and the coordinate variances, which fall as
#' # 1 / (phi + 1).
#' set.seed(3)
#' Z <- distrib_rng(d, 3e5, th)
#' rbind(sample = c(mean(Z[, 1]), var(Z[, 1])),
#'       theoretical = c(mv_location(d, th)[1], mv_sigma(d, th)[1, 1]))
S7::method(distrib_rng, DirichletDistrib) <- function(distrib, n, theta, ...) {
  p <- dir_parts(distrib, theta)
  d <- length(p$alpha)
  g <- matrix(stats::rgamma(n * d, rep(p$alpha, each = n)), nrow = n)
  g / rowSums(g)
}

#' @title Dirichlet Score
#' @name distrib_gradient.DirichletDistrib
#' @description
#' Computes the first derivatives of the Dirichlet log-density with respect to
#' the mean's free values and the concentration, one value per observation, in
#' closed form. Writing \eqn{g_j = \log y_j - \psi(\alpha_j)} and
#' \eqn{A = \partial\mu/\partial\eta} for the simplex's own Jacobian,
#' \deqn{\dfrac{\partial\ell}{\partial\eta_k} = \phi \sum_{j} g_j A_{jk},
#'       \qquad
#'       \dfrac{\partial\ell}{\partial\phi} = \psi(\phi) + \sum_{j} g_j \mu_j.}
#' The data enter only through \eqn{\log y}, the family being an exponential
#' family in that statistic, so \eqn{g} carries the whole of it.
#'
#' With `scale = "link"` the generic applies the chain rule for the links the
#' family carries before returning. Every mean coordinate rides the identity,
#' its free value being unconstrained already, so only the concentration's
#' component changes.
#'
#' @param distrib A `DirichletDistrib` object, from [dirichlet_distrib()].
#' @param y A numeric matrix with one row per observation and \eqn{p} columns,
#'   each row on the open simplex. A single observation may be given as a plain
#'   vector.
#' @param theta A named list of parameters on the parameter scale: the mean's
#'   free values followed by `phi`, each of length 1.
#' @param scale One of `"parameter"` (the default) or `"link"`, matched by
#'   [base::match.arg()]. Read by the generic, not by this method.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A named list with one numeric vector per parameter, in the order
#'   `distrib@params` gives, each of length `nrow(y)`.
#'
#' @section Notation:
#' \eqn{\ell} is the log-density of one observation, \eqn{\mu} the mean vector,
#' \eqn{\phi > 0} the concentration, \eqn{\alpha = \phi\mu} the shapes,
#' \eqn{\eta} the free vector of the simplex chart and \eqn{\psi} the digamma
#' function.
#'
#' @seealso [distrib_hessian.DirichletDistrib()] for the second derivatives,
#'   [distrib_expected_hessian.DirichletDistrib()] for their expectation,
#'   [dir_parts()] for the shared pieces, and [distrib_gradient()] for the
#'   generic.
#'
#' @examples
#' d <- dirichlet_distrib(3)
#' th <- list(mean_alr1 = 0.3, mean_alr2 = -0.2, phi = 12)
#' set.seed(1)
#' Y <- distrib_rng(d, 4, th)
#' g <- distrib_gradient(d, Y, th)
#' names(g)
#'
#' # It is the derivative of the log-density, so numDeriv on the summed
#' # log-density reproduces the summed score.
#' fn <- function(v)
#'   sum(distrib_pdf(d, Y, as.list(stats::setNames(v, d@params)), log = TRUE))
#' rbind(numeric = numDeriv::grad(fn, unlist(th)),
#'       closed = vapply(g, sum, numeric(1)))
#'
#' # The summed score vanishes at the maximum likelihood estimate.
#' set.seed(9)
#' Z <- distrib_rng(d, 800, th)
#' mle <- as.list(coef(fit_distrib(d, Z)))
#' vapply(distrib_gradient(d, Z, mle), sum, numeric(1))
S7::method(distrib_gradient, DirichletDistrib) <- function(distrib, y, theta,
                                                            scale = c("parameter", "link"), ...) {
  p <- dir_parts(distrib, theta)
  y <- if (is.matrix(y)) y else matrix(y, nrow = 1L)
  g <- sweep(t(log(y)), 1L, digamma(p$alpha), "-")
  out <- vector("list", distrib@n_params)
  for (k in seq_len(p$k)) out[[k]] <- p$phi * as.numeric(crossprod(p$A[, k], g))
  out[[p$k + 1L]] <- digamma(p$phi) + as.numeric(crossprod(p$mu, g))
  stats::setNames(out, distrib@params)
}

#' @title Dirichlet Observed Hessian
#' @name distrib_hessian.DirichletDistrib
#' @description
#' Computes the distinct second derivatives of the Dirichlet log-density with
#' respect to the mean's free values and the concentration, one value per
#' observation, in closed form. With \eqn{g_j = \log y_j - \psi(\alpha_j)},
#' \eqn{t_j = \psi'(\alpha_j)}, \eqn{A = \partial\mu/\partial\eta} and
#' \eqn{B_{\cdot,kl} = \partial^2\mu/\partial\eta_k\partial\eta_l},
#' \deqn{\ell^{(\eta_k\eta_l)} = \phi\sum_{j}\left(-t_j\phi A_{jk}A_{jl}
#'         + g_j B_{j,kl}\right), \qquad
#'       \ell^{(\eta_k\phi)} = -\phi\sum_{j} t_j \mu_j A_{jk}
#'         + \sum_{j} g_j A_{jk},}
#' \deqn{\ell^{(\phi\phi)} = \psi'(\phi) - \sum_{j} t_j \mu_j^2.}
#' The pure-concentration entry is **free of the data**, the family being an
#' exponential family in \eqn{\log y} whose concentration multiplies a
#' statistic linearly; it therefore equals its own expectation at every
#' observation. The other two carry the data through \eqn{g}.
#'
#' @param distrib A `DirichletDistrib` object, from [dirichlet_distrib()].
#' @param y A numeric matrix with one row per observation and \eqn{p} columns,
#'   each row on the open simplex. A single observation may be given as a plain
#'   vector.
#' @param theta A named list of parameters on the parameter scale: the mean's
#'   free values followed by `phi`, each of length 1.
#' @param scale One of `"parameter"` (the default) or `"link"`, matched by
#'   [base::match.arg()]. Read by the generic, not by this method.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A named list of numeric vectors, one per distinct entry of the
#'   symmetric matrix and named as [hess_names()] names them, each of length
#'   `nrow(y)`. For \eqn{p} coordinates there are \eqn{p(p+1)/2} of them.
#'
#' @section Notation:
#' \eqn{\ell} is the log-density of one observation, \eqn{\mu} the mean vector,
#' \eqn{\phi > 0} the concentration, \eqn{\alpha = \phi\mu} the shapes,
#' \eqn{\eta} the free vector of the simplex chart, and \eqn{\psi}, \eqn{\psi'}
#' the digamma and trigamma functions.
#'
#' @seealso [distrib_gradient.DirichletDistrib()] for the score,
#'   [distrib_expected_hessian.DirichletDistrib()] for the expectation of this
#'   quantity, [distrib_deriv3.DirichletDistrib()] for the order above, and
#'   [distrib_hessian()] for the generic.
#'
#' @examples
#' d <- dirichlet_distrib(3)
#' th <- list(mean_alr1 = 0.3, mean_alr2 = -0.2, phi = 12)
#' set.seed(1)
#' Y <- distrib_rng(d, 4, th)
#' h <- distrib_hessian(d, Y, th)
#' names(h)
#'
#' # The pure-concentration entry is the same number at every observation.
#' h$phi_phi
#'
#' # numDeriv on the summed log-density reproduces the summed matrix.
#' fn <- function(v)
#'   sum(distrib_pdf(d, Y, as.list(stats::setNames(v, d@params)), log = TRUE))
#' H <- numDeriv::hessian(fn, unlist(th))
#' rbind(numeric = c(H[1, 1], H[2, 2], H[3, 3], H[1, 2], H[1, 3], H[2, 3]),
#'       closed = vapply(h, sum, numeric(1)))
S7::method(distrib_hessian, DirichletDistrib) <- function(distrib, y, theta,
                                                           scale = c("parameter", "link"), ...) {
  p <- dir_parts(distrib, theta)
  y <- if (is.matrix(y)) y else matrix(y, nrow = 1L)
  n <- nrow(y)
  g <- sweep(t(log(y)), 1L, digamma(p$alpha), "-")
  tt <- trigamma(p$alpha)
  kk <- p$k
  nm <- hess_names(distrib@params)
  pr <- hess_pairs(distrib@params)
  ix <- match(distrib@params, distrib@params)
  pos <- stats::setNames(seq_along(distrib@params), distrib@params)

  stats::setNames(lapply(seq_along(nm), function(m) {
    a <- pos[[pr[[m]][1]]]
    b <- pos[[pr[[m]][2]]]
    if (a <= kk && b <= kk) {
      bi <- dir_b_index(p$idx, a, b)
      base <- -p$phi * p$phi * sum(tt * p$A[, a] * p$A[, b])
      rep(base, n) + p$phi * as.numeric(crossprod(p$B[[bi]], g))
    } else if (a > kk && b > kk) {
      rep(trigamma(p$phi) - sum(tt * p$mu^2), n)
    } else {
      j <- if (a > kk) b else a
      base <- -p$phi * sum(tt * p$mu * p$A[, j])
      rep(base, n) + as.numeric(crossprod(p$A[, j], g))
    }
  }), nm)
}

#' @title Dirichlet Expected Hessian
#' @name distrib_expected_hessian.DirichletDistrib
#' @description
#' Returns the expectation of the observed Hessian under the model, in closed
#' form and with no quadrature or simulation:
#' \deqn{\mathbb{E}[\ell^{(\eta_k\eta_l)}] = -\phi^2\sum_{j} t_j A_{jk}A_{jl},
#'       \qquad
#'       \mathbb{E}[\ell^{(\eta_k\phi)}] = -\phi\sum_{j} t_j \mu_j A_{jk},
#'       \qquad
#'       \mathbb{E}[\ell^{(\phi\phi)}] = \psi'(\phi) - \sum_{j} t_j \mu_j^2,}
#' with \eqn{t_j = \psi'(\alpha_j)} and \eqn{A = \partial\mu/\partial\eta}.
#' `approx` and `nsim` are ignored: every strategy returns the same matrix.
#'
#' @details
#' # Why the data-carrying terms drop out
#'
#' Under the model \eqn{\mathbb{E}[\log y_j] = \psi(\alpha_j) - \psi(\phi)}, so
#' \eqn{\mathbb{E}[g_j] = -\psi(\phi)}, **the same constant for every**
#' \eqn{j}. Each term of the observed Hessian that carries the data is that
#' constant times a sum over \eqn{j} of a column of \eqn{A} or of a
#' second-derivative vector of the simplex, and both of those sum to zero:
#' differentiating \eqn{\sum_j \mu_j = 1} once gives the first, twice the
#' second. Every such term is therefore exactly zero under expectation, which
#' is why a family with no location and scale to separate still has a
#' closed-form information matrix.
#'
#' @param distrib A `DirichletDistrib` object, from [dirichlet_distrib()].
#' @param y A numeric matrix of observations. Only its row count is read
#'   through [n_obs()], the expectation not depending on the data; the values
#'   are ignored.
#' @param theta A named list of parameters on the parameter scale: the mean's
#'   free values followed by `phi`, each of length 1.
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
#' \eqn{\ell} is the log-density of one observation, \eqn{\mu} the mean vector,
#' \eqn{\phi > 0} the concentration, \eqn{\alpha = \phi\mu} the shapes,
#' \eqn{\eta} the free vector of the simplex chart, and \eqn{\psi}, \eqn{\psi'}
#' the digamma and trigamma functions.
#'
#' @seealso [distrib_hessian.DirichletDistrib()] for the quantity this is the
#'   expectation of, [distrib_gradient.DirichletDistrib()] for the score, and
#'   [distrib_expected_hessian()] for the generic.
#'
#' @examples
#' d <- dirichlet_distrib(3)
#' th <- list(mean_alr1 = 0.3, mean_alr2 = -0.2, phi = 12)
#' set.seed(1)
#' Y <- distrib_rng(d, 4, th)
#' eh <- distrib_expected_hessian(d, Y, th)
#' vapply(eh, function(v) v[1], numeric(1))
#'
#' # Averaging the observed Hessian over draws reaches the same matrix.
#' set.seed(3)
#' Z <- distrib_rng(d, 3e5, th)
#' vapply(distrib_hessian(d, Z, th), mean, numeric(1))
#'
#' # The two zero sums the derivation rests on, read off the simplex itself.
#' eta <- c(0.3, -0.2)
#' colSums(do.call(cbind, parameters7::param_d1(d@param, eta)))
#' vapply(parameters7::param_d2(d@param, eta), sum, numeric(1))
#'
#' # The strategy argument is inert, the expectation being exact.
#' identical(eh, distrib_expected_hessian(d, Y, th, approx = "mc", nsim = 50))
S7::method(distrib_expected_hessian, DirichletDistrib) <- function(distrib, y, theta,
                                                                    scale = c("parameter", "link"),
                                                                    approx = c("bartlett", "integrate", "mc", "opg"),
                                                                    nsim = 10000, ...) {
  p <- dir_parts(distrib, theta)
  n <- n_obs(distrib, y)
  tt <- trigamma(p$alpha)
  kk <- p$k
  nm <- hess_names(distrib@params)
  pr <- hess_pairs(distrib@params)
  pos <- stats::setNames(seq_along(distrib@params), distrib@params)

  stats::setNames(lapply(seq_along(nm), function(m) {
    a <- pos[[pr[[m]][1]]]
    b <- pos[[pr[[m]][2]]]
    val <- if (a <= kk && b <= kk) {
      -p$phi * p$phi * sum(tt * p$A[, a] * p$A[, b])
    } else if (a > kk && b > kk) {
      trigamma(p$phi) - sum(tt * p$mu^2)
    } else {
      j <- if (a > kk) b else a
      -p$phi * sum(tt * p$mu * p$A[, j])
    }
    rep(val, n)
  }), nm)
}

#' @title Dirichlet Mean Vector
#' @name mv_location.DirichletDistrib
#' @description
#' Returns the point of the simplex the mean parameter carries at the free
#' values in `theta`. For this family that point **is** the mean,
#' \eqn{\mathbb{E}[Y_j] = \mu_j}, so [mean.DirichletDistrib()] delegates here.
#' The value comes from [parameters7::param_value()] on the object's own
#' simplex chart, so it lies on the open simplex by construction whatever the
#' free values are.
#'
#' @param distrib A `DirichletDistrib` object, from [dirichlet_distrib()].
#' @param theta A named list of parameters on the parameter scale: the mean's
#'   free values followed by `phi`, each of length 1. The concentration is read
#'   but not used here.
#'
#' @return A numeric vector of length \eqn{p}, strictly positive and summing to
#'   one.
#'
#' @seealso [mv_sigma.DirichletDistrib()] for the covariance,
#'   [mean.DirichletDistrib()], which calls this, and [mv_location()] for the
#'   generic.
#'
#' @examples
#' d <- dirichlet_distrib(3)
#' mu <- mv_location(d, list(mean_alr1 = 0.3, mean_alr2 = -0.2, phi = 12))
#' mu
#' sum(mu)
#'
#' # Whatever the free values, the result is on the simplex.
#' m2 <- mv_location(d, list(mean_alr1 = -8, mean_alr2 = 5, phi = 12))
#' c(min = min(m2), sum = sum(m2))
S7::method(mv_location, DirichletDistrib) <- function(distrib, theta) {
  dir_parts(distrib, theta)$mu
}

#' @title Dirichlet Covariance Matrix
#' @name mv_sigma.DirichletDistrib
#' @description
#' Returns the covariance matrix of the family,
#' \deqn{\operatorname{Cov}(Y_i, Y_j)
#'       = \dfrac{\delta_{ij}\mu_i - \mu_i\mu_j}{\phi + 1},}
#' so the coordinate variances are \eqn{\mu_j(1-\mu_j)/(\phi+1)} and every
#' covariance is negative, a rise in one coordinate having to be paid for by
#' the others.
#'
#' The matrix is **singular by construction**: the coordinates sum to one, so
#' the vector of ones is in its null space and the rank is \eqn{p-1}. Anything
#' that inverts a covariance must use the marginals or the free vector instead.
#' The concentration acts as a precision, every entry falling as
#' \eqn{1/(\phi+1)}.
#'
#' @param distrib A `DirichletDistrib` object, from [dirichlet_distrib()].
#' @param theta A named list of parameters on the parameter scale: the mean's
#'   free values followed by `phi`, each of length 1.
#'
#' @return A symmetric \eqn{p \times p} numeric matrix of rank \eqn{p-1}.
#'
#' @seealso [mv_location.DirichletDistrib()] for the mean,
#'   [variance.DirichletDistrib()], which calls this,
#'   [mv_marginal.DirichletDistrib()] for a coordinate's law, and
#'   [mv_sigma()] for the generic.
#'
#' @examples
#' d <- dirichlet_distrib(3)
#' th <- list(mean_alr1 = 0.3, mean_alr2 = -0.2, phi = 12)
#' S <- mv_sigma(d, th)
#' round(S, 6)
#'
#' # Singular by construction: the rank is p - 1 and the ones vector is in
#' # the null space.
#' c(rank = qr(S)$rank, dim = ncol(S))
#' round(S %*% rep(1, 3), 12)
#'
#' # The concentration is a precision: every entry falls as 1 / (phi + 1).
#' vapply(c(2, 12, 100), function(p)
#'   mv_sigma(d, list(mean_alr1 = 0.3, mean_alr2 = -0.2, phi = p))[1, 1],
#'   numeric(1))
S7::method(mv_sigma, DirichletDistrib) <- function(distrib, theta) {
  p <- dir_parts(distrib, theta)
  (diag(p$mu, nrow = length(p$mu)) - tcrossprod(p$mu)) / (p$phi + 1)
}

#' @title Mean of a Dirichlet
#' @name mean.DirichletDistrib
#' @description
#' Returns the mean vector \eqn{\mathbb{E}[Y] = \mu}, which for this family is
#' a parameter and needs no computation: the method delegates to
#' [mv_location.DirichletDistrib()]. The result sums to one, being a point of
#' the simplex.
#'
#' @param x A `DirichletDistrib` object, from [dirichlet_distrib()]. The
#'   argument is named `x` because the generic is [base::mean()].
#' @param theta A named list of parameters on the parameter scale: the mean's
#'   free values followed by `phi`, each of length 1.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of length \eqn{p}, strictly positive and summing to
#'   one.
#'
#' @seealso [variance.DirichletDistrib()] for the covariance,
#'   [mv_location.DirichletDistrib()], which this calls, and
#'   [mean.distrib()] for the generic's default.
#' @keywords internal
S7::method(mean, DirichletDistrib) <- function(x, theta, ...) {
  mv_location(x, theta)
}

#' @title Variance of a Dirichlet
#' @name variance.DirichletDistrib
#' @description
#' Returns the covariance matrix of the family by delegating to
#' [mv_sigma.DirichletDistrib()]. It is singular, the coordinates summing to
#' one, so a caller wanting something to invert should take the marginals or
#' work on the simplex's free vector.
#'
#' @param x A `DirichletDistrib` object, from [dirichlet_distrib()]. The
#'   argument is named `x` for consistency with [base::mean()], whose signature
#'   the moment generics follow.
#' @param theta A named list of parameters on the parameter scale: the mean's
#'   free values followed by `phi`, each of length 1.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A symmetric \eqn{p \times p} numeric matrix of rank \eqn{p-1}.
#'
#' @seealso [mean.DirichletDistrib()] for the mean vector,
#'   [mv_sigma.DirichletDistrib()], which this calls, and [variance()] for the
#'   generic.
#' @keywords internal
S7::method(variance, DirichletDistrib) <- function(x, theta, ...) {
  mv_sigma(x, theta)
}

#' @title A Uniform Proposal on the Simplex
#' @name mv_reference_draw.DirichletDistrib
#' @description
#' Supplies the importance-sampling proposal [check_distrib()] integrates the
#' density against: the uniform distribution on the simplex, which is the
#' Dirichlet with every shape equal to one. Its density is the constant
#' \eqn{\Gamma(p)} with respect to the same dominating measure the family's own
#' density is written against, so the normalization check is a plain average.
#'
#' The base class's proposal is an inflated gaussian on \eqn{\mathbb{R}^p},
#' which places no mass on the simplex at all. It does not fail loudly there:
#' [base::chol()] accepts the singular covariance, and the estimate of an
#' integral that is 1 comes back at about `2e-08`. Overriding the proposal is
#' what makes that check mean anything for this family.
#'
#' @param distrib A `DirichletDistrib` object, from [dirichlet_distrib()], read
#'   for its `n_dim` alone.
#' @param theta A named list of parameters. Ignored: the proposal is uniform
#'   and does not depend on where the family sits.
#' @param n A single positive integer, the number of draws.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A named list with `y`, an `n` by \eqn{p} matrix whose rows are
#'   uniform on the simplex, and `logd`, a numeric vector of length `n` holding
#'   the constant \eqn{\log\Gamma(p)}.
#'
#' @seealso [check_distrib()], which consumes this,
#'   [distrib_rng.DirichletDistrib()] for draws from the family itself, and
#'   [mv_reference_draw()] for the generic.
#' @keywords internal
S7::method(mv_reference_draw, DirichletDistrib) <- function(distrib, theta, n, ...) {
  p <- distrib@n_dim
  g <- matrix(stats::rexp(n * p), nrow = n)
  list(y = g / rowSums(g), logd = rep(lgamma(p), n))
}

#' @title Dirichlet Marginal
#' @name mv_marginal.DirichletDistrib
#' @description
#' Returns the law of one coordinate, which is
#' \eqn{\mathrm{Beta}(\alpha_j, \phi - \alpha_j)}. In this package's mean and
#' precision parametrization of the beta that is mean \eqn{\mu_j} and precision
#' \eqn{\phi}, so no reparametrization is needed at all and **the
#' concentration is shared by every marginal**, as the multivariate t's degrees
#' of freedom are.
#'
#' This is one of the few families for which the generic returns an object
#' rather than signaling an error. A panel of a pairs plot is therefore a real
#' distribution here.
#'
#' Several coordinates at once are refused. A sub-vector of a Dirichlet is
#' again Dirichlet, but only after the remaining mass is collapsed into a
#' coordinate of its own, and returning that object under this name would
#' mislead. The error says so.
#'
#' @param distrib A `DirichletDistrib` object, from [dirichlet_distrib()].
#' @param theta A named list of parameters on the parameter scale: the mean's
#'   free values followed by `phi`, each of length 1.
#' @param which A single integer in \eqn{1, \dots, p}, the coordinate wanted.
#'   A vector of length other than 1 signals an error explaining why.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A named list with `distrib`, a [beta1_distrib()] object, and
#'   `theta`, a list with the marginal's `mu` and `phi`.
#'
#' @seealso [beta1_distrib()] for the family returned,
#'   [mv_sigma.DirichletDistrib()] for the coordinate variances, and
#'   [mv_marginal()] for the generic.
#'
#' @examples
#' d <- dirichlet_distrib(3)
#' th <- list(mean_alr1 = 0.3, mean_alr2 = -0.2, phi = 12)
#' m <- mv_marginal(d, th, which = 1)
#' m$distrib
#' m$theta
#'
#' # The marginal's mean and variance are the first coordinate's.
#' rbind(marginal = c(mean(m$distrib, m$theta), variance(m$distrib, m$theta)),
#'       joint = c(mv_location(d, th)[1], mv_sigma(d, th)[1, 1]))
#'
#' # Every marginal carries the same concentration.
#' vapply(1:3, function(j) mv_marginal(d, th, which = j)$theta$phi, numeric(1))
S7::method(mv_marginal, DirichletDistrib) <- function(distrib, theta, which, ...) {
  if (length(which) != 1L) {
    stop(paste0(
      "A Dirichlet marginal is returned one coordinate at a time. Several\n",
      "  coordinates are again Dirichlet, but only after the remaining mass is\n",
      "  collapsed into a coordinate of its own, which is a different object\n",
      "  from a sub-vector and would be misleading to return under this name."
    ), call. = FALSE)
  }
  p <- dir_parts(distrib, theta)
  # The Beta here is written in a mean and a precision, exactly as this family
  # is, so the marginal needs no reparametrization at all: coordinate j has
  # shapes (alpha_j, phi - alpha_j), hence mean mu_j and precision phi. The
  # CONCENTRATION IS SHARED by every marginal, as the multivariate t's degrees
  # of freedom are.
  list(distrib = beta1_distrib(),
       theta = list(mu = p$mu[which], phi = p$phi))
}

# --- CONSTRUCTOR WRAPPER ---

#' Dirichlet Distribution, Mean Vector and Concentration
#'
#' @description
#' Builds the distribution object for the Dirichlet family on the open simplex
#' of \eqn{p} coordinates, parametrized by a mean vector \eqn{\mu} and a
#' concentration \eqn{\phi > 0}. The returned object carries closed-form
#' derivatives of the log-density to fourth order and a **closed-form expected
#' information**, which is unusual for a family with no location and scale to
#' separate.
#'
#' The mean is carried on a `parameters7` simplex and flattened into scalar
#' parameters, so every generic of the package indexes it as it does any other
#' family.
#'
#' @param n_dim The number of coordinates \eqn{p}, a single integer of at least
#'   2. Anything else signals an error naming the argument.
#' @param mean A `parameters7` parameter producing \eqn{p} coordinates that sum
#'   to one, normally a [parameters7::simplex()]. Defaults to
#'   `parameters7::simplex(n_dim)`. An object that is not a `parameters7`
#'   parameter, or that produces a different number of coordinates, signals an
#'   error naming both counts.
#' @param link_phi A `link` object from `linkfunctions7` for the concentration
#'   \eqn{\phi}. Defaults to [linkfunctions7::log_link()], which maps
#'   \eqn{(0, \infty)} onto the line and so keeps every fitted value positive.
#'   The mean's free values are unconstrained already and carry the identity.
#'
#' @details
#' # The parametrization
#'
#' Writing \eqn{\alpha = \phi\mu} for the shapes, the density on the simplex
#' \eqn{\{y_j > 0, \sum_j y_j = 1\}} is
#' \deqn{f(y) = \frac{\Gamma(\phi)}{\prod_{j=1}^{p}\Gamma(\alpha_j)}
#'       \prod_{j=1}^{p} y_j^{\alpha_j - 1},}
#' with
#' \deqn{\mathbb{E}[Y_j] = \mu_j, \qquad
#'       \operatorname{Cov}(Y_i, Y_j) = \frac{\delta_{ij}\mu_i - \mu_i\mu_j}{\phi + 1}.}
#' So \eqn{\phi} is a **precision**: the larger it is, the tighter the draws
#' sit about the mean. The covariance is singular, the coordinates summing to
#' one, and every off-diagonal entry is negative.
#'
#' The parametrization follows the design of the multivariate gaussian's. The
#' constrained object, here a point of the simplex, is carried by a
#' `parameters7` parameter and **flattened into scalars** with identity links,
#' so `align_theta()`, the derivative names, the link scale and
#' [fit_distrib()] need no special case. The parameter names are the simplex's
#' own free names prefixed by `mean_`, followed by `phi`.
#'
#' # A multivariate family that is not elliptical
#'
#' This is the first family of the package with no location and scale to
#' separate: the support is a set of dimension \eqn{p-1}, not a Euclidean
#' space, and the covariance is singular by construction. It is therefore the
#' real test of the multivariate layer, and two of its methods exist to answer
#' that test. [mv_marginal.DirichletDistrib()] returns a beta object rather
#' than an error, and [mv_reference_draw.DirichletDistrib()] replaces the base
#' class's gaussian proposal with the uniform on the simplex, without which
#' [check_distrib()] would estimate an integral that is 1 at about `2e-08`.
#'
#' # Why the expected information is closed form
#'
#' Two identities do the work. Differentiating \eqn{\sum_j \mu_j = 1} once
#' shows that the columns of \eqn{A = \partial\mu/\partial\eta} sum to zero,
#' and twice that every second-derivative vector of the simplex does; and
#' \eqn{\mathbb{E}[\log y_j] = \psi(\alpha_j) - \psi(\phi)} makes
#' \eqn{\mathbb{E}[g_j] = -\psi(\phi)} the same constant for every \eqn{j}.
#' Every term of the observed Hessian that carries the data is that constant
#' times one of the two zero sums, and vanishes.
#'
#' # What is refused
#'
#' The distribution function and the quantile are refused by
#' [multivariate_distrib()], as for every family of that class: a distribution
#' function on \eqn{\mathbb{R}^p} is an orthant probability and a quantile
#' needs an ordering. The response derivatives are refused too. A marginal over
#' several coordinates at once is refused with the reason.
#'
#' # Estimation
#'
#' [fit_distrib()] maximizes the log-likelihood on the link scale, over the
#' simplex's free values and \eqn{\log\phi}. No estimate is closed form.
#'
#' @section Notation:
#' \eqn{\ell} is the log-density of one observation, \eqn{\mu} the mean vector,
#' \eqn{\phi > 0} the concentration, \eqn{\alpha = \phi\mu} the shapes,
#' \eqn{p} the number of coordinates, \eqn{\eta} the free vector of the simplex
#' chart and \eqn{\psi} the digamma function.
#'
#' @return An S7 object of class `DirichletDistrib`, inheriting from
#'   `multivariate_distrib`, with `param` the simplex given here,
#'   `distrib_name` `"dirichlet [pd, mean=<chart>]"`, `dimension`
#'   `"multivariate"`, `n_dim` \eqn{p}, `bounds` `c(0, 1)`, `params` the
#'   simplex's free names prefixed by `mean_` followed by `"phi"`, `n_params`
#'   \eqn{p}, and `link_params` the identity for each mean coordinate and
#'   `link_phi` for the concentration.
#'
#' @references
#' Kotz, S., Balakrishnan, N. and Johnson, N. L. (2000).
#' *Continuous Multivariate Distributions*, Volume 1, 2nd edition, Chapter 49.
#' Wiley, New York.
#'
#' Aitchison, J. (1986). *The Statistical Analysis of Compositional Data*.
#' Chapman and Hall, London.
#'
#' @importFrom linkfunctions7 log_link identity_link
#' @importFrom stats rgamma
#'
#' @examples
#' d <- dirichlet_distrib(3)
#' d
#'
#' # Two free mean values and a concentration, all on the parameter scale.
#' th <- list(mean_alr1 = 0.3, mean_alr2 = -0.2, phi = 12)
#' mv_location(d, th)
#' round(mv_sigma(d, th), 6)
#'
#' # The covariance is singular: rank p - 1, with the ones vector in the null
#' # space.
#' c(rank = qr(mv_sigma(d, th))$rank, dim = 3)
#'
#' # The marginals are beta with the same concentration, so a panel of a pairs
#' # plot is a real distribution.
#' mv_marginal(d, th, which = 1)$theta
#'
#' # The density and the sample agree on the mean and the coordinate variance.
#' set.seed(3)
#' Z <- distrib_rng(d, 3e5, th)
#' rbind(sample = c(mean(Z[, 1]), var(Z[, 1])),
#'       theoretical = c(mv_location(d, th)[1], mv_sigma(d, th)[1, 1]))
#'
#' # Fitting recovers the parameters.
#' set.seed(9)
#' coef(fit_distrib(d, distrib_rng(d, 800, th)))
#'
#' @seealso
#' [beta1_distrib()] for a coordinate's marginal and the two-coordinate case
#' seen on the line; [multinomial_distrib()] for the discrete family on the
#' same simplex; [mvgaussian1_distrib()] for the elliptical multivariate
#' family; [parameters7::simplex()] for the chart the mean rides;
#' [fit_distrib()] to estimate the parameters; [check_distrib()] to validate a
#' family of your own against the same battery this one passes;
#' [DirichletDistrib] for the class.
#' @export
dirichlet_distrib <- function(n_dim, mean = parameters7::simplex(n_dim),
                              link_phi = log_link()) {
  p <- as.integer(n_dim)
  if (length(p) != 1L || is.na(p) || p < 2L) {
    stop("'n_dim' must be a single integer of at least 2.", call. = FALSE)
  }
  if (!S7::S7_inherits(mean, parameters7::parameter)) {
    stop("'mean' must be a parameters7 parameter.", call. = FALSE)
  }
  if (length(parameters7::param_value(mean, numeric(mean@n_free))) != p) {
    stop(sprintf(
      "The mean parameter produces %d coordinates but the distribution has %d.",
      length(parameters7::param_value(mean, numeric(mean@n_free))), p
    ), call. = FALSE)
  }

  params <- c(paste0("mean_", mean@free_names), "phi")
  n_par <- length(params)

  DirichletDistrib(
    distrib_name = sprintf("dirichlet [%dd, mean=%s]", p, mean@param_name),
    dimension = "multivariate", n_dim = p, bounds = c(0, 1),
    params = params,
    params_interpretation = stats::setNames(
      c(rep("mean", mean@n_free), "concentration"), params
    ),
    n_params = n_par,
    params_bounds = stats::setNames(
      c(rep(list(c(-Inf, Inf)), mean@n_free), list(c(0, Inf))), params
    ),
    link_params = stats::setNames(
      c(rep(list(linkfunctions7::identity_link()), mean@n_free),
        list(link_phi)), params
    ),
    param = mean
  )
}

