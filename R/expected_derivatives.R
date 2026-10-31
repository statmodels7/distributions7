#' @include distrib.R generics.R utility_functions.R numerical_functions.R
NULL

#' Methods for the Expectation of a Derivative
#'
#' @name expected_derivative_methods
#'
#' @description
#' An expected derivative of the log-density, such as the expected Hessian (the
#' negative of the Fisher information) or the expected third derivative, is the
#' expectation \eqn{\mathbb{E}[\,\cdot\,]} of a function of the response taken
#' under the distribution at the parameters. A family may write that
#' expectation in closed form. When it does not, or when a different route is
#' wanted, the expectation is obtained by one of three methods, which are
#' objects: [expect_opg()], [expect_integrate()] and [expect_mc()]. All of them
#' inherit from [ExpectMethod] and are passed through the argument `method` of
#' [distrib_expected_hessian()], [distrib_deriv3()], [distrib_deriv4()], the
#' derivatives of the expected information ([distrib_dexpected_hessian()] and
#' the functions of its family), [fisher_scoring()] and [expected_info()].
#'
#' The argument `method` has one rule. If it is `NULL`, which is the default,
#' the function returns the closed form of the family where the family has one,
#' and otherwise the default route of that function (the table below). If it is
#' an object, the object is always used, also for a family that has the
#' expectation in closed form; the closed form is then bypassed, and the result
#' is the expectation as the method computes it.
#'
#' @details
#' # Two identities and three ways
#'
#' Let \eqn{\ell} be the log-density and \eqn{\ell_{i}, \ell_{ij}, \dots} its
#' derivatives with respect to the parameters. The expected derivative of order
#' \eqn{k} has two expressions inside the expectation, which are equal in
#' expectation for a regular model.
#'
#' \describe{
#'   \item{The Hessian identity}{The expectation of the observed derivative of
#'     order \eqn{k} itself,
#'     \deqn{\mathbb{E}[\ell_{i_1 \dots i_k}].}
#'     At order 2 this is the expected Hessian read literally.}
#'   \item{The score identity}{The Bartlett partition identity, which writes
#'     the same expectation through products of derivatives of lower order. The
#'     identity states that
#'     \deqn{\sum_{\pi} \mathbb{E}\left[\prod_{B \in \pi} \ell_B\right] = 0,}
#'     where \eqn{\pi} runs over the partitions of the index set
#'     \eqn{\{i_1,\dots,i_k\}} and \eqn{\ell_B} is the derivative with respect
#'     to the parameters indexed by the block \eqn{B}. The single-block
#'     partition is the target, so that
#'     \deqn{\mathbb{E}[\ell_{i_1 \dots i_k}] = -\sum_{\pi,\ |\pi| \ge 2}
#'     \mathbb{E}\left[\prod_{B \in \pi} \ell_B\right].}
#'     At order 2 this is \eqn{-\mathbb{E}[\ell_i \ell_j]}, at order 3
#'     \eqn{-\mathbb{E}[\ell_{ij}\ell_k + \ell_{ik}\ell_j + \ell_{jk}\ell_i +
#'     \ell_i\ell_j\ell_k]}, and at order 4 the sum over the partitions of four
#'     indices into at least two blocks. The target order is never needed.}
#' }
#'
#' The expectation of either quantity is then obtained in one of three ways:
#' computed by quadrature (or by a series over the support of a discrete
#' family), simulated, or estimated from the sample in hand. The methods fill
#' the table as follows.
#'
#' \tabular{lccc}{
#'   \tab computed by quadrature \tab simulated \tab estimated from the sample \cr
#'   Hessian identity \tab `expect_integrate(identity = "hessian")` \tab
#'     `expect_mc(identity = "hessian")` \tab not available \cr
#'   score identity \tab `expect_integrate(identity = "score")` \tab
#'     `expect_mc(identity = "score")` \tab `expect_opg()`
#' }
#'
#' The estimate from the sample has no Hessian-identity entry because the
#' observed derivative at the data, averaged over the observations, is the
#' observed information and not an expectation. Reading the score identity at
#' the data gives the empirical outer product of the scores, which is the
#' BHHH estimator at order 2.
#'
#' \describe{
#'   \item{[expect_integrate()]}{Computes the expectation by batched adaptive
#'     quadrature for a continuous family and by series summation for a
#'     discrete one. It is deterministic, and it is normally the most accurate
#'     method. With `identity = "hessian"` it costs one integral per component
#'     and needs the derivative of the target order, which every family of the
#'     package supplies or approximates by finite differences. With
#'     `identity = "score"` it costs one integral per partition and component
#'     and needs only the lower orders, which is the reason to choose it where
#'     only those are available in closed form. The tolerances and limits are
#'     arguments of the constructor. For a family with several parameters
#'     varying by observation, one integral is computed per distinct
#'     parameter configuration.}
#'   \item{[expect_mc()]}{Draws `nsim` observations from the distribution, with
#'     [distrib_rng()], at each parameter configuration and averages the
#'     quantity of the identity over them. It estimates the same quantity as
#'     [expect_integrate()] with the same `identity`. The error decreases as
#'     \eqn{1/\sqrt{n_{sim}}}, and the result is reproducible only if the seed
#'     is fixed. The cost is that of [distrib_rng()] times `nsim`, per
#'     parameter configuration. It is the simplest method where quadrature
#'     performs poorly (heavy tails, awkward supports); for a distribution
#'     whose generator inverts a slow quantile function, [expect_integrate()]
#'     is cheaper.}
#'   \item{[expect_opg()]}{Reads the score identity at each observation and
#'     takes no expectation: nothing is integrated and nothing is simulated.
#'     At order 2 the quantity is \eqn{-\ell_i \ell_j}, evaluated with one
#'     call to [distrib_gradient()]; at orders 3 and 4 it is the sum of
#'     products of the identity, evaluated with the derivatives of orders 1 to
#'     \eqn{k-1}. Each observation contributes a noisy rank-one-like term, and
#'     the sum over \eqn{n} observations converges to \eqn{n} times the
#'     expectation at the rate \eqn{n^{-1/2}}. At order 2 the matrix is
#'     positive semidefinite by construction, which makes it adequate for the
#'     step of Fisher scoring, and it has the same fixed point as the
#'     expected information. The variance grows with the order, because the
#'     terms are products of up to \eqn{k} derivatives. The method is never
#'     routed to quadrature at any order.}
#' }
#'
#' # Defaults of each function
#'
#' When `method` is `NULL`, the function returns the closed form of the family
#' where there is one. The route used otherwise is:
#'
#' \tabular{ll}{
#'   [distrib_expected_hessian()] \tab [expect_opg()] \cr
#'   [distrib_deriv3()] and [distrib_deriv4()] with `expected = TRUE` \tab
#'     [expect_integrate()] (`identity = "hessian"`) \cr
#'   [fisher_scoring()] and [expected_info()] \tab [expect_opg()] (order 2) \cr
#'   [distrib_dexpected_hessian()], [distrib_d2expected_hessian()],
#'     [distrib_d3expected_hessian()] and [distrib_d4expected_hessian()] \tab
#'     none: a family without a closed form signals an error
#' }
#'
#' These defaults reproduce the values that the functions returned before
#' `method` was introduced.
#'
#' # Derivatives of the expected information
#'
#' The functions that return the derivatives of the expected information with
#' respect to the parameters have the closed form of the family as their only
#' default. With a method, the \eqn{m}th derivative of the matrix
#' \eqn{E_{ab}(\theta)} that the method defines is returned.
#'
#' For [expect_integrate()] and [expect_mc()] the matrix is the expectation
#' \eqn{E_{ab}(\theta) = \mathbb{E}_\theta[f_{ab}]} of the integrand \eqn{f_{ab}}
#' of the identity (\eqn{\ell_{ab}} for the Hessian identity and
#' \eqn{-\ell_a\ell_b} for the score identity), and the derivative moves the
#' density as well as the integrand,
#' \deqn{\frac{\partial}{\partial\theta_c}\mathbb{E}_\theta[f] =
#' \mathbb{E}_\theta\left[\frac{\partial f}{\partial\theta_c} + f\,\ell_c\right],}
#' repeated for higher orders. The expectation of the right-hand side is taken
#' by the method.
#'
#' For [expect_opg()] the matrix is the sum over the observations of
#' \eqn{-\ell_a\ell_b} at the data, which is the matrix whose inverse the
#' scoring step uses, and the density does not move: the exact derivative of
#' that matrix with respect to \eqn{\theta_c} is, per observation,
#' \deqn{-\left(\ell_{ac}\,\ell_b + \ell_a\,\ell_{bc}\right),}
#' and the second derivative is
#' \deqn{-\left(\ell_{acd}\,\ell_b + \ell_{ac}\,\ell_{bd} + \ell_{ad}\,\ell_{bc}
#' + \ell_a\,\ell_{bcd}\right).}
#' These are not estimates of \eqn{\partial \mathbb{E}[\ell_{ab}]/\partial
#' \theta_c}, which contains the moved density; they are the derivatives of the
#' matrix that [expect_opg()] returns, and a criterion built on that matrix
#' needs these. The derivatives of order 3 and 4 follow the same rule by
#' Leibniz's formula.
#'
#' A derivative of order \eqn{m} of the matrix needs the observed derivatives
#' of the log-density up to order \eqn{m+1} for the score identity, which is
#' always the case for \eqn{m \le 4}, and up to order \eqn{m+2} for the Hessian
#' identity, which is available for \eqn{m \le 3}; for \eqn{m = 4} the Hessian
#' identity signals an error that names `identity = "score"`.
#'
#' # A kink in the log-density
#'
#' Where the log-density has a kink in a parameter, the observed second
#' derivative has a point mass at the kink, and an average of its pointwise
#' values does not contain that mass. For the Laplace location,
#' \eqn{\ell_{\mu\mu}} is zero almost everywhere, so its pointwise average is
#' zero, while the information is \eqn{1/\sigma^2}. The information of such a
#' family is the variance of the score, which the score identity evaluates.
#' For a family that declares a non-smooth parameter (see
#' [param_smoothness()]), the Hessian identity at order 2, whether by
#' [expect_integrate()] or by [expect_mc()], signals an error that names
#' `identity = "score"`. [expect_opg()], which uses the score identity, and
#' `identity = "score"` remain valid. When the family declares where its kink
#' is, through [kink_decomposition()], the quadrature of [expectation()] places
#' a knot there.
#'
#' @return Nothing. This page documents the argument `method` shared by the
#'   functions named above; the value returned is theirs.
#'
#' @examples
#' # A family defined by its density alone has no closed-form information.
#' Dens <- S7::new_class("Dens", parent = continuous_distrib)
#' S7::method(distrib_pdf, Dens) <- function(distrib, y, theta, log = FALSE, ...) {
#'   stats::dnorm(y, theta[[1]], theta[[2]], log = log)
#' }
#' d <- Dens(distrib_name = "dens", dimension = "univariate",
#'           bounds = c(-Inf, Inf), params = c("mu", "sigma"),
#'           params_interpretation = c(mu = "location", sigma = "scale"),
#'           n_params = 2L,
#'           params_bounds = list(mu = c(-Inf, Inf), sigma = c(0, Inf)),
#'           link_params = list(mu = linkfunctions7::identity_link(),
#'                              sigma = linkfunctions7::log_link()),
#'           params_smooth = c(mu = TRUE, sigma = TRUE))
#' th <- list(mu = 0, sigma = 2)
#' y <- c(-1, 0.5, 2)
#'
#' # The sigma-sigma entry of the expected Hessian is -2 / sigma^2 = -0.5.
#' # The two quadratures agree; the empirical outer product of the scores is
#' # the value at one observation, and Monte Carlo carries sampling error.
#' methods <- list(
#'   opg = expect_opg(),
#'   integrate_hessian = expect_integrate(identity = "hessian"),
#'   integrate_score = expect_integrate(identity = "score"),
#'   mc_hessian = expect_mc(nsim = 2000, identity = "hessian"),
#'   mc_score = expect_mc(nsim = 2000, identity = "score")
#' )
#' set.seed(3)
#' vapply(methods, function(m) {
#'   distrib_expected_hessian(d, y, th, method = m)$sigma_sigma[1]
#' }, numeric(1))
#'
#' # An object is used also by a family that has a closed form.
#' g <- gaussian1_distrib()
#' distrib_expected_hessian(g, 0, list(mu = 0, sigma = 2))$sigma_sigma
#' distrib_expected_hessian(g, 0, list(mu = 0, sigma = 2),
#'                          method = expect_integrate())$sigma_sigma
#'
#' @seealso [expect_opg()], [expect_integrate()] and [expect_mc()] for the
#'   constructors; [ExpectMethod] for the classes;
#'   [distrib_expected_hessian()], [distrib_deriv3()] and
#'   [distrib_deriv4()], the generics that take `method`;
#'   [distrib_dexpected_hessian()] for the derivatives of the expected
#'   information; [fisher_scoring()], where it is set for a fit;
#'   [expected_info()] and [observed_info()], the objects that select the
#'   information of a fit, and [information_objects] for their comparison;
#'   [expected_by_bartlett()], [expected_by_integrate()] and
#'   [expected_by_mc()] for the implementations.
NULL

#' @title Observed Derivatives of a Given Order
#'
#' @description
#' Routes to [distrib_gradient()], [distrib_hessian()], [distrib_deriv3()],
#' [distrib_deriv4()], [distrib_deriv5()] or [distrib_deriv6()] according to
#' `order`, so that code
#' working at an order fixed only at run time does not have to branch. The
#' methods of [expected_derivative_methods] read it: [expect_integrate()] and
#' [expect_mc()] for the quantity that they average, and the score identity for
#' the lower orders that it multiplies together.
#'
#' @param distrib An object inheriting from `distrib`.
#' @param y A numeric vector of observations.
#' @param theta A named list of parameters, aligned to `distrib@params`.
#' @param order The derivative order, a single integer from 1 to 6. Anything
#'   else signals an error.
#' @param scale `"parameter"` (the default) or `"link"`, passed to the generic.
#'
#' @return A named list of derivative component vectors; each has length
#'   `length(y)`. The keys are `distrib@params` at order 1, [hess_names()] at
#'   order 2 (diagonal first) and [deriv_names()] above it (lexicographic), so
#'   a caller pairing two orders must match by name.
#'
#' @seealso [expected_derivative_methods], which reads it;
#'   [deriv_names()] and [hess_names()] for the two keyings;
#'   [expected_derivative()], the router above it.
#' @keywords internal
observed_deriv <- function(distrib, y, theta, order, scale = "parameter") {
  switch(as.character(order),
    "1" = distrib_gradient(distrib, y, theta, scale = scale),
    "2" = distrib_hessian(distrib, y, theta, scale = scale),
    "3" = distrib_deriv3(distrib, y, theta, scale = scale),
    "4" = distrib_deriv4(distrib, y, theta, scale = scale),
    "5" = distrib_deriv5(distrib, y, theta, scale = scale),
    "6" = distrib_deriv6(distrib, y, theta, scale = scale),
    stop("Unsupported derivative order: ", order, call. = FALSE)
  )
}

#' Stop at the Hessian Identity of a Family With a Kink
#'
#' @description
#' Signals an error when the Hessian identity is requested at order 2 for a
#' family that has a non-smooth parameter, because the observed second
#' derivative then lacks the point mass at the kink.
#'
#' @param distrib An object inheriting from `distrib`.
#' @param order The derivative order.
#' @param who The name of the constructor, used in the message.
#'
#' @return `NULL`, invisibly.
#'
#' @seealso [expected_derivative_methods]
#' @keywords internal
check_kink_identity <- function(distrib, order, who) {
  if (order == 2L && !all(param_smoothness(distrib))) {
    stop(sprintf(paste0(
      "'%s' has a parameter in which the log-density is not smooth, so the\n",
      "  observed second derivative lacks the point mass at the kink and its\n",
      "  expectation is not the information. Use %s(identity = \"score\")."),
      distrib@distrib_name, who), call. = FALSE)
  }
  invisible(NULL)
}

#' @title Expected Derivatives by Numerical Integration
#'
#' @description
#' Integrates each observed derivative component directly against the density
#' through [expectation()], component by component: the Hessian identity of
#' [expected_derivative_methods]. It is deterministic and is normally the most
#' accurate route where the observed derivative is available in closed form.
#'
#' @details
#' The function estimates \eqn{\mathbb{E}[\partial^k \ell]} literally. For a
#' regular model that is the required quantity. For a family with a non-smooth
#' parameter it is not at order 2: the observed \eqn{\ell_{ij}} has a point
#' mass at the kink, which an integral of its pointwise values does not contain
#' (on a Laplace the location component would be 0 instead of
#' \eqn{-1/\sigma^2}). The function signals an error there, and the score
#' identity of [expected_by_bartlett()] is the route that applies.
#'
#' Quadrature is also unreliable where the observed derivative is itself a
#' finite difference, since it then integrates numerical noise. The function
#' catches that error and signals it again, naming the component and pointing
#' at the alternatives, so that the caller receives a statement about the route
#' rather than an opaque failure from the integrator.
#'
#' @param distrib An object inheriting from `distrib`.
#' @param y A numeric vector of observations. Only its length is used, to
#'   recycle the result: an expectation does not depend on the sample.
#' @param theta A named list of parameters, aligned to `distrib@params`.
#' @param order The derivative order, a single integer from 2 to 4.
#' @param method An [ExpectIntegrate] object, whose tolerances and limits
#'   govern the quadrature and the series summation. The default is
#'   [expect_integrate()].
#'
#' @return A named list of expected derivative component vectors; each has
#'   length `length(y)`, and the names are the keys that [observed_deriv()]
#'   uses at that order.
#'
#' @seealso [expected_derivative_methods()] for the methods;
#'   [expected_by_bartlett()], the route that survives a kink;
#'   [expected_by_mc()] where quadrature performs poorly;
#'   [expectation()], which does the integration.
#' @keywords internal
expected_by_integrate <- function(distrib, y, theta, order,
                                  method = expect_integrate()) {
  old <- set_integration_settings(method)
  on.exit(set_integration_settings(old))
  check_kink_identity(distrib, order, "expect_integrate")
  n <- length(y)
  nms <- if (order == 2L) hess_names(distrib@params) else deriv_names(distrib@params, order)
  out <- lapply(nms, function(nm) {
    v <- tryCatch(
      expectation(
        distrib,
        function(y, theta) observed_deriv(distrib, y, theta, order)[[nm]],
        theta
      ),
      error = function(e) {
        stop(sprintf(
          paste0(
            "expect_integrate() failed on component '%s' of the order-%d ",
            "expected derivative:\n  %s\n",
            "Quadrature is unreliable here, which typically happens when the ",
            "observed derivative is itself obtained by finite differences. ",
            "Try method = expect_integrate(identity = \"score\") ",
            "(deterministic, uses only lower-order derivatives) or ",
            "method = expect_mc()."
          ),
          nm, order, conditionMessage(e)
        ), call. = FALSE)
      }
    )
    rep(v, length.out = n)
  })
  names(out) <- nms
  out
}

#' Index Tuples and Partitions of the Score Identity
#'
#' @description
#' Returns, for an order, the component names, the parameter index tuple
#' behind each name and the partitions of the index set into at least two
#' blocks, which are the ingredients of the score identity of
#' [expected_derivative_methods].
#'
#' @details
#' The index tuples come from the enumeration that produced the names.
#' Splitting the names on `"_"` fails for a parameter whose own name contains
#' one; see [deriv_indices()]. At order 2 the names come from [hess_names()],
#' which orders the diagonal first, so the tuples come from its inverse and not
#' from the lexicographic enumeration.
#'
#' @param params The parameter names of the family.
#' @param order The derivative order, 2 to 4.
#'
#' @return A list with the elements `nms`, `idx_of` and `parts`.
#'
#' @seealso [expected_by_bartlett()], [score_identity_terms()]
#' @keywords internal
score_identity_plan <- function(params, order) {
  nms <- if (order == 2L) hess_names(params) else deriv_names(params, order)
  idx_of <- if (order == 2L) {
    lapply(hess_pairs(params), function(pr) match(pr, params))
  } else {
    deriv_indices(params, order)
  }
  parts <- Filter(function(pp) length(pp) > 1L, numericals7::set_partitions(order))
  list(nms = nms, idx_of = idx_of, parts = parts)
}

#' The Score Identity Evaluated at Observations
#'
#' @description
#' Returns, for each component of the expected derivative of an order, the
#' quantity \eqn{-\sum_{\pi,\,|\pi| \ge 2} \prod_{B \in \pi} \ell_B} evaluated
#' at the observations `y` and the parameters `theta`, with no expectation. It
#' is the integrand of the score identity of [expected_derivative_methods], and
#' it is what [expect_opg()] reads and what [expect_mc()] averages when
#' `identity = "score"`.
#'
#' @param distrib An object inheriting from `distrib`.
#' @param y A numeric vector of observations.
#' @param theta A named list of parameters.
#' @param order The derivative order, 2 to 4.
#'
#' @return A named list of component vectors, each of the length of `y`,
#'   named as [observed_deriv()] names that order.
#'
#' @seealso [expected_by_opg()], [expected_by_mc()]
#' @keywords internal
score_identity_terms <- function(distrib, y, theta, order) {
  params <- distrib@params
  plan <- score_identity_plan(params, order)
  obs <- lapply(seq_len(order - 1L), function(k) observed_deriv(distrib, y, theta, k))
  out <- lapply(seq_along(plan$nms), function(t) {
    idx <- plan$idx_of[[t]]
    acc <- 0
    for (pp in plan$parts) {
      prod_val <- 1
      for (B in pp) {
        blk <- sort(idx[B])
        key <- paste(params[blk], collapse = "_")
        prod_val <- prod_val * obs[[length(blk)]][[key]]
      }
      acc <- acc + prod_val
    }
    -acc
  })
  names(out) <- plan$nms
  out
}

#' Expected Derivatives by Monte Carlo
#'
#' @description
#' Simulates `nsim` draws and averages over them the quantity of the identity
#' that is requested: the observed derivative for `identity = "hessian"`, and
#' the integrand of the score identity for `identity = "score"`.
#'
#' @details
#' One simulation is run for each row of `theta` after recycling.
#' Observations whose parameters coincide are not merged, so a `theta` that
#' varies over the observations costs one simulation per observation, and a
#' scalar `theta` costs one simulation however long `y` is.
#'
#' Estimates the same quantity as [expected_by_integrate()] (Hessian identity)
#' or [expected_by_bartlett()] (score identity), has an error that falls as
#' \eqn{1/\sqrt{n_{sim}}}, and is stochastic unless the seed is fixed. Its cost
#' is the cost of [distrib_rng()], so it is expensive for a distribution whose
#' generator inverts a slow quantile function. The Hessian identity at order 2
#' for a family with a non-smooth parameter signals an error, as in
#' [expected_by_integrate()].
#'
#' @param distrib An object inheriting from class `"distrib"`.
#' @param y A numeric vector of observations; only its length is used.
#' @param theta A named list of parameters.
#' @param order The derivative order, 2 to 4.
#' @param nsim The number of draws per parameter configuration.
#' @param identity `"hessian"` or `"score"`.
#'
#' @return A named list of expected derivative component vectors; each has
#'   length `length(y)`.
#'
#' @seealso [expected_derivative_methods()]
#' @keywords internal
expected_by_mc <- function(distrib, y, theta, order, nsim,
                           identity = "hessian") {
  score <- identical(identity, "score")
  if (!score) check_kink_identity(distrib, order, "expect_mc")
  n <- length(y)
  nms <- if (order == 2L) hess_names(distrib@params) else deriv_names(distrib@params, order)
  p <- distrib@n_params

  # one simulation per distinct parameter configuration
  th_rows <- transpose_params(expand_params(theta[seq_len(p)]))
  vals <- lapply(th_rows, function(r) {
    th <- as.list(r)
    ys <- distrib_rng(distrib, nsim, th)
    d <- if (score) score_identity_terms(distrib, ys, th, order)
         else observed_deriv(distrib, ys, th, order)
    vapply(nms, function(nm) base::mean(d[[nm]]), numeric(1))
  })
  m <- do.call(rbind, vals)

  out <- lapply(seq_along(nms), function(k) rep(unname(m[, k]), length.out = n))
  names(out) <- nms
  out
}

#' The Empirical Outer Product of the Scores
#'
#' @description
#' Returns the score identity of [expected_derivative_methods] evaluated at
#' each observation, with no expectation. At order 2 this is
#' \eqn{-\ell^{(i)}\ell^{(j)}}, the per-observation outer product of the score
#' (the BHHH estimator of the information); at orders 3 and 4 it is the sum of
#' products
#' \eqn{-\sum_{\pi,\,|\pi| \ge 2} \prod_{B \in \pi} \ell_B}. Unlike every other
#' route it takes no expectation: it integrates nothing and sums nothing over
#' the support, so it costs the derivatives of orders 1 to `order - 1` at the
#' observations.
#'
#' @details
#' This route differs from the identity from which it comes. The second
#' Bartlett identity states that
#' \eqn{\mathbb{E}[\ell^{(ij)}] = -\mathbb{E}[\ell^{(i)}\ell^{(j)}]}, and
#' [expected_by_bartlett()] evaluates the right-hand side as written, which
#' for a discrete family is a sum over the whole support and for a continuous
#' one a quadrature. This function drops the expectation and reads the
#' integrand at the observation in hand. The two agree in expectation and
#' differ in cost: for a discrete family the identity evaluates the derivative
#' kernel over the whole support for every observation, whereas this route
#' evaluates it once per observation.
#'
#' The estimate has more variance than the exact expectation and remains
#' consistent: each term is an unbiased estimate of the corresponding
#' component of the expected derivative at that observation, so the sum over
#' observations estimates it consistently. Component by component at order 2 it
#' is a rank-one matrix and correlates poorly with the exact expectation;
#' summed over the observations it approximates the exact information, and a
#' fit driven by it reaches the same maximum. At orders 3 and 4 the terms are
#' products of up to four derivatives, and the variance is larger still.
#'
#' It is the right default for a scoring step because it is a sum of outer
#' products and is therefore positive semidefinite by construction, whereas the
#' observed Hessian need not be positive semidefinite away from the optimum, and
#' because the fixed point is unchanged, since the score is exact and any
#' positive definite matrix takes a scoring iteration to the same stationary
#' point. The route is not suited to a reported standard error, which is read
#' off the observed information instead; see [fit_distrib()] and its `vcov()`
#' method.
#'
#' @param distrib An object inheriting from class `"distrib"`.
#' @param y A numeric vector of observations.
#' @param theta A named list of parameters.
#' @param order The derivative order, 2 to 4. The default is 2.
#'
#' @return A named list of expected derivative component vectors, named as
#'   [observed_deriv()] names that order; each has the length of `y`.
#'
#' @seealso [expected_by_bartlett()] for the identity that this function
#'   estimates, [expected_derivative_methods()] for the catalog of routes, and
#'   [expected_hessian_exact()] for the predicate that tests whether a family
#'   needs any of them.
#' @keywords internal
expected_by_opg <- function(distrib, y, theta, order = 2L) {
  n <- length(y)
  if (order != 2L) {
    out <- score_identity_terms(distrib, y, theta, order)
    return(lapply(out, rep_len, length.out = n))
  }
  g <- distrib_gradient(distrib, y, theta)
  pairs <- hess_pairs(distrib@params)
  out <- lapply(pairs, function(pr) {
    rep_len(-g[[pr[[1L]]]] * g[[pr[[2L]]]], n)
  })
  stats::setNames(out, hess_names(distrib@params))
}

#' Expected Derivatives by the Bartlett Identity
#'
#' @description
#' Obtains the expected derivative of a given order from expectations of
#' *products of lower-order* derivatives, by summing over set partitions of
#' the index set: the score identity of [expected_derivative_methods], computed
#' by quadrature or by summation over the support.
#'
#' @details
#' The identity of order \eqn{k} states that
#' \deqn{\sum_{\pi} \mathbb{E}\left[\prod_{B \in \pi} \ell_B\right] = 0,}
#' where the sum runs over every partition \eqn{\pi} of the index set. The
#' single-block partition is the target, so it equals minus the sum of all the
#' others, so that [numericals7::set_partitions()] is the whole algorithm and
#' the top-order derivative is never needed.
#'
#' At order 2 this reduces to
#' \eqn{\mathbb{E}[\ell_{ij}] = -\mathbb{E}[\ell_i \ell_j]}, which needs only
#' the score and remains the information where the log-likelihood has a kink:
#' there the pointwise average of \eqn{\partial^2 \ell} is not the information,
#' while the score variance still is.
#'
#' @param distrib An object inheriting from class `"distrib"`.
#' @param y A numeric vector of observations; only its length is used.
#' @param theta A named list of parameters.
#' @param order The derivative order, 2 to 4.
#'
#' @return A named list of expected derivative component vectors; each has
#'   length `length(y)`.
#'
#' @seealso [expected_derivative_methods()], [numericals7::set_partitions()]
#' @keywords internal
expected_by_bartlett <- function(distrib, y, theta, order) {
  params <- distrib@params
  n <- length(y)
  plan <- score_identity_plan(params, order)
  nms <- plan$nms
  idx_of <- plan$idx_of
  parts <- plan$parts

  out <- vector("list", length(nms))
  for (t in seq_along(nms)) {
    idx <- idx_of[[t]]
    acc <- 0
    for (pp in parts) {
      integrand <- function(y, theta) {
        prod_val <- 1
        for (B in pp) {
          blk <- sort(idx[B])
          key <- paste(params[blk], collapse = "_")
          prod_val <- prod_val * observed_deriv(distrib, y, theta, length(blk))[[key]]
        }
        prod_val
      }
      acc <- acc + expectation(distrib, integrand, theta)
    }
    out[[t]] <- rep(-acc, length.out = n)
  }

  names(out) <- nms
  out
}

#' Compute an Expected Derivative by a Method
#'
#' @description
#' Validates `method` and routes to the computation that the object selects. It
#' is the shared entry point for the expected derivatives that a family does not
#' write out, and for those that a caller requests with an explicit method. A
#' `NULL` method is resolved to the default of the order: [expect_opg()] at
#' order 2 and [expect_integrate()] above it.
#'
#' @details
#' The computation of each class is the method of the internal generic
#' `expected_by_method()`; see [expected_derivative_methods] for what the three
#' classes compute, their costs and their defaults.
#'
#' @param distrib An object inheriting from class `"distrib"`.
#' @param y A numeric vector of observations.
#' @param theta A named list of parameters.
#' @param order The derivative order, 2 to 4.
#' @param method An object inheriting from [ExpectMethod], or `NULL` (the
#'   default). Any other value, a character string included, signals an error.
#'
#' @return A named list of expected derivative component vectors.
#'
#' @seealso [expected_derivative_methods]
#' @keywords internal
expected_derivative <- function(distrib, y, theta, order, method = NULL) {
  check_method(method)
  expected_by_method(resolve_method(method, order), distrib, y, theta, order)
}
