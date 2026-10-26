#' @include distrib.R generics.R utility_functions.R numerical_functions.R
NULL

#' Approximations of an Expected Derivative
#'
#' @name expected_derivative_methods
#'
#' @description
#' When a distribution does not supply a closed-form expected derivative, the
#' expectation has to be approximated. The generics
#' [distrib_expected_hessian()], [distrib_deriv3()] and [distrib_deriv4()], the
#' derivatives of the expected information, and [fisher_scoring()] therefore
#' take an `approx` argument, an object that selects how the expectation is
#' taken. The argument is ignored when the distribution provides an analytical
#' method, which is always preferred. A character string is not accepted.
#'
#' The four approximations are built by [expect_opg()], [expect_bartlett()],
#' [expect_integrate()] and [expect_mc()], and they inherit from the class
#' [ExpectApprox].
#'
#' The expected information of a fit, the matrix whose inverse gives the
#' standard errors when [expected_info()] is selected, is built from the
#' expected Hessian. An [expected_info()] object carries one of the four
#' approximations in its `approx` property, and the object [observed_info()]
#' takes none; [information_objects] compares the two information objects.
#'
#' @details
#' Let \eqn{\ell} be the log-density and \eqn{\ell_{i}, \ell_{ij}, \dots} its
#' derivatives with respect to the parameters.
#'
#' \describe{
#'   \item{[expect_opg()]}{At order 2, reads \eqn{-\ell_i \ell_j} at each
#'     observation and takes no expectation, at the cost of one score
#'     evaluation; see [expected_by_opg()]. Above order 2 it is routed to the
#'     computation of [expect_bartlett()]. The matrix is positive semidefinite,
#'     and it is the least accurate of the four as an estimate of the
#'     information at one observation.}
#'
#'   \item{[expect_bartlett()]}{Uses the Bartlett
#'     identity of the corresponding order, which expresses the expected
#'     derivative through expectations of *products of lower-order*
#'     derivatives. Writing the identity as a sum over set partitions of the
#'     index set,
#'     \deqn{\sum_{\pi} \mathbb{E}\left[\prod_{B \in \pi} \ell_B\right] = 0,}
#'     the target term (the single-block partition) is obtained from all the
#'     others. At order 2 this is exactly the familiar outer product of
#'     gradients, \eqn{\mathbb{E}[\ell_{ij}] = -\mathbb{E}[\ell_i \ell_j]};
#'     at order 3, \eqn{\mathbb{E}[\ell_{ijk}] = -\left(\mathbb{E}[\ell_{ij}\ell_k] +
#'     \mathbb{E}[\ell_{ik}\ell_j] + \mathbb{E}[\ell_{jk}\ell_i] +
#'     \mathbb{E}[\ell_i\ell_j\ell_k]\right)}, and so on.
#'     At order 2 it needs only the score, whose product it integrates or
#'     sums, so it is slower than [expect_opg()]. It stays valid when the
#'     log-likelihood is not differentiable in a parameter, where
#'     \eqn{\mathbb{E}[H]} degenerates but the score variance still gives the
#'     information (see [laplace_distrib()]).
#'     It is deterministic. At higher orders it needs several integrals, so it
#'     is usually slower than [expect_integrate()]. It never requires the
#'     top-order derivative itself, which is useful where only the lower orders
#'     are available in closed form.}
#'
#'   \item{[expect_integrate()]}{Integrates the observed derivative of that
#'     order directly against the density (numerical quadrature for continuous
#'     distributions, series summation for discrete ones). It is deterministic
#'     and normally the most accurate when the observed derivative is available
#'     in closed form. It estimates \eqn{\mathbb{E}[\partial^k \ell]} literally,
#'     with one exception: at order 2, for a family with a non-smooth
#'     parameter (see [param_smoothness()]), it evaluates
#'     \eqn{-\mathbb{E}[\ell_i \ell_j]} as [expect_bartlett()] does, because the
#'     observed \eqn{\ell_{ij}} then lacks the point mass at the kink. The
#'     tolerances and limits of the quadrature and of the series are
#'     arguments of the constructor, and their defaults are those of
#'     [numericals7::quad_vec()] and [numericals7::series_vec()].}
#'
#'   \item{[expect_mc()]}{Simulates `nsim` observations from the distribution
#'     and averages the observed derivative over them. The simplest and most
#'     robust option when quadrature performs poorly (heavy tails, awkward
#'     supports), but it is stochastic: the error decreases only as
#'     \eqn{1/\sqrt{n_{sim}}}, so results are not exactly reproducible unless
#'     the seed is fixed. Estimates
#'     the same quantity as [expect_integrate()].
#'     The cost of this option is the cost of simulating from the
#'     distribution. A native [distrib_rng()] is best; failing that, the
#'     default generator uses inverse transform when the family has a quantile
#'     method of its own and [rng_grou()], which only needs the density,
#'     otherwise. It becomes expensive only for a distribution that supplies a
#'     quantile function slow enough to make inverse transform sampling the
#'     bottleneck; prefer [expect_bartlett()] or [expect_integrate()] in that
#'     case.}
#' }
#'
#' # Defaults
#' [distrib_expected_hessian()] and the generics for the derivatives of the
#' expected information default to [expect_opg()], which needs one score
#' evaluation per observation and no expectation. [distrib_deriv3()] and
#' [distrib_deriv4()] default to [expect_integrate()], since at those orders
#' direct integration of the available derivative is usually cheaper and more
#' accurate. [fisher_scoring()] defaults to [expect_opg()], and
#' [fit_distrib()] rejects an approximation other than that default for a
#' family with a closed-form expected information. [expected_info()]
#' defaults to [expect_opg()], and the `vcov()` method of a fit defaults to
#' [observed_info()], which takes no approximation. Passing [expect_opg()] to [distrib_deriv3()]
#' or [distrib_deriv4()] gives the computation of [expect_bartlett()], because
#' the outer product of the scores has no counterpart above order 2.
#'
#' # A kink in the log-density
#' Where the log-density has a kink in a
#' parameter, the observed second derivative has a point mass at the kink, and
#' an average of its pointwise values does not contain that mass. For the
#' Laplace location, \eqn{\ell_{\mu\mu}} is zero almost everywhere, so its
#' pointwise average is zero, while the information is \eqn{1/\sigma^2}. The
#' information of such a family is the variance of the score, and the second
#' Bartlett identity evaluates that variance. [expect_bartlett()],
#' [expect_opg()] and, at order 2, [expect_integrate()] therefore return
#' \eqn{-1/\sigma^2} per observation; [expect_mc()] averages the observed
#' \eqn{\ell_{\mu\mu}} and returns 0. When the family declares where its kink
#' is, through [kink_decomposition()], the quadrature of [expectation()] places
#' a knot there.
#'
#' @return Nothing. This page documents the `approx` argument shared by the
#'   generics named above; the value returned is theirs.
#'
#' @examples
#' # A family defined by its density alone reads 'approx'. The two
#' # deterministic approximations agree on the exact -2 / sigma^2; the outer
#' # product of the scores reads one observation's score product, and Monte
#' # Carlo carries sampling error.
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
#' approximations <- list(opg = expect_opg(), bartlett = expect_bartlett(),
#'                        integrate = expect_integrate())
#' vapply(approximations, function(a) {
#'   distrib_expected_hessian(d, y, th, approx = a)$sigma_sigma[1]
#' }, numeric(1))
#'
#' set.seed(3)
#' distrib_expected_hessian(d, y, th, approx = expect_mc(nsim = 2000))$sigma_sigma[1]
#'
#' # A family that writes its expected information out ignores the argument,
#' # and fit_distrib() rejects one given there rather than dropping it.
#' g <- gaussian1_distrib()
#' identical(distrib_expected_hessian(g, c(-1, 0, 1), list(mu = 0, sigma = 1)),
#'           distrib_expected_hessian(g, c(-1, 0, 1), list(mu = 0, sigma = 1),
#'                                    approx = expect_mc()))
#'
#' @seealso [expect_opg()], [expect_bartlett()], [expect_integrate()] and
#'   [expect_mc()] for the constructors; [ExpectApprox] for the classes;
#'   [distrib_expected_hessian()], [distrib_deriv3()] and
#'   [distrib_deriv4()], the generics that take `approx`;
#'   [fisher_scoring()], where it is set for a fit;
#'   [expected_info()] and [observed_info()], the objects that select the
#'   information of a fit, and [information_objects] for their comparison;
#'   [expected_by_bartlett()], [expected_by_integrate()] and
#'   [expected_by_mc()] for the implementations.
NULL

#' @title Observed Derivatives of a Given Order
#'
#' @description
#' Routes to [distrib_gradient()], [distrib_hessian()], [distrib_deriv3()] or
#' [distrib_deriv4()] according to `order`, so that code working at an order
#' fixed only at run time does not have to branch. The approximations of
#' [expected_derivative_methods] other than [expect_opg()] read it:
#' [expect_integrate()] and [expect_mc()] for the quantity that they average,
#' and [expect_bartlett()] for the lower orders that its identity multiplies
#' together.
#'
#' @param distrib An object inheriting from `distrib`.
#' @param y A numeric vector of observations.
#' @param theta A named list of parameters, aligned to `distrib@params`.
#' @param order The derivative order, a single integer from 1 to 4. Anything
#'   else signals an error.
#'
#' @return A named list of derivative component vectors; each has length
#'   `length(y)`. The keys are `distrib@params` at order 1, [hess_names()] at
#'   order 2 (diagonal first) and [deriv_names()] above it (lexicographic), so
#'   a caller pairing two orders must match by name.
#'
#' @seealso [expected_derivative_methods()], which reads it;
#'   [deriv_names()] and [hess_names()] for the two keyings;
#'   [expected_derivative()], the router above it.
#' @keywords internal
observed_deriv <- function(distrib, y, theta, order) {
  switch(as.character(order),
    "1" = distrib_gradient(distrib, y, theta),
    "2" = distrib_hessian(distrib, y, theta),
    "3" = distrib_deriv3(distrib, y, theta),
    "4" = distrib_deriv4(distrib, y, theta),
    stop("Unsupported derivative order: ", order, call. = FALSE)
  )
}

#' @title Expected Derivatives by Numerical Integration
#'
#' @description
#' Integrates each observed derivative component directly against the density
#' through [expectation()], component by component. It is deterministic and is
#' normally the most accurate route where the observed derivative is available
#' in closed form.
#'
#' @details
#' The function estimates \eqn{\mathbb{E}[\partial^k \ell]} literally. For a
#' regular model that is the required quantity. For a family with a non-smooth
#' parameter it is not: the observed \eqn{\ell_{ij}} has a point mass at the
#' kink, which an integral of its pointwise values does not contain (on a
#' Laplace the location component would be 0 instead of
#' \eqn{-1/\sigma^2}). At order 2 such a family is therefore passed to
#' [expected_by_bartlett()], which needs only the score.
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
#' @param approx An [ExpectIntegrate] object, whose tolerances and limits
#'   govern the quadrature and the series summation. The default is
#'   [expect_integrate()].
#'
#' @return A named list of expected derivative component vectors; each has
#'   length `length(y)`, and the names are the keys that [observed_deriv()]
#'   uses at that order.
#'
#' @seealso [expected_derivative_methods()] for the four approximations;
#'   [expected_by_bartlett()], the route that survives a kink;
#'   [expected_by_mc()] where quadrature performs poorly;
#'   [expectation()], which does the integration.
#' @keywords internal
expected_by_integrate <- function(distrib, y, theta, order,
                                  approx = expect_integrate()) {
  old <- set_integration_settings(approx)
  on.exit(set_integration_settings(old))
  if (order == 2L && !all(param_smoothness(distrib))) {
    return(expected_by_bartlett(distrib, y, theta, 2L))
  }
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
            "Try approx = expect_bartlett() (deterministic, uses only lower-order ",
            "derivatives) or approx = expect_mc()."
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

#' Expected Derivatives by Monte Carlo
#'
#' @description
#' Simulates `nsim` draws and averages the observed derivative over them.
#'
#' @details
#' One simulation is run for each row of `theta` after recycling.
#' Observations whose parameters coincide are not merged, so a `theta` that
#' varies over the observations costs one simulation per observation, and a
#' scalar `theta` costs one simulation however long `y` is.
#'
#' Estimates the same quantity as [expected_by_integrate()], has an
#' error that falls as \eqn{1/\sqrt{n_{sim}}}, and is stochastic unless the seed
#' is fixed. Its cost is the cost of [distrib_rng()], so it is expensive for a
#' distribution whose generator inverts a slow quantile function.
#'
#' @param distrib An object inheriting from class `"distrib"`.
#' @param y A numeric vector of observations; only its length is used.
#' @param theta A named list of parameters.
#' @param order The derivative order, 2 to 4.
#' @param nsim The number of draws per parameter configuration.
#'
#' @return A named list of expected derivative component vectors; each has
#'   length `length(y)`.
#'
#' @seealso [expected_derivative_methods()]
#' @keywords internal
expected_by_mc <- function(distrib, y, theta, order, nsim) {
  n <- length(y)
  nms <- if (order == 2L) hess_names(distrib@params) else deriv_names(distrib@params, order)
  p <- distrib@n_params

  # one simulation per distinct parameter configuration
  th_rows <- transpose_params(expand_params(theta[seq_len(p)]))
  vals <- lapply(th_rows, function(r) {
    th <- as.list(r)
    ys <- distrib_rng(distrib, nsim, th)
    d <- observed_deriv(distrib, ys, th, order)
    vapply(nms, function(nm) base::mean(d[[nm]]), numeric(1))
  })
  m <- do.call(rbind, vals)

  out <- lapply(seq_along(nms), function(k) rep(unname(m[, k]), length.out = n))
  names(out) <- nms
  out
}

#' The Information as the Outer Product of the Observed Scores
#'
#' @description
#' Returns \eqn{\ell^{(ij)} \approx -\ell^{(i)}\ell^{(j)}} evaluated at each
#' observation, the per-observation outer product of the score. It is the
#' BHHH estimator of the information, and unlike every other route here it
#' takes no expectation: it integrates nothing and sums nothing over the
#' support, so it costs one call to [distrib_gradient()].
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
#' consistent.
#' \eqn{-\ell^{(i)}\ell^{(j)}} is an unbiased estimate of the corresponding
#' component of \eqn{\mathbb{E}[\ell^{(ij)}]} at that observation, so the sum
#' over observations, which a scoring step aggregates into
#' \eqn{X^\top W X}, estimates the information consistently. Component by
#' component it is a rank-one matrix and correlates poorly with the exact
#' expectation; summed over the observations it approximates the exact
#' information, and a fit driven by it reaches the same maximum.
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
#' The function applies at order 2 only. The outer product of scores is the
#' second-order identity and
#' has no counterpart above it, so [expected_derivative()] routes a higher
#' order to [expected_by_bartlett()].
#'
#' @param distrib An object inheriting from class `"distrib"`.
#' @param y A numeric vector of observations.
#' @param theta A named list of parameters.
#'
#' @return A named list of expected Hessian component vectors, named by
#'   [hess_names()]; each has the length of `y`.
#'
#' @seealso [expected_by_bartlett()] for the identity that this function
#'   estimates, [expected_derivative_methods()] for the catalog of routes, and
#'   [expected_hessian_exact()] for the predicate that tests whether a family
#'   needs any of them.
#' @keywords internal
expected_by_opg <- function(distrib, y, theta) {
  g <- distrib_gradient(distrib, y, theta)
  n <- length(y)
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
#' the index set.
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
#' while the score variance still is. [expected_by_integrate()] calls this
#' function at order 2 for such a family.
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
  nms <- if (order == 2L) hess_names(params) else deriv_names(params, order)

  # Index tuple behind each component name, taken from the enumeration that
  # produced the names. Splitting the names on "_" fails for a parameter whose
  # own name contains one; see deriv_indices(). At order 2 the names come from
  # hess_names(), which orders the diagonal first, so the tuples must come from
  # its inverse rather than from the lexicographic enumeration.
  idx_of <- if (order == 2L) {
    lapply(hess_pairs(params), function(pr) match(pr, params))
  } else {
    deriv_indices(params, order)
  }

  # every partition except the single-block one (that is the target term)
  parts <- Filter(function(pp) length(pp) > 1L, numericals7::set_partitions(order))

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

#' Dispatch an Expected-Derivative Approximation
#'
#' @description
#' Validates `approx` and routes to the computation that the object selects. It
#' is the shared entry point for the fallback expected-derivative methods.
#'
#' @details
#' [expect_opg()] and [expect_bartlett()] are two readings of the second
#' Bartlett identity and are not the same computation. The identity equates
#' \eqn{\mathbb{E}[\ell^{(ij)}]} with \eqn{-\mathbb{E}[\ell^{(i)}\ell^{(j)}]};
#' the Bartlett route evaluates the expectation, which is a sum over the
#' support or a quadrature, while the outer product of the scores reads the
#' integrand at the observation and takes no expectation at all. They agree in
#' expectation and differ by orders of magnitude in cost, so [expect_opg()] is
#' the default at order 2.
#'
#' Above order 2 the outer product of scores is not an identity, so
#' [expect_opg()] is routed to the Bartlett computation and the function does
#' not signal an error. The argument is one setting for a whole call, and a
#' method of order 4 can receive the value that a method of order 2 received.
#' The computation for each class is the method of the internal generic
#' `expected_by_approx()`.
#'
#' @param distrib An object inheriting from class `"distrib"`.
#' @param y A numeric vector of observations.
#' @param theta A named list of parameters.
#' @param order The derivative order, 2 to 4.
#' @param approx An object inheriting from [ExpectApprox]: [expect_opg()] (the
#'   default), [expect_bartlett()], [expect_integrate()] or [expect_mc()]. Any
#'   other value, a character string included, signals an error.
#'
#' @return A named list of expected derivative component vectors.
#'
#' @seealso [expected_derivative_methods]
#' @keywords internal
expected_derivative <- function(distrib, y, theta, order,
                                approx = expect_opg()) {
  check_approx(approx)
  expected_by_approx(approx, distrib, y, theta, order)
}
