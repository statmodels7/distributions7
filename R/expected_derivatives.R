#' @include distrib.R generics.R utility_functions.R numerical_functions.R
NULL

#' Strategies for Expected Derivatives
#'
#' @name expected_derivative_methods
#'
#' @description
#' When a distribution does not supply a closed-form expected derivative, the
#' expectation has to be approximated. The generics
#' [distrib_expected_hessian()], [distrib_deriv3()] and
#' [distrib_deriv4()] therefore accept an `approx` argument
#' selecting *how* the expectation is taken. The argument is **ignored**
#' when the distribution provides an analytical method (which is always preferred).
#'
#' @details
#' Let \eqn{\ell} be the log-density and \eqn{\ell_{i}, \ell_{ij}, \dots} its
#' derivatives with respect to the parameters. The three strategies are:
#'
#' \describe{
#'   \item{`"bartlett"` (also spelled `"opg"`)}{Uses the Bartlett
#'     identity of the corresponding order, which expresses the expected
#'     derivative through expectations of *products of lower-order*
#'     derivatives. Writing the identity as a sum over set partitions of the
#'     index set,
#'     \deqn{\sum_{\pi} \mathbb{E}\left[\prod_{B \in \pi} \ell_B\right] = 0,}
#'     the target term (the single-block partition) is obtained from all the
#'     others. At order 2 this is exactly the familiar **outer product of
#'     gradients**, \eqn{\mathbb{E}[\ell_{ij}] = -\mathbb{E}[\ell_i \ell_j]};
#'     at order 3, \eqn{\mathbb{E}[\ell_{ijk}] = -\left(\mathbb{E}[\ell_{ij}\ell_k] +
#'     \mathbb{E}[\ell_{ik}\ell_j] + \mathbb{E}[\ell_{jk}\ell_i] +
#'     \mathbb{E}[\ell_i\ell_j\ell_k]\right)}, and so on.
#'     *Fastest at order 2* (only the score is needed) and the only variant
#'     that stays valid when the log-likelihood is not differentiable in a
#'     parameter, where \eqn{\mathbb{E}[H]} degenerates but the score variance
#'     still gives the information (see [laplace_distrib()]).
#'     Deterministic. At higher orders it needs several integrals, so it is
#'     usually slower than `"integrate"`. It never requires the top-order
#'     derivative itself, which is useful where only the lower orders are
#'     available in closed form.}
#'
#'   \item{`"integrate"`}{Integrates the observed derivative of that order
#'     directly against the density (numerical quadrature for continuous
#'     distributions, series summation for discrete ones). Deterministic and
#'     normally the most accurate when the observed derivative is available in
#'     closed form. Estimates \eqn{\mathbb{E}[\partial^k \ell]} literally, which
#'     for a non-regular model is *not* the information.}
#'
#'   \item{`"mc"`}{Simulates `nsim` observations from the distribution
#'     and averages the observed derivative over them. The simplest and most
#'     robust option when quadrature struggles (heavy tails, awkward supports),
#'     but stochastic: the error decreases only as \eqn{1/\sqrt{n_{sim}}}, so
#'     results are not exactly reproducible unless the seed is fixed. Estimates
#'     the same quantity as `"integrate"`.
#'     **Note:** the cost of this option is the cost of simulating from the
#'     distribution. A native [distrib_rng()] is best, and failing that
#'     the default RNG falls back to [rng_grou()], which only needs the
#'     density. It becomes expensive only for a distribution that supplies a
#'     quantile function slow enough to make inverse transform sampling the
#'     bottleneck; prefer `"bartlett"` or `"integrate"` in that case.}
#' }
#'
#' **Defaults.** `distrib_expected_hessian` defaults to
#' `"bartlett"`, because at order 2 it is both the cheapest (only first
#' derivatives) and the most broadly valid. `distrib_deriv3` and
#' `distrib_deriv4` default to `"integrate"`, since at those orders
#' direct integration of the available derivative is usually cheaper and more
#' accurate.
#'
#' **What the kink costs, measured.** On a Laplace carrying a density, a score
#' and a Hessian but no expected method, at \eqn{\sigma = 1} over 200
#' observations: `"bartlett"` returns \eqn{-200}, which is \eqn{-n/\sigma^2}
#' and agrees with the shipped family's closed form to the digit, while
#' `"integrate"` and `"mc"` both return **exactly 0**. Neither is wrong about
#' what it computes. The observed \eqn{\ell_{\mu\mu}} really is zero almost
#' everywhere, so its expectation is zero; what fails is the identification of
#' that expectation with \eqn{-\mathcal{I}(\theta)}, which is the second
#' Bartlett identity. Only the score-based route survives, and the information
#' of a non-regular family is *defined* as the variance of the score.
#'
#' @return Nothing. This page documents the `approx` argument shared by the
#'   three generics named above; the value returned is theirs.
#'
#' @examples
#' # A family with no closed-form expected information reads 'approx'. The two
#' # deterministic strategies agree; Monte Carlo is the same quantity with
#' # sampling error on it.
#' sn <- skewnormal1_distrib()
#' th <- list(mu = 0, sigma = 1, alpha = 3)
#' set.seed(2)
#' y <- distrib_rng(sn, 40, th)
#'
#' vapply(c("bartlett", "integrate"), function(a) {
#'   sum(distrib_expected_hessian(sn, y, th, approx = a)$alpha_alpha)
#' }, numeric(1))
#'
#' set.seed(3)
#' sum(distrib_expected_hessian(sn, y, th, approx = "mc", nsim = 500)$alpha_alpha)
#'
#' # A family that writes its expected information out ignores the argument,
#' # and fit_distrib() rejects one given there rather than dropping it.
#' g <- gaussian1_distrib()
#' identical(distrib_expected_hessian(g, c(-1, 0, 1), list(mu = 0, sigma = 1)),
#'           distrib_expected_hessian(g, c(-1, 0, 1), list(mu = 0, sigma = 1),
#'                                    approx = "mc"))
#'
#' @seealso [distrib_expected_hessian()], [distrib_deriv3()] and
#'   [distrib_deriv4()], the three generics that take `approx`;
#'   [fisher_scoring()], where it is set for a fit;
#'   [expected_by_bartlett()], [expected_by_integrate()] and
#'   [expected_by_mc()] for the three implementations.
NULL

#' @title Observed Derivatives of a Given Order
#'
#' @description
#' Routes to [distrib_gradient()], [distrib_hessian()], [distrib_deriv3()] or
#' [distrib_deriv4()] according to `order`, so that code working at an order
#' fixed only at run time does not have to branch. Every strategy in
#' [expected_derivative_methods()] reads it, `"integrate"` and `"mc"` for the
#' quantity they average and `"bartlett"` for the lower orders its identity
#' multiplies together.
#'
#' @param distrib An object inheriting from `distrib`.
#' @param y A numeric vector of observations.
#' @param theta A named list of parameters, aligned to `distrib@params`.
#' @param order The derivative order, a single integer from 1 to 4. Anything
#'   else signals an error.
#'
#' @return A named list of derivative component vectors, each of length
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
#' What this estimates is \eqn{\mathbb{E}[\partial^k \ell]} literally. For a
#' regular model that is the quantity wanted. For a non-regular one it is not
#' the information, and the difference is total rather than small: on a Laplace
#' with no closed-form expected method, at \eqn{\sigma = 1} over 200
#' observations, this returns **exactly 0** for the location component while
#' [expected_by_bartlett()] returns \eqn{-200 = -n/\sigma^2}, which is the
#' information. The observed \eqn{\ell_{\mu\mu}} is zero almost everywhere, so
#' its integral against the density is zero and the point mass at the kink is
#' invisible to it.
#'
#' Quadrature is also unreliable where the observed derivative is itself a
#' finite difference, since it then integrates numerical noise. That error is
#' caught and re-raised naming the component and pointing at the alternatives,
#' so it reaches the caller as a statement about the route rather than as an
#' opaque failure from the integrator.
#'
#' @param distrib An object inheriting from `distrib`.
#' @param y A numeric vector of observations. Only its length is used, to
#'   recycle the result: an expectation does not depend on the sample.
#' @param theta A named list of parameters, aligned to `distrib@params`.
#' @param order The derivative order, a single integer from 2 to 4.
#'
#' @return A named list of expected derivative component vectors, each of
#'   length `length(y)`, keyed as [observed_deriv()] keys that order.
#'
#' @seealso [expected_derivative_methods()] for the three strategies;
#'   [expected_by_bartlett()], the route that survives a kink;
#'   [expected_by_mc()] where quadrature struggles;
#'   [expectation()], which does the integration.
#' @keywords internal
expected_by_integrate <- function(distrib, y, theta, order) {
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
            "approx = \"integrate\" failed on component '%s' of the order-%d ",
            "expected derivative:\n  %s\n",
            "Quadrature is unreliable here, which typically happens when the ",
            "observed derivative is itself obtained by finite differences. ",
            "Try approx = \"bartlett\" (deterministic, uses only lower-order ",
            "derivatives) or approx = \"mc\"."
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
#' One simulation is run per *distinct* parameter configuration rather than
#' per observation: the expectation depends on \eqn{\theta} alone, so observations
#' sharing a \eqn{\theta} share an answer, and a model with a scalar \eqn{\theta}
#' costs one simulation however long `y` is.
#'
#' Estimates the same quantity as [expected_by_integrate()], with an
#' error falling as \eqn{1/\sqrt{n_{sim}}}, and is stochastic unless the seed is
#' fixed. Its cost is the cost of sampling, so it is the wrong choice for a
#' distribution whose RNG falls back to inverse transform on a numerical quantile.
#'
#' @param distrib An object inheriting from class `"distrib"`.
#' @param y A numeric vector of observations; only its length is used.
#' @param theta A named list of parameters.
#' @param order The derivative order, 2 to 4.
#' @param nsim The number of draws per parameter configuration.
#'
#' @return A named list of expected derivative component vectors, each of length
#'   `length(y)`.
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
#' takes no expectation: nothing is integrated and nothing is summed over the
#' support, so it costs one call to [distrib_gradient()].
#'
#' @details
#' The name is worth keeping apart from the identity it comes from. The second
#' Bartlett identity states that
#' \eqn{\mathbb{E}[\ell^{(ij)}] = -\mathbb{E}[\ell^{(i)}\ell^{(j)}]}, and
#' [expected_by_bartlett()] evaluates the right-hand side as written, which
#' for a discrete family is a sum over the whole support and for a continuous
#' one a quadrature. This function drops the expectation and reads the
#' integrand at the observation in hand. The two agree in expectation and not
#' in cost: measured on `pig1_distrib()` with \eqn{\mu} varying by
#' observation, one call to the identity at \eqn{n = 1000} evaluated 1.67
#' million rows of the family's derivative kernel and took 18.4 seconds, where
#' this route evaluates \eqn{n} rows.
#'
#' What is given up is variance, not consistency.
#' \eqn{-\ell^{(i)}\ell^{(j)}} is an unbiased estimate of the corresponding
#' component of \eqn{\mathbb{E}[\ell^{(ij)}]} at that observation, so the sum
#' over observations, which is what a scoring step aggregates into
#' \eqn{X^\top W X}, estimates the information consistently. Component by
#' component it is a rank-one matrix and correlates poorly with the exact
#' expectation; summed it agreed with the exact route to within 10 per cent on
#' the same measurement, and a fit driven by it reached the same maximum.
#'
#' Two properties make it the right default for a scoring step. It is positive
#' semidefinite by construction, being a sum of outer products, where the
#' observed Hessian need not be away from the optimum; and the fixed point is
#' unchanged, because the score is exact and any positive definite matrix
#' takes a scoring iteration to the same stationary point. What it is not
#' suited to is a reported standard error, which is read off the observed
#' information instead -- see [fit_distrib()] and its `vcov()` method.
#'
#' Order 2 only. The outer product of scores is the second-order identity and
#' has no counterpart above it, so [expected_derivative()] routes a higher
#' order to [expected_by_bartlett()].
#'
#' @param distrib An object inheriting from class `"distrib"`.
#' @param y A numeric vector of observations.
#' @param theta A named list of parameters.
#'
#' @return A named list of expected Hessian component vectors, named by
#'   [hess_names()], each of the length of `y`.
#'
#' @seealso [expected_by_bartlett()] for the identity this estimates,
#'   [expected_derivative_methods()] for the catalog of routes, and
#'   [expected_hessian_exact()] for the predicate that says whether a family
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
#' the sum running over every partition \eqn{\pi} of the index set. The
#' single-block partition is the target, so it equals minus the sum of all the
#' others -- which is why [numericals7::set_partitions()] is the whole algorithm and
#' why the top-order derivative is never needed.
#'
#' At order 2 this reduces to the outer product of gradients,
#' \eqn{\mathbb{E}[\ell_{ij}] = -\mathbb{E}[\ell_i \ell_j]}, which is both the
#' cheapest route and the only one that survives a model where the
#' log-likelihood has a kink: there \eqn{\mathbb{E}[\partial^2 \ell]} genuinely
#' is not the information, while the score variance still is. That is why it is
#' the default at order 2 and why `"opg"` is accepted as a spelling of it.
#'
#' @param distrib An object inheriting from class `"distrib"`.
#' @param y A numeric vector of observations; only its length is used.
#' @param theta A named list of parameters.
#' @param order The derivative order, 2 to 4.
#'
#' @return A named list of expected derivative component vectors, each of length
#'   `length(y)`.
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

#' Dispatch an Expected-Derivative Strategy
#'
#' @description
#' Shared entry point for the fallback expected-derivative methods: validates
#' `approx` and routes to the chosen strategy.
#'
#' @details
#' `"opg"` and `"bartlett"` are two readings of the second Bartlett identity
#' and are NOT the same computation. The identity equates
#' \eqn{\mathbb{E}[\ell^{(ij)}]} with \eqn{-\mathbb{E}[\ell^{(i)}\ell^{(j)}]};
#' `"bartlett"` evaluates the expectation, which is a sum over the support or a
#' quadrature, while `"opg"` reads the integrand at the observation and takes
#' no expectation at all. They agree in expectation and differ by orders of
#' magnitude in cost, which is why `"opg"` is the default at order 2.
#'
#' Above order 2 the outer product of scores is not an identity for anything,
#' so `"opg"` there is routed to `"bartlett"` rather than refused: the
#' argument is one setting for a whole call and an order-4 method may ask for
#' the same strategy an order-2 one was given.
#'
#' @param distrib An object inheriting from class `"distrib"`.
#' @param y A numeric vector of observations.
#' @param theta A named list of parameters.
#' @param order The derivative order, 2 to 4.
#' @param approx One of `"bartlett"`, `"integrate"`, `"mc"` or
#'   `"opg"`.
#' @param nsim Number of draws, used only by `"mc"`.
#'
#' @return A named list of expected derivative component vectors.
#'
#' @seealso [expected_derivative_methods()]
#' @keywords internal
expected_derivative <- function(distrib, y, theta, order,
                                approx = c("opg", "bartlett", "integrate", "mc"),
                                nsim = 10000) {
  approx <- match.arg(approx)
  # the outer product of scores is the SECOND-order identity and has no
  # counterpart above it
  if (approx == "opg" && order != 2L) approx <- "bartlett"
  switch(approx,
    opg       = expected_by_opg(distrib, y, theta),
    bartlett  = expected_by_bartlett(distrib, y, theta, order),
    integrate = expected_by_integrate(distrib, y, theta, order),
    mc        = expected_by_mc(distrib, y, theta, order, nsim)
  )
}
