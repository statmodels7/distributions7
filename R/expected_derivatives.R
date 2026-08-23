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
#'     usually slower than `"integrate"`, but it never requires the
#'     top-order derivative itself -- useful when only the lower orders are
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
#' @return Nothing. This page documents the `approx` argument shared by
#'   the three generics named above; the value returned is theirs.
#'
#' @seealso [distrib_expected_hessian()], [distrib_deriv3()],
#'   [distrib_deriv4()]
NULL

#' Observed Derivatives of a Given Order
#'
#' @description
#' Routes to [distrib_gradient()], [distrib_hessian()],
#' [distrib_deriv3()] or [distrib_deriv4()] according to
#' `order`, so that code working at an order fixed only at run time does not
#' have to branch.
#'
#' @param distrib An object inheriting from class `"distrib"`.
#' @param y A numeric vector of observations.
#' @param theta A named list of parameters.
#' @param order The derivative order, 1 to 4.
#'
#' @return A named list of derivative component vectors, keyed as
#'   [hess_names()] at order 2 and [deriv_names()] above it.
#'
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

#' Expected Derivatives by Numerical Integration
#'
#' @description
#' Integrates each observed derivative component directly against the density,
#' component by component, through [expectation()].
#'
#' @details
#' This estimates \eqn{\mathbb{E}[\partial^k \ell]} literally. For a regular
#' model that is the quantity wanted; for a non-regular one it is not the
#' information, and [expected_by_bartlett()] is the route that stays
#' valid (see [laplace_distrib()]).
#'
#' Quadrature is unreliable when the observed derivative is itself a finite
#' difference, since it then integrates numerical noise, so the error is caught
#' and re-raised naming the component and pointing at the alternatives rather
#' than surfacing as an opaque failure from the integrator.
#'
#' @param distrib An object inheriting from class `"distrib"`.
#' @param y A numeric vector of observations; only its length is used, to recycle
#'   the result.
#' @param theta A named list of parameters.
#' @param order The derivative order, 2 to 4.
#'
#' @return A named list of expected derivative component vectors, each of length
#'   `length(y)`.
#'
#' @seealso [expected_derivative_methods()]
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
#' `"opg"` is folded into `"bartlett"` here rather than implemented
#' separately, because at order 2 the outer product of gradients *is* the
#' Bartlett identity; the two names describe the same computation from different
#' traditions.
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
                                approx = c("bartlett", "integrate", "mc", "opg"),
                                nsim = 10000) {
  approx <- match.arg(approx)
  if (approx == "opg") approx <- "bartlett"
  switch(approx,
    bartlett  = expected_by_bartlett(distrib, y, theta, order),
    integrate = expected_by_integrate(distrib, y, theta, order),
    mc        = expected_by_mc(distrib, y, theta, order, nsim)
  )
}
