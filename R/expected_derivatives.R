#' @include distrib.R generics.R utility_functions.R numerical_functions.R
NULL

#' Strategies for Expected Derivatives
#'
#' @name expected_derivative_methods
#'
#' @description
#' When a distribution does not supply a closed-form expected derivative, the
#' expectation has to be approximated. The generics
#' \code{\link{distrib_expected_hessian}}, \code{\link{distrib_deriv3}} and
#' \code{\link{distrib_deriv4}} therefore accept an \code{approx} argument
#' selecting \emph{how} the expectation is taken. The argument is \strong{ignored}
#' when the distribution provides an analytical method (which is always preferred).
#'
#' @details
#' Let \eqn{\ell} be the log-density and \eqn{\ell_{i}, \ell_{ij}, \dots} its
#' derivatives with respect to the parameters. The three strategies are:
#'
#' \describe{
#'   \item{\code{"bartlett"} (also spelled \code{"opg"})}{Uses the Bartlett
#'     identity of the corresponding order, which expresses the expected
#'     derivative through expectations of \emph{products of lower-order}
#'     derivatives. Writing the identity as a sum over set partitions of the
#'     index set,
#'     \deqn{\sum_{\pi} \mathbb{E}\left[\prod_{B \in \pi} \ell_B\right] = 0,}
#'     the target term (the single-block partition) is obtained from all the
#'     others. At order 2 this is exactly the familiar \strong{outer product of
#'     gradients}, \eqn{\mathbb{E}[\ell_{ij}] = -\mathbb{E}[\ell_i \ell_j]};
#'     at order 3, \eqn{\mathbb{E}[\ell_{ijk}] = -\left(\mathbb{E}[\ell_{ij}\ell_k] +
#'     \mathbb{E}[\ell_{ik}\ell_j] + \mathbb{E}[\ell_{jk}\ell_i] +
#'     \mathbb{E}[\ell_i\ell_j\ell_k]\right)}, and so on.
#'     \emph{Fastest at order 2} (only the score is needed) and the only variant
#'     that stays valid when the log-likelihood is not differentiable in a
#'     parameter, where \eqn{\mathbb{E}[H]} degenerates but the score variance
#'     still gives the information (see \code{\link{laplace_distrib}}).
#'     Deterministic. At higher orders it needs several integrals, so it is
#'     usually slower than \code{"integrate"}, but it never requires the
#'     top-order derivative itself -- useful when only the lower orders are
#'     available in closed form.}
#'
#'   \item{\code{"integrate"}}{Integrates the observed derivative of that order
#'     directly against the density (numerical quadrature for continuous
#'     distributions, series summation for discrete ones). Deterministic and
#'     normally the most accurate when the observed derivative is available in
#'     closed form. Estimates \eqn{\mathbb{E}[\partial^k \ell]} literally, which
#'     for a non-regular model is \emph{not} the information.}
#'
#'   \item{\code{"mc"}}{Simulates \code{nsim} observations from the distribution
#'     and averages the observed derivative over them. The simplest and most
#'     robust option when quadrature struggles (heavy tails, awkward supports),
#'     but stochastic: the error decreases only as \eqn{1/\sqrt{n_{sim}}}, so
#'     results are not exactly reproducible unless the seed is fixed. Estimates
#'     the same quantity as \code{"integrate"}.
#'     \strong{Note:} the cost of this option is the cost of simulating from the
#'     distribution. A native \code{\link{distrib_rng}} is best, and failing that
#'     the default RNG falls back to \code{\link{rng_grou}}, which only needs the
#'     density. It becomes expensive only for a distribution that supplies a
#'     quantile function slow enough to make inverse transform sampling the
#'     bottleneck; prefer \code{"bartlett"} or \code{"integrate"} in that case.}
#' }
#'
#' \strong{Defaults.} \code{distrib_expected_hessian} defaults to
#' \code{"bartlett"}, because at order 2 it is both the cheapest (only first
#' derivatives) and the most broadly valid. \code{distrib_deriv3} and
#' \code{distrib_deriv4} default to \code{"integrate"}, since at those orders
#' direct integration of the available derivative is usually cheaper and more
#' accurate.
#'
#' @seealso \code{\link{distrib_expected_hessian}}, \code{\link{distrib_deriv3}},
#'   \code{\link{distrib_deriv4}}
NULL

# All set partitions of {1, ..., n}, as a list of lists of integer vectors.
set_partitions <- function(n) {
  if (n <= 1L) return(list(list(1L)))
  prev <- set_partitions(n - 1L)
  out <- vector("list", 0L)
  for (p in prev) {
    for (b in seq_along(p)) {
      q <- p
      q[[b]] <- c(q[[b]], n)
      out[[length(out) + 1L]] <- q
    }
    q <- p
    q[[length(q) + 1L]] <- n
    out[[length(out) + 1L]] <- q
  }
  out
}

# Observed derivative components of a given order, as a named list.
observed_deriv <- function(distrib, y, theta, order) {
  switch(as.character(order),
    "1" = distrib_gradient(distrib, y, theta),
    "2" = distrib_hessian(distrib, y, theta),
    "3" = distrib_deriv3(distrib, y, theta),
    "4" = distrib_deriv4(distrib, y, theta),
    stop("Unsupported derivative order: ", order, call. = FALSE)
  )
}

# --- strategy: numerical integration of the observed derivative -------------
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

# --- strategy: Monte Carlo average of the observed derivative ---------------
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

# --- strategy: Bartlett identity (order 2 == outer product of gradients) ----
expected_by_bartlett <- function(distrib, y, theta, order) {
  params <- distrib@params
  n <- length(y)
  nms <- if (order == 2L) hess_names(params) else deriv_names(params, order)

  # index tuple behind each component name
  idx_of <- lapply(nms, function(nm) match(strsplit(nm, "_", fixed = TRUE)[[1]], params))

  # every partition except the single-block one (that is the target term)
  parts <- Filter(function(pp) length(pp) > 1L, set_partitions(order))

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

# Dispatcher shared by the expected-derivative fallbacks.
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
