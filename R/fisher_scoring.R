#' @include distrib.R generics.R
NULL

#' @title Fisher Scoring as an Object
#' @name FisherScoring
#' @description
#' The S7 class of Fisher scoring specifications, returned by
#' \code{\link{fisher_scoring}} and passed to \code{\link{fit_distrib}} through
#' its \code{method} argument.
#' @param approx How the expectation is approximated when the distribution has
#'   no closed-form expected information.
#' @param nsim Monte Carlo sample size, used when \code{approx = "mc"}.
#' @param criterion The stopping rule, or \code{NULL} to take the fit's.
#' @param maxit The iteration limit, or \code{NULL} to take the fit's.
#' @return An object of class \code{FisherScoring}.
#' @seealso \code{\link{fisher_scoring}}
#' @examples
#' fs <- fisher_scoring(approx = "mc", nsim = 2000)
#' S7::S7_inherits(fs, FisherScoring)
#' fs@approx
#' @export
FisherScoring <- S7::new_class("FisherScoring",
  properties = list(
    approx = S7::class_character,
    nsim = S7::class_numeric,
    criterion = S7::class_any,
    maxit = S7::class_any
  )
)

#' Fisher Scoring, With Its Own Settings
#'
#' @description
#' Returns a specification of Fisher scoring for \code{\link{fit_distrib}},
#' carrying how the expected information is to be obtained when the
#' distribution has no closed form for it.
#'
#' @details
#' \code{fit_distrib()} takes one argument saying how to optimise, and it takes
#' either an optimiser of \pkg{optimizers7} or this. The two are the same kind
#' of thing said in the same place:
#'
#' \tabular{ll}{
#'   \code{method = fisher_scoring()} \tab Newton's method with the
#'     \strong{expected} information \cr
#'   \code{method = optimizers7::newton()} \tab Newton's method with the
#'     \strong{observed} Hessian \cr
#'   \code{method = optimizers7::lbfgs()} \tab whatever that optimiser does
#' }
#'
#' Fisher scoring is not a separate algorithm, which is why it has no
#' implementation of its own: it is a Newton step with one matrix replaced by
#' another. What it does need, and an optimiser cannot carry, is a statement of
#' how that matrix is to be obtained when the family does not supply it in
#' closed form --- and that is what this object holds. A family that does
#' supply one ignores \code{approx} entirely, and \code{fit_distrib()} refuses
#' the argument in that case rather than accepting something it will not use.
#'
#' Fisher scoring is Newton's method with one matrix replaced, so how the run
#' stops and how long it may take are set here, as they would be on any other
#' optimiser. Both default to \code{NULL}, meaning the defaults of
#' \code{\link[optimizers7]{newton}} and \code{\link[optimizers7]{crit_grad}}
#' stand.
#'
#' @param approx How the expectation is approximated when the distribution has
#'   no closed-form expected information: \code{"bartlett"} (the default, the
#'   outer product of the score, equivalently \code{"opg"}), \code{"integrate"}
#'   for quadrature of the observed information, or \code{"mc"} for Monte
#'   Carlo. See \code{\link{expected_derivative_methods}}.
#' @param nsim Monte Carlo sample size, used when \code{approx = "mc"}.
#' @param criterion A stopping rule from \pkg{optimizers7}, or \code{NULL} for
#'   the default of \code{\link[optimizers7]{newton}}.
#' @param maxit An iteration limit, or \code{NULL} for the same default.
#'
#' @return An object of class \code{\link{FisherScoring}}.
#'
#' @seealso \code{\link{fit_distrib}}, \code{\link{expected_derivative_methods}}
#'
#' @examples
#' set.seed(1)
#' d <- gaussian_distrib()
#' y <- distrib_rng(d, 200, list(mu = 1, sigma = 2))
#'
#' # the default, and the same thing said explicitly
#' coef(fit_distrib(d, y))
#' coef(fit_distrib(d, y, method = fisher_scoring()))
#'
#' # A family whose expected information has no closed form takes a strategy.
#' # The same argument on a family that HAS one is refused rather than
#' # silently ignored.
#' sn <- skewnormal_distrib()
#' set.seed(2)
#' ys <- distrib_rng(sn, 300, list(mu = 0, sigma = 1, alpha = 3))
#' coef(fit_distrib(sn, ys, method = fisher_scoring(approx = "opg")))
#'
#' try(fit_distrib(d, y, method = fisher_scoring(approx = "mc")))
#'
#' @export
fisher_scoring <- function(approx = c("bartlett", "integrate", "mc", "opg"),
                           nsim = 10000, criterion = NULL, maxit = NULL) {
  approx <- match.arg(approx)
  if (!is.numeric(nsim) || length(nsim) != 1L || !is.finite(nsim) || nsim < 1) {
    stop("'nsim' must be a single positive number.", call. = FALSE)
  }
  if (!is.null(criterion) &&
      !S7::S7_inherits(criterion, optimizers7::criterion)) {
    stop("'criterion' must be a optimizers7 stopping rule, or NULL.",
         call. = FALSE)
  }
  if (!is.null(maxit) &&
      (!is.numeric(maxit) || length(maxit) != 1L || maxit < 1)) {
    stop("'maxit' must be a single positive number, or NULL.", call. = FALSE)
  }
  FisherScoring(
    approx = approx, nsim = nsim,
    criterion = criterion, maxit = maxit
  )
}

#' @title Print a Fisher Scoring Specification
#' @name print.FisherScoring
#' @description Reports the strategy and any settings that override the fit's.
#' @param x A \code{\link{FisherScoring}} object.
#' @param ... Unused.
#' @return \code{x}, invisibly.
#' @keywords internal
S7::method(print, FisherScoring) <- function(x, ...) {
  cat("Fisher scoring\n")
  cat("  expected information: ", x@approx,
      if (identical(x@approx, "mc")) paste0(" (nsim = ", x@nsim, ")") else "",
      "  [ignored when the family has a closed form]\n", sep = "")
  if (!is.null(x@criterion)) cat("  criterion: set here\n")
  if (!is.null(x@maxit)) cat("  maxit: ", x@maxit, "\n", sep = "")
  invisible(x)
}
