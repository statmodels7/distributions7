#' @include distrib.R generics.R
NULL

#' @title S7 Class for Fisher Scoring Specifications
#' @name FisherScoring
#'
#' @description
#' The class [fisher_scoring()] returns and [fit_distrib()] accepts through its
#' `method` argument. It holds the two things a Fisher scoring run needs that
#' an ordinary optimizer cannot carry: how the expected information is to be
#' obtained when the family has no closed form for it, and the stopping rule
#' and iteration budget the run is to use.
#'
#' This page documents the raw S7 constructor, which validates nothing. Build
#' one with [fisher_scoring()], which checks all four arguments.
#'
#' @param approx How the expectation is approximated when the family has no
#'   closed-form expected information: one of `"opg"`, `"bartlett"`,
#'   `"integrate"` or `"mc"`. Ignored by a family that has one, and
#'   [fit_distrib()] rejects any value other than the default `"opg"` in that
#'   case.
#' @param nsim Monte Carlo sample size, a single positive number, read only
#'   when `approx = "mc"` is in force.
#' @param criterion A stopping rule from \pkg{optimizers7}, or `NULL` to take
#'   [optimizers7::crit_grad()] at its default tolerance.
#' @param maxit An iteration limit, a single positive number, or `NULL` for the
#'   default of [optimizers7::newton()].
#'
#' @return An S7 object of class `FisherScoring` with those four properties.
#'
#' @examples
#' fs <- fisher_scoring(approx = "mc", nsim = 2000)
#' fs
#' c(approx = fs@approx, nsim = fs@nsim)
#'
#' # The defaults leave the stopping rule and the budget to the fit.
#' d <- fisher_scoring()
#' c(criterion = is.null(d@criterion), maxit = is.null(d@maxit))
#'
#' @seealso [fisher_scoring()], which builds and validates one;
#'   [fit_distrib()], which consumes it.
#' @export
FisherScoring <- S7::new_class("FisherScoring",
  properties = list(
    approx = S7::class_character,
    nsim = S7::class_numeric,
    criterion = S7::class_any,
    maxit = S7::class_any
  )
)

#' @title Fisher Scoring, With Its Own Settings
#'
#' @description
#' Returns a specification of Fisher scoring for [fit_distrib()]'s `method`
#' argument, carrying how the expected information is to be obtained when the
#' family has no closed form for it, and optionally a stopping rule and an
#' iteration limit of its own. It is [fit_distrib()]'s default method.
#'
#' Fisher scoring is Newton's method with the observed Hessian replaced by
#' minus the expected information, so it needs no implementation of its own and
#' is not an optimizer. What it does need, and an optimizer has no slot for, is
#' the statement of how that matrix is to be obtained, and that is what this
#' object holds.
#'
#' @details
#' # Values accepted by the method argument
#' [fit_distrib()] takes one argument saying how to optimize, and it takes an
#' optimizer of \pkg{optimizers7} or this:
#'
#' \tabular{ll}{
#'   `method = fisher_scoring()` \tab Newton's method with the
#'     **expected** information \cr
#'   `method = optimizers7::newton()` \tab Newton's method with the
#'     **observed** Hessian \cr
#'   `method = optimizers7::lbfgs()` \tab whatever that optimizer does
#' }
#'
#' # The step
#' Writing \eqn{s(\eta) = \partial \ell / \partial \eta} for the score on the
#' unconstrained scale and \eqn{\mathcal{I}(\eta) = \mathbb{E}[-\partial^{2}
#' \ell / \partial \eta \partial \eta']} for the expected information there,
#' \deqn{\eta^{(t+1)} = \eta^{(t)}
#'   + \mathcal{I}(\eta^{(t)})^{-1} s(\eta^{(t)}).}
#' Under the diagonal reparametrization \eqn{\theta_i = h_i(\eta_i)} the
#' information transforms by congruence,
#' \eqn{\mathcal{I}(\eta) = D\,\mathcal{I}(\theta)\,D} with
#' \eqn{D = \operatorname{diag}(h_i'(\eta_i))}, the first-order term of the
#' chain rule vanishing in expectation because the score has mean zero. The
#' matrix inverted is therefore positive definite wherever
#' \eqn{\mathcal{I}(\theta)} is, and the step is a direction of ascent. A fit
#' therefore succeeds on a family whose observed Hessian is indefinite, the
#' Laplace being the standard case.
#'
#' # Strategies for the expected information
#' `approx` is read only where the family has no closed-form expected
#' information, as for a family defined by its density alone or a truncated
#' family. The strategies differ in cost and in variance: `"opg"` costs one
#' gradient, `"bartlett"` and `"integrate"` are deterministic and dearer, and
#' `"mc"` adds simulation noise that decreases with `nsim`.
#'
#' A family that supplies the expectation in closed form ignores `approx`, and
#' [fit_distrib()] signals an error when `approx` differs from the default
#' `"opg"` for such a family, so that a fit never reports a strategy that was
#' not used.
#'
#' @param approx How the expectation is approximated when the family has no
#'   closed-form expected information. `"opg"` (the default) is the outer
#'   product of the observed scores, which costs one gradient; `"bartlett"`
#'   evaluates the right-hand side of the second Bartlett identity, a sum over
#'   the support for a discrete family and a quadrature for a continuous one;
#'   `"integrate"` is quadrature of the observed information; `"mc"` is Monte
#'   Carlo. Matched by [base::match.arg()], so
#'   anything else signals an error listing the four. See
#'   [expected_derivative_methods()].
#' @param nsim Monte Carlo sample size, a single positive finite number. Read
#'   only when `approx = "mc"`. Defaults to 10000. A vector, a non-finite
#'   value or a number below 1 signals an error.
#' @param criterion A stopping rule from \pkg{optimizers7}, or `NULL` (the
#'   default), in which case [fit_distrib()] stops on
#'   [optimizers7::crit_grad()] at its default tolerance. An object that is not
#'   a `criterion` signals an error.
#' @param maxit An iteration limit, a single number at least 1, or `NULL` (the
#'   default) for the limit of [optimizers7::newton()]. When the limit is reached the fit reports
#'   `converged = FALSE` and keeps the point.
#'
#' @return An S7 object of class [FisherScoring()], with properties `approx`,
#'   `nsim`, `criterion` and `maxit`.
#'
#' @examples
#' set.seed(1)
#' d <- gaussian1_distrib()
#' y <- distrib_rng(d, 100, list(mu = 1, sigma = 2))
#'
#' # The default, and the same thing said explicitly.
#' coef(fit_distrib(d, y))
#' coef(fit_distrib(d, y, method = fisher_scoring()))
#'
#' # The congruence the step rests on: the information on the link scale is
#' # D I(theta) D, with D the diagonal of the inverse link's derivative.
#' th <- list(mu = 1, sigma = 2)
#' It <- distrib_expected_hessian(d, y, th)
#' Ie <- distrib_expected_hessian(d, y, th, scale = "link")
#' c(theta = sum(It$sigma_sigma), link = sum(Ie$sigma_sigma))   # ratio 2^2
#'
#' # 'approx' is read only where the family has no closed form (see the
#' # examples of expected_derivative_methods). On a family that has one, a
#' # value other than the default is rejected.
#' try(fit_distrib(d, y, method = fisher_scoring(approx = "mc")))
#'
#' # The stopping rule and the budget belong to the method.
#' fit <- fit_distrib(d, y, method = fisher_scoring(maxit = 1))
#' c(iterations = fit@iterations, converged = fit@converged)
#'
#' @seealso [fit_distrib()], which consumes this;
#'   [FisherScoring()] for the class;
#'   [expected_derivative_methods()] for what each `approx` computes;
#'   [distrib_expected_hessian()] for the matrix itself;
#'   [optimizers7::newton()], whose iteration limit stands where this sets
#'   none.
#' @export
fisher_scoring <- function(approx = c("opg", "bartlett", "integrate", "mc"),
                           nsim = 10000, criterion = NULL, maxit = NULL) {
  approx <- match.arg(approx)
  if (!is.numeric(nsim) || length(nsim) != 1L || !is.finite(nsim) || nsim < 1) {
    stop("'nsim' must be a single positive number.", call. = FALSE)
  }
  if (!is.null(criterion) &&
      !S7::S7_inherits(criterion, optimizers7::criterion)) {
    stop("'criterion' must be an optimizers7 stopping rule, or NULL.",
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
#'
#' @description
#' Writes the strategy for the expected information, with the Monte Carlo
#' sample size when that is the strategy, and then any setting that overrides
#' the fit's: a stopping rule set here is reported as set, and an iteration
#' limit by its value. The line about the strategy carries the reminder that a
#' family with a closed-form expectation ignores it.
#'
#' @param x A [FisherScoring()] object.
#' @param ... Unused, accepted for compatibility with [base::print()].
#'
#' @return `x`, invisibly. Called for the output it writes.
#'
#' @examples
#' fisher_scoring()
#' fisher_scoring(approx = "mc", nsim = 2000)
#' fisher_scoring(criterion = optimizers7::crit_grad(1e-9), maxit = 50)
#'
#' @seealso [fisher_scoring()] for what each line reports.
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
