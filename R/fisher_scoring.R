#' @include distrib.R generics.R expect_approx.R
NULL

#' @title S7 Class for Fisher Scoring Specifications
#' @name FisherScoring
#'
#' @description
#' This class is the one that [fisher_scoring()] returns and that
#' [fit_distrib()] accepts through its `method` argument. It holds the two
#' things that a Fisher scoring run needs and that an ordinary optimizer cannot
#' carry: the approximation of the expected information, an object of class
#' [ExpectApprox], used when the family has no closed form for it, and the
#' stopping rule and iteration budget that the run uses.
#'
#' This page documents the raw S7 constructor, which validates nothing. Build
#' one with [fisher_scoring()], which checks all four arguments.
#'
#' @param approx An object inheriting from [ExpectApprox] that gives the
#'   approximation of the expected information when the family has no
#'   closed-form expected information. It is ignored by a family that has one,
#'   and [fit_distrib()] rejects any approximation other than [expect_opg()]
#'   in that case.
#' @param criterion A stopping rule from \pkg{optimizers7}, or `NULL` to take
#'   [optimizers7::crit_grad()] at its default tolerance.
#' @param maxit An iteration limit, a single positive number, or `NULL` for the
#'   default of [optimizers7::newton()].
#'
#' @return An S7 object of class `FisherScoring` with the properties `approx`,
#'   `criterion` and `maxit`.
#'
#' @examples
#' fs <- fisher_scoring(approx = expect_mc(nsim = 2000))
#' fs
#' fs@approx@nsim
#'
#' # The defaults leave the stopping rule and the budget to the fit.
#' d <- fisher_scoring()
#' c(criterion = is.null(d@criterion), maxit = is.null(d@maxit))
#'
#' @seealso [fisher_scoring()], which builds and validates one;
#'   [fit_distrib()], which consumes it; [ExpectApprox] for the class of the
#'   `approx` property.
#' @export
FisherScoring <- S7::new_class("FisherScoring",
  properties = list(
    approx = ExpectApprox,
    criterion = S7::class_any,
    maxit = S7::class_any
  )
)

#' @title Fisher Scoring, With Its Own Settings
#'
#' @description
#' Returns a specification of Fisher scoring for [fit_distrib()]'s `method`
#' argument. The specification records how the expected information is
#' obtained when the family has no closed form for it, and optionally a
#' stopping rule and an iteration limit of its own. It is [fit_distrib()]'s
#' default method.
#'
#' Fisher scoring is Newton's method in which the observed Hessian is replaced
#' by minus the expected information, so it needs no implementation of its own
#' and is not an optimizer. It needs one thing for which an optimizer has no
#' slot, the statement of how that matrix is obtained, and this object holds
#' that statement.
#'
#' @details
#' # Values accepted by the method argument
#' [fit_distrib()] has one argument that selects how to optimize, and that
#' argument accepts an optimizer of \pkg{optimizers7} or this object:
#'
#' \tabular{ll}{
#'   `method = fisher_scoring()` \tab Newton's method with the
#'     expected information \cr
#'   `method = optimizers7::newton()` \tab Newton's method with the
#'     observed Hessian \cr
#'   `method = optimizers7::lbfgs()` \tab whatever that optimizer does
#' }
#'
#' # The step
#' Writing \eqn{s(\eta) = \partial \ell / \partial \eta} for the score on the
#' unconstrained scale and \eqn{\mathcal{I}(\eta) = \mathbb{E}[-\partial^{2}
#' \ell / \partial \eta \partial \eta']} for the expected information there,
#' the step is
#' \deqn{\eta^{(t+1)} = \eta^{(t)}
#'   + \mathcal{I}(\eta^{(t)})^{-1} s(\eta^{(t)}).}
#' Under the diagonal reparametrization \eqn{\theta_i = h_i(\eta_i)} the
#' information transforms by congruence,
#' \eqn{\mathcal{I}(\eta) = D\,\mathcal{I}(\theta)\,D} with
#' \eqn{D = \operatorname{diag}(h_i'(\eta_i))}, and the first-order term of the
#' chain rule vanishes in expectation because the score has mean zero. The
#' matrix inverted is therefore positive definite wherever
#' \eqn{\mathcal{I}(\theta)} is, and the step is a direction of ascent. A fit
#' therefore succeeds on a family whose observed Hessian is indefinite; the
#' Laplace is the standard case.
#'
#' # Strategies for the expected information
#' `approx` is read only where the family has no closed-form expected
#' information, as for a family defined by its density alone or a truncated
#' family. The approximations differ in cost and in variance:
#' [expect_opg()] costs one gradient, [expect_bartlett()] and
#' [expect_integrate()] are deterministic and dearer, and [expect_mc()] adds
#' simulation noise that decreases with its `nsim`.
#'
#' A family that supplies the expectation in closed form ignores `approx`, and
#' [fit_distrib()] signals an error when `approx` differs from the default
#' [expect_opg()] for such a family, so that a fit never reports an
#' approximation that was not used.
#'
#' @param approx An object that gives the approximation of the expected
#'   information when the family has no closed form for it, one of:
#'   - [expect_opg()], the default: the outer product of the observed scores,
#'     which costs one gradient;
#'   - [expect_bartlett()]: the right-hand side of the second Bartlett
#'     identity, a sum over the support for a discrete family and a quadrature
#'     for a continuous one;
#'   - [expect_integrate()]: quadrature of the observed information;
#'   - [expect_mc()]: Monte Carlo, with the sample size `nsim` set in that
#'     constructor.
#'
#'   Any other value, a character string included, signals an error.
#'   [expected_derivative_methods] compares the four.
#' @param criterion A stopping rule from \pkg{optimizers7}, or `NULL` (the
#'   default), in which case [fit_distrib()] stops on
#'   [optimizers7::crit_grad()] at its default tolerance. An object that is not
#'   a `criterion` signals an error.
#' @param maxit An iteration limit, a single number at least 1, or `NULL` (the
#'   default) for the limit of [optimizers7::newton()]. When the limit is reached the fit reports
#'   `converged = FALSE` and keeps the point.
#'
#' @return An S7 object of class [FisherScoring()], with properties `approx`
#'   (an object inheriting from [ExpectApprox]), `criterion` and `maxit`.
#'
#' @examples
#' set.seed(1)
#' d <- gaussian1_distrib()
#' y <- distrib_rng(d, 100, list(mu = 1, sigma = 2))
#'
#' # The default and the explicit specification give the same fit.
#' coef(fit_distrib(d, y))
#' coef(fit_distrib(d, y, method = fisher_scoring()))
#'
#' # The step rests on a congruence: the information on the link scale is
#' # D I(theta) D, with D the diagonal of the inverse link's derivative.
#' th <- list(mu = 1, sigma = 2)
#' It <- distrib_expected_hessian(d, y, th)
#' Ie <- distrib_expected_hessian(d, y, th, scale = "link")
#' c(theta = sum(It$sigma_sigma), link = sum(Ie$sigma_sigma))   # ratio 2^2
#'
#' # 'approx' is read only where the family has no closed form (see the
#' # examples of expected_derivative_methods). On a family that has one, a
#' # value other than the default is rejected.
#' try(fit_distrib(d, y, method = fisher_scoring(approx = expect_mc())))
#'
#' # The stopping rule and the budget belong to the method.
#' fit <- fit_distrib(d, y, method = fisher_scoring(maxit = 1))
#' c(iterations = fit@iterations, converged = fit@converged)
#'
#' @seealso [fit_distrib()], which consumes this;
#'   [FisherScoring()] for the class;
#'   [expected_derivative_methods] for what each `approx` computes;
#'   [distrib_expected_hessian()] for the matrix itself;
#'   [optimizers7::newton()], whose iteration limit stands where this sets
#'   none.
#' @export
fisher_scoring <- function(approx = expect_opg(), criterion = NULL,
                           maxit = NULL) {
  check_approx(approx)
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
    approx = approx,
    criterion = criterion, maxit = maxit
  )
}

#' @title Print a Fisher Scoring Specification
#' @name print.FisherScoring
#'
#' @description
#' Writes the approximation of the expected information, with the Monte Carlo
#' sample size when that is the approximation, and then any setting that overrides
#' the fit's: a stopping rule set here is reported as set, and an iteration
#' limit by its value. The line about the approximation carries the reminder that a
#' family with a closed-form expectation ignores it.
#'
#' @param x A [FisherScoring()] object.
#' @param ... Unused, accepted for compatibility with [base::print()].
#'
#' @return `x`, invisibly. The function is called for the output that it writes.
#'
#' @examples
#' fisher_scoring()
#' fisher_scoring(approx = expect_mc(nsim = 2000))
#' fisher_scoring(criterion = optimizers7::crit_grad(1e-9), maxit = 50)
#'
#' @seealso [fisher_scoring()] for what each line reports.
#' @keywords internal
S7::method(print, FisherScoring) <- function(x, ...) {
  cat("Fisher scoring\n")
  cat("  expected information: ", sub("^Expectation ", "", approx_label(x@approx)),
      "  [ignored when the family has a closed form]\n", sep = "")
  if (!is.null(x@criterion)) cat("  criterion: set here\n")
  if (!is.null(x@maxit)) cat("  maxit: ", x@maxit, "\n", sep = "")
  invisible(x)
}
