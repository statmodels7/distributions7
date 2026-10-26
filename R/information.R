#' @include expect_approx.R
NULL

#' @title S7 Classes for the Information Matrix
#' @name Information
#' @aliases ObservedInfo ExpectedInfo
#'
#' @description
#' `Information` is the abstract parent of the objects that select which
#' matrix is taken as the information of a likelihood: the observed
#' information, which is the negative Hessian of the log-likelihood at the
#' data, or the expected (Fisher) information. Its two subclasses are
#' `ObservedInfo` and `ExpectedInfo`. They are built by [observed_info()] and
#' [expected_info()], which are the functions to call; the classes are
#' exported so that `S7::S7_inherits()` can test an object, and so that a
#' package can define a subclass of `Information` for a further choice of
#' matrix. [information_objects] compares the two.
#'
#' `ObservedInfo` has no properties. `ExpectedInfo` has one, `approx`, an
#' object that inherits from [ExpectApprox] and selects how the expected
#' information is computed for a family that has no closed form for it. The raw
#' constructors do not validate their properties; the functions
#' [observed_info()] and [expected_info()] do.
#'
#' @param approx An object inheriting from [ExpectApprox], for `ExpectedInfo`.
#'
#' @return An S7 object of the class named, which inherits from `Information`.
#'
#' @examples
#' a <- expected_info(approx = expect_mc(nsim = 500))
#' S7::S7_inherits(a, Information)
#' S7::S7_inherits(a, ExpectedInfo)
#' a@approx
#'
#' @seealso [information_objects] for the comparison of the two constructors;
#'   [observed_info()] and [expected_info()] for the constructors;
#'   [vcov.distrib_fit()] for a function that takes the object.
#' @export
Information <- S7::new_class("Information", abstract = TRUE)

#' @rdname Information
#' @export
ObservedInfo <- S7::new_class("ObservedInfo", parent = Information)

#' @rdname Information
#' @export
ExpectedInfo <- S7::new_class("ExpectedInfo", parent = Information,
  properties = list(approx = ExpectApprox)
)


#' @title Information Matrices
#' @name information_objects
#'
#' @description
#' The information of a likelihood is the matrix of second derivatives that
#' measures how sharply the log-likelihood \eqn{\ell(\eta)} peaks at the
#' estimate \eqn{\hat\eta} of the parameters on the link scale. Its inverse is
#' the estimated variance matrix of the estimates, so the choice of the matrix
#' is the choice of the standard errors. Two objects select it, [observed_info()]
#' and [expected_info()], and both inherit from [Information]. A function that
#' takes such an object, such as the `vcov()` method of a fit (see
#' [vcov.distrib_fit()]), is told the matrix by the object and not by a
#' character string.
#'
#' @details
#' \describe{
#'   \item{[observed_info()]}{The negative Hessian of the log-likelihood at the
#'     data,
#'     \deqn{-\frac{\partial^2 \ell}{\partial\eta\,\partial\eta^\top}\bigg|_{\hat\eta}.}
#'     It is computed from the second derivatives that the family supplies in
#'     closed form, or from their numerical fallbacks when the family
#'     supplies only a density. It has no approximation and takes no
#'     argument. It is the curvature of the likelihood in the sample at hand,
#'     and every family has it.}
#'   \item{[expected_info()]}{The Fisher information,
#'     \deqn{\mathbb{E}\left[-\frac{\partial^2 \ell}{\partial\eta\,\partial\eta^\top}
#'     \right]_{\hat\eta},}
#'     the average of the observed information over samples drawn from the
#'     model at the estimates. Where the family writes the expectation in
#'     closed form, that form is used and `approx` is not read. Where it does
#'     not (see [expected_hessian_exact()]), the expectation is approximated by
#'     the object in `approx`, one of [expect_opg()], [expect_bartlett()],
#'     [expect_integrate()] and [expect_mc()]; [expected_derivative_methods]
#'     compares the four.}
#' }
#'
#' The two matrices agree asymptotically and differ in a finite sample. The
#' expected information is a function of the parameters alone and is positive
#' definite wherever the family is regular. The observed information depends
#' on the data, and it can be indefinite away from the estimate or at an
#' outlier. Standard errors from the observed information describe the
#' likelihood that was actually observed, which is the usual recommendation
#' (Efron and Hinkley, 1978); for this reason [observed_info()] is the default
#' of [vcov.distrib_fit()]. For a family that has no closed-form expected
#' information, [expect_opg()] is a noisier estimate of the information than the
#' observed matrix.
#'
#' @references
#' Efron, B. and Hinkley, D. V. (1978). Assessing the accuracy of the maximum
#' likelihood estimator: observed versus expected Fisher information.
#' *Biometrika*, 65, 457-483.
#'
#' @return No value. The topic describes the two objects, of classes
#'   [ObservedInfo][Information] and [ExpectedInfo][Information], that
#'   [observed_info()] and [expected_info()] return.
#'
#' @seealso [observed_info()] and [expected_info()] for the constructors;
#'   [Information] for the classes; [vcov.distrib_fit()], [summary.distrib_fit()]
#'   and [confint.distrib_fit()] for the functions that take the object;
#'   [expected_derivative_methods] for the approximations of the expected
#'   information.
#'
#' @examples
#' observed_info()
#' expected_info()
#' expected_info(approx = expect_bartlett())
#'
#' # The standard errors of one fit under the two matrices.
#' set.seed(1)
#' d <- gaussian1_distrib()
#' fit <- fit_distrib(d, distrib_rng(d, 300, list(mu = 1, sigma = 2)))
#' sqrt(diag(vcov(fit, information = observed_info())))
#' sqrt(diag(vcov(fit, information = expected_info())))
NULL


#' @title The Observed Information
#' @name observed_info
#'
#' @description
#' Selects the observed information, the negative Hessian of the
#' log-likelihood at the data,
#' \deqn{-\frac{\partial^2 \ell}{\partial\eta\,\partial\eta^\top},}
#' evaluated at the estimates. The matrix is computed from the second
#' derivatives of the family, which are analytic where the family writes them
#' and numerical otherwise (see [distrib_hessian()]). It has no approximation
#' and the function takes no argument.
#'
#' It is the default of the `information` argument of [vcov.distrib_fit()],
#' [summary.distrib_fit()], [confint.distrib_fit()] and [print.distrib_fit()],
#' so the standard errors that a fit reports are the inverse of this matrix
#' unless another object is passed. [information_objects] compares it with
#' [expected_info()].
#'
#' @return An S7 object of class `ObservedInfo`, which inherits from
#'   [Information], with no properties.
#'
#' @examples
#' observed_info()
#'
#' set.seed(1)
#' d <- gamma1_distrib()
#' fit <- fit_distrib(d, distrib_rng(d, 300, list(mu = 3, phi = 0.6)))
#' vcov(fit, information = observed_info())
#'
#' @seealso [expected_info()] for the other object; [information_objects] for
#'   the comparison; [vcov.distrib_fit()] for the function that takes the
#'   object.
#' @export
observed_info <- function() ObservedInfo()


#' @title The Expected (Fisher) Information
#' @name expected_info
#'
#' @description
#' Selects the expected information, the Fisher information
#' \deqn{\mathbb{E}\left[-\frac{\partial^2 \ell}{\partial\eta\,\partial\eta^\top}
#' \right],}
#' evaluated at the estimates. Where the family writes the expected information
#' out, the closed form is used and `approx` is not read. Where it does not
#' (see [expected_hessian_exact()]), the expectation is approximated by
#' `approx`.
#'
#' @param approx An object that selects the approximation of the expectation,
#'   read only for a family with no closed-form expected information:
#'   - [expect_opg()], the default: the outer product of the observed scores,
#'     which costs one gradient;
#'   - [expect_bartlett()]: the right-hand side of the second Bartlett
#'     identity, a sum over the support for a discrete family and a quadrature
#'     for a continuous one;
#'   - [expect_integrate()]: quadrature of the observed information against
#'     the density;
#'   - [expect_mc()]: Monte Carlo, with its own sample size `nsim`.
#'
#'   [expected_derivative_methods] compares the four. Any other value signals
#'   an error.
#'
#' @return An S7 object of class `ExpectedInfo`, which inherits from
#'   [Information], with the property `approx`.
#'
#' @examples
#' expected_info()
#' expected_info(approx = expect_mc(nsim = 2000))
#'
#' set.seed(1)
#' d <- gamma1_distrib()
#' fit <- fit_distrib(d, distrib_rng(d, 300, list(mu = 3, phi = 0.6)))
#' vcov(fit, information = expected_info())
#'
#' @seealso [observed_info()] for the other object; [information_objects] for
#'   the comparison; [expected_derivative_methods] for the approximations;
#'   [vcov.distrib_fit()] for the function that takes the object.
#' @export
expected_info <- function(approx = expect_opg()) {
  check_approx(approx)
  ExpectedInfo(approx = approx)
}


#' Check That an Argument Is an Information Object
#'
#' @description
#' Returns `x` invisibly when it is an [observed_info()] or an
#' [expected_info()] object and signals an error otherwise.
#'
#' @param x The value to check.
#' @param arg The argument name used in the message.
#'
#' @return `x`, invisibly.
#'
#' @seealso [information_objects]
#' @keywords internal
check_information <- function(x, arg = "information") {
  if (S7::S7_inherits(x, ObservedInfo) || S7::S7_inherits(x, ExpectedInfo)) {
    return(invisible(x))
  }
  stop(sprintf("'%s' must be observed_info() or expected_info().", arg),
       call. = FALSE)
}


#' @title Print an Information Object
#' @name print.Information
#'
#' @description
#' Writes one line naming the matrix and, for [expected_info()], the
#' approximation of the expectation.
#'
#' @param x An object inheriting from [Information].
#' @param ... Unused, accepted for compatibility with [base::print()].
#'
#' @return `x`, invisibly. The function is called for the output that it writes.
#'
#' @examples
#' observed_info()
#' expected_info(approx = expect_mc(nsim = 500))
#'
#' @seealso [information_objects]
#' @keywords internal
S7::method(print, Information) <- function(x, ...) {
  if (S7::S7_inherits(x, ExpectedInfo)) {
    lab <- sub("^Expectation by ", "expectation by ", approx_label(x@approx))
    cat(sprintf("Expected information (%s)\n", lab))
  } else if (S7::S7_inherits(x, ObservedInfo)) {
    cat("Observed information\n")
  } else {
    cat("Information of class", class(x)[[1L]], "\n")
  }
  invisible(x)
}
