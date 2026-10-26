#' @include expect_method.R
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
#' `ObservedInfo` has no properties. `ExpectedInfo` has one, `method`, which is
#' `NULL` or an object that inherits from [ExpectMethod] and selects how the
#' expected information is obtained. The raw constructors do not validate their
#' properties; the functions [observed_info()] and [expected_info()] do.
#'
#' @param method `NULL` or an object inheriting from [ExpectMethod], for
#'   `ExpectedInfo`.
#'
#' @return An S7 object of the class named, which inherits from `Information`.
#'
#' @examples
#' a <- expected_info(method = expect_mc(nsim = 500))
#' S7::S7_inherits(a, Information)
#' S7::S7_inherits(a, ExpectedInfo)
#' a@method
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
  properties = list(method = S7::new_union(NULL, ExpectMethod))
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
#' character string. The argument `information` of those functions is `NULL` by
#' default, which selects the matrix described under "The default" below.
#'
#' @details
#' \describe{
#'   \item{[observed_info()]}{The negative Hessian of the log-likelihood at the
#'     data,
#'     \deqn{-\frac{\partial^2 \ell}{\partial\eta\,\partial\eta^\top}\bigg|_{\hat\eta}.}
#'     It is computed from the second derivatives that the family supplies in
#'     closed form, or from their numerical fallbacks when the family
#'     supplies only a density. It takes no argument. It is the curvature of the
#'     likelihood in the sample at hand, and every family has it. It describes
#'     the uncertainty only where the log-likelihood is twice differentiable in
#'     every parameter at the estimate.}
#'   \item{[expected_info()]}{The Fisher information,
#'     \deqn{\mathbb{E}\left[-\frac{\partial^2 \ell}{\partial\eta\,\partial\eta^\top}
#'     \right]_{\hat\eta},}
#'     the average of the observed information over samples drawn from the
#'     model at the estimates. With `method = NULL`, the default, the closed
#'     form of the family is used where the family writes the expectation out
#'     (see [expected_hessian_exact()]), and [expect_opg()] otherwise. An object
#'     in `method` is used for every family, one of [expect_opg()],
#'     [expect_integrate()] and [expect_mc()]; [expected_derivative_methods]
#'     compares the three.}
#' }
#'
#' The two matrices agree asymptotically and differ in a finite sample. The
#' expected information is a function of the parameters alone and is positive
#' definite wherever the family is regular. The observed information depends
#' on the data, and it can be indefinite away from the estimate or at an
#' outlier. Standard errors from the observed information describe the
#' likelihood that was actually observed, which is the usual recommendation
#' (Efron and Hinkley, 1978). For a family that has no closed-form expected
#' information, [expect_opg()] is a noisier estimate of the information than the
#' observed matrix.
#'
#' # The default
#' The argument `information` of [vcov.distrib_fit()], [confint.distrib_fit()],
#' [print.distrib_fit()] and [mv_summary()] is `NULL`. It resolves to
#' [observed_info()] when every parameter of the family has
#' [params_order()] of at least 2, that is, when the log-likelihood is twice
#' differentiable in every parameter, and to [expected_info()] with
#' `method = NULL` otherwise, for the whole matrix. A family with a kink in a
#' parameter (the Laplace, the elastic net) has a singular or undefined
#' observed information in that parameter, and the closed-form expected
#' information gives finite standard errors. A fit states which information it
#' used whenever it is not the observed one.
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
#'   [Information] for the classes; [vcov.distrib_fit()],
#'   [confint.distrib_fit()] and [mv_summary()] for the functions that take the
#'   object;
#'   [expected_derivative_methods] for the methods of the expected
#'   information.
#'
#' @examples
#' observed_info()
#' expected_info()
#' expected_info(method = expect_integrate(identity = "score"))
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
#' and numerical otherwise (see [distrib_hessian()]). The function takes no
#' argument.
#'
#' It is the default of the `information` argument of [vcov.distrib_fit()],
#' [confint.distrib_fit()], [print.distrib_fit()] and [mv_summary()] for a family whose log-likelihood is twice differentiable in every
#' parameter (see [information_objects]). [information_objects] compares it with
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
#' evaluated at the estimates. With `method = NULL`, the default, the closed
#' form is used where the family writes the expected information out (see
#' [expected_hessian_exact()]), and [expect_opg()] otherwise. An object in
#' `method` is used for every family, also for one that has the closed form.
#'
#' @param method `NULL` (the default) or an object that selects how the
#'   expectation is obtained:
#'   - `NULL`: the closed form of the family where it has one, and
#'     [expect_opg()] otherwise;
#'   - [expect_opg()]: the outer product of the observed scores, which costs
#'     one gradient;
#'   - [expect_integrate()]: quadrature (or summation over the support) of the
#'     observed information against the density with `identity = "hessian"`,
#'     or of minus the outer product of the scores with `identity = "score"`;
#'   - [expect_mc()]: Monte Carlo, with its own sample size `nsim` and
#'     `identity`.
#'
#'   [expected_derivative_methods] compares the methods. Any other value
#'   signals an error.
#'
#' @return An S7 object of class `ExpectedInfo`, which inherits from
#'   [Information], with the property `method`.
#'
#' @examples
#' expected_info()
#' expected_info(method = expect_mc(nsim = 2000))
#'
#' set.seed(1)
#' d <- gamma1_distrib()
#' fit <- fit_distrib(d, distrib_rng(d, 300, list(mu = 3, phi = 0.6)))
#' vcov(fit, information = expected_info())
#'
#' @seealso [observed_info()] for the other object; [information_objects] for
#'   the comparison; [expected_derivative_methods] for the methods;
#'   [vcov.distrib_fit()] for the function that takes the object.
#' @export
expected_info <- function(method = NULL) {
  check_method(method)
  ExpectedInfo(method = method)
}


#' Check That an Argument Is an Information Object
#'
#' @description
#' Returns `x` invisibly when it is `NULL`, an [observed_info()] or an
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
  if (is.null(x) || S7::S7_inherits(x, ObservedInfo) ||
      S7::S7_inherits(x, ExpectedInfo)) {
    return(invisible(x))
  }
  stop(sprintf("'%s' must be NULL, observed_info() or expected_info().", arg),
       call. = FALSE)
}


#' The Default Information of a Family
#'
#' @description
#' Returns the information object that a function takes when its `information`
#' argument is `NULL`: [observed_info()] when every parameter of the family
#' has [params_order()] of at least 2, that is, when the log-likelihood is
#' twice differentiable in every parameter, and [expected_info()] with
#' `method = NULL` otherwise, for the whole matrix. A family with a kink in a
#' parameter (the Laplace, the elastic net) has a singular or undefined observed
#' information in that parameter, and the closed-form expected information gives
#' finite standard errors. A parameter whose order is not established (`NA`)
#' does not satisfy the condition. The rule is the one of
#' [vcov.distrib_fit()], and it is available to packages that take an
#' `information` argument with the same default.
#'
#' @param distrib An object inheriting from `distrib`.
#'
#' @return An [Information] object: an `ObservedInfo` or an `ExpectedInfo` with
#'   `method = NULL`.
#'
#' @examples
#' default_information(gaussian1_distrib())
#' default_information(laplace_distrib())
#'
#' @seealso [information_objects], [observed_info()], [expected_info()],
#'   [params_order()]
#' @export
default_information <- function(distrib) {
  ord <- params_order(distrib)
  if (isTRUE(all(ord >= 2))) observed_info() else expected_info()
}


#' Whether Two Information Objects Select the Same Matrix
#'
#' @description
#' Returns `TRUE` when both objects are [observed_info()], or both are
#' [expected_info()] with the same `method` (two `NULL` values included).
#'
#' @param a,b [Information] objects.
#'
#' @return A single logical.
#'
#' @seealso [default_information()]
#' @keywords internal
same_information <- function(a, b) {
  if (S7::S7_inherits(a, ObservedInfo) && S7::S7_inherits(b, ObservedInfo)) {
    return(TRUE)
  }
  if (S7::S7_inherits(a, ExpectedInfo) && S7::S7_inherits(b, ExpectedInfo)) {
    return(identical(a@method, b@method))
  }
  FALSE
}


#' Resolve the Information Argument
#'
#' @description
#' Returns `information` when it is an object and [default_information()] of the
#' family when it is `NULL`.
#'
#' @param distrib An object inheriting from `distrib`.
#' @param information `NULL`, an [observed_info()] or an [expected_info()]
#'   object.
#'
#' @return An [Information] object.
#'
#' @seealso [information_objects], [default_information()]
#' @keywords internal
resolve_information <- function(distrib, information) {
  check_information(information)
  if (is.null(information)) default_information(distrib) else information
}


#' @title Print an Information Object
#' @name print.Information
#'
#' @description
#' Writes one line naming the matrix and, for [expected_info()], the method
#' that obtains the expectation.
#'
#' @param x An object inheriting from [Information].
#' @param ... Unused, accepted for compatibility with [base::print()].
#'
#' @return `x`, invisibly. The function is called for the output that it writes.
#'
#' @examples
#' observed_info()
#' expected_info(method = expect_mc(nsim = 500))
#'
#' @seealso [information_objects]
#' @keywords internal
S7::method(print, Information) <- function(x, ...) {
  if (S7::S7_inherits(x, ExpectedInfo)) {
    lab <- if (is.null(x@method)) {
      "closed form where the family has one, otherwise the outer product of the scores"
    } else {
      sub("^Expectation by ", "expectation by ", method_label(x@method))
    }
    cat(sprintf("Expected information (%s)\n", lab))
  } else if (S7::S7_inherits(x, ObservedInfo)) {
    cat("Observed information\n")
  } else {
    cat("Information of class", class(x)[[1L]], "\n")
  }
  invisible(x)
}
