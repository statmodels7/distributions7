#' @include distrib.R
NULL

#' @title S7 Classes for the Method of an Expectation
#' @name ExpectMethod
#' @aliases ExpectOPG ExpectIntegrate ExpectMC
#'
#' @description
#' `ExpectMethod` is the abstract parent of the objects that select how an
#' expected derivative is obtained. Its three subclasses are `ExpectOPG`,
#' `ExpectIntegrate` and `ExpectMC`. They are built by [expect_opg()],
#' [expect_integrate()] and [expect_mc()], which are the functions to call; the
#' classes are exported so that `S7::S7_inherits()` can test an object.
#' [expected_derivative_methods] compares the three.
#'
#' `ExpectOPG` has no properties. `ExpectMC` has two, `nsim`, the number of
#' draws, and `identity`. `ExpectIntegrate` has seven, `identity` and the
#' tolerances and limits of the quadrature and of the series summation (see
#' [expect_integrate()]). The raw constructors `ExpectMC()` and
#' `ExpectIntegrate()` do not validate their properties, whereas [expect_mc()]
#' and [expect_integrate()] do.
#'
#' @param nsim The number of Monte Carlo draws, a single positive number, for
#'   `ExpectMC`.
#' @param identity One of `"hessian"` and `"score"`, the identity that the
#'   expectation is taken of, for `ExpectIntegrate` and `ExpectMC`. It is
#'   documented in [expected_derivative_methods].
#' @param quad_atol,quad_rtol,max_panels,series_atol,series_rtol,max_terms
#'   The properties of `ExpectIntegrate` that set the quadrature and the series
#'   summation, documented in [expect_integrate()].
#'
#' @return An S7 object of the class named, which inherits from `ExpectMethod`.
#'
#' @examples
#' a <- expect_mc(nsim = 500)
#' S7::S7_inherits(a, ExpectMethod)
#' S7::S7_inherits(a, ExpectMC)
#' a@nsim
#' a@identity
#'
#' @seealso [expected_derivative_methods] for the comparison of the three
#'   methods; [expect_opg()], [expect_integrate()] and [expect_mc()] for the
#'   constructors.
#' @export
ExpectMethod <- S7::new_class("ExpectMethod", abstract = TRUE)

#' @rdname ExpectMethod
#' @export
ExpectOPG <- S7::new_class("ExpectOPG", parent = ExpectMethod)

#' @rdname ExpectMethod
#' @export
ExpectIntegrate <- S7::new_class("ExpectIntegrate", parent = ExpectMethod,
  properties = list(
    identity = S7::new_property(S7::class_character, default = "hessian"),
    quad_atol = S7::new_property(S7::class_numeric, default = 1e-10),
    quad_rtol = S7::new_property(S7::class_numeric, default = 1e-8),
    max_panels = S7::new_property(S7::class_numeric, default = 4096L),
    series_atol = S7::new_property(S7::class_numeric, default = 1e-12),
    series_rtol = S7::new_property(S7::class_numeric, default = 1e-10),
    max_terms = S7::new_property(S7::class_numeric, default = 100000L)
  )
)

#' @rdname ExpectMethod
#' @export
ExpectMC <- S7::new_class("ExpectMC", parent = ExpectMethod,
  properties = list(
    nsim = S7::class_numeric,
    identity = S7::new_property(S7::class_character, default = "hessian")
  )
)

#' @title The Outer Product of the Scores as the Expected Information
#' @name expect_opg
#'
#' @description
#' Selects the empirical outer product of the observed scores (OPG, the
#' estimator of Berndt, Hall, Hall and Hausman) as the estimate of an expected
#' derivative. No expectation is taken and nothing is simulated: the quantity
#' that the score identity (see [expected_derivative_methods]) averages over
#' the model is read at the observations in hand. At order 2 the expected
#' Hessian is replaced, observation by observation, by
#' \deqn{-\,\ell_i(y)\,\ell_j(y),}
#' where \eqn{\ell_i} is the derivative of the log-density with respect to the
#' \eqn{i}th parameter. At order 3 the expected third derivative is replaced by
#' \deqn{-\left(\ell_{ij}\ell_k + \ell_{ik}\ell_j + \ell_{jk}\ell_i
#' + \ell_i\ell_j\ell_k\right),}
#' and at order 4 by the corresponding sum over the partitions of the four
#' indices into at least two blocks, with the sign changed. Every term is a
#' product of observed derivatives of lower order than the target, so the cost
#' is one evaluation of each of those orders and the target order itself is
#' never needed.
#'
#' The matrix at order 2 is positive semidefinite and has the same fixed point
#' as the expected information, which makes it adequate for the step of Fisher
#' scoring. As an estimate of the information it has the variance of a single
#' product of scores per observation, and at orders 3 and 4 the variance of a
#' product of up to four derivatives; the sum over \eqn{n} observations
#' converges to \eqn{n} times the expectation at the rate \eqn{n^{-1/2}}. A
#' standard error should be read from the observed information (see
#' [fit_distrib()]).
#'
#' This is the default of [distrib_expected_hessian()] for a family without a
#' closed form, and of [fisher_scoring()] and [expected_info()]. At orders 3 and
#' 4 it is never the default and is never replaced by another method: a call
#' with `method = expect_opg()` returns the empirical sums of products
#' described above. [expected_derivative_methods] compares it with the other
#' methods.
#'
#' @return An S7 object of class `ExpectOPG`, which inherits from
#'   [ExpectMethod], with no properties.
#'
#' @examples
#' # A family defined by its density alone has no closed-form information.
#' Dens <- S7::new_class("DensOPG", parent = continuous_distrib)
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
#' distrib_expected_hessian(d, c(-1, 0.5, 2), th, method = expect_opg())$mu_mu
#'
#' # The same object is accepted at the higher orders.
#' distrib_deriv3(d, c(-1, 0.5, 2), th, expected = TRUE,
#'                method = expect_opg())$mu_mu_sigma
#'
#' # And it selects the step of Fisher scoring.
#' fisher_scoring(method = expect_opg())
#'
#' @seealso [expect_integrate()] and [expect_mc()] for the other methods;
#'   [expected_derivative_methods] for the comparison;
#'   [distrib_expected_hessian()], [distrib_deriv3()], [distrib_deriv4()] and
#'   [fisher_scoring()] for the functions that take the object.
#' @export
expect_opg <- function() ExpectOPG()

#' @title Integration Against the Density as the Expected Derivative
#' @name expect_integrate
#'
#' @description
#' Selects numerical integration over the response as the method that obtains
#' an expected derivative. The quantity integrated against the density \eqn{f}
#' is the one that `identity` names:
#' \deqn{\mathbb{E}[\partial^k \ell] = \int \partial^k \ell(y;\theta)\,
#' f(y;\theta)\,dy \quad \text{for } \texttt{identity = "hessian"},}
#' the observed derivative of the requested order, or
#' \deqn{-\sum_{\pi} \int \prod_{B \in \pi} \ell_B(y;\theta)\,f(y;\theta)\,dy
#' \quad \text{for } \texttt{identity = "score"},}
#' where \eqn{\pi} runs over the partitions of the index set into at least two
#' blocks and \eqn{\ell_B} is the derivative of the log-density with respect to
#' the parameters indexed by the block \eqn{B} (the Bartlett partition
#' identity; at order 2 it reads \eqn{-\int \ell_i \ell_j f\,dy}). The
#' integrals are computed by batched adaptive quadrature for a continuous
#' family and by series summation for a discrete one. The computation is
#' deterministic. The accuracy of the quadrature and of the series is set by
#' the arguments of the constructor, whose defaults are those of
#' [numericals7::quad_vec()] and [numericals7::series_vec()].
#'
#' With `identity = "hessian"` the integral needs the derivative of the target
#' order and costs one integral per component; it is normally the most
#' accurate method when that derivative is available in closed form. With
#' `identity = "score"` the target order is never needed, only the lower ones,
#' and the cost is one integral per partition and component. The score identity
#' is the one that remains valid when the log-likelihood has a kink in a
#' parameter: at order 2, for a family with a non-smooth parameter (see
#' [param_smoothness()]), the observed second derivative lacks the point mass
#' at the kink, and `identity = "hessian"` signals an error that names
#' `identity = "score"`.
#'
#' With the defaults this is the method that [distrib_deriv3()] and
#' [distrib_deriv4()] use for an expected derivative when `method` is `NULL`.
#'
#' @param identity One of `"hessian"` (the default) and `"score"`. See the
#'   description and [expected_derivative_methods].
#' @param quad_atol,quad_rtol The absolute and the relative tolerance of the
#'   quadrature of a continuous family, single positive finite numbers. A
#'   panel of the integration is accepted when its error estimate is below
#'   `quad_atol + quad_rtol * abs(integral)`. The defaults are 1e-10 and 1e-8.
#' @param max_panels The largest number of panels that the quadrature of one
#'   parameter configuration may reach, a single positive finite number. A row
#'   that needs more panels is not accepted by the batched quadrature and is
#'   integrated by the scalar integrator of \pkg{stats}. The default is 4096.
#' @param series_atol,series_rtol The absolute and the relative tolerance of
#'   the series summation of a discrete family, single positive finite
#'   numbers. The defaults are 1e-12 and 1e-10.
#' @param max_terms The largest number of terms of a series, a single positive
#'   finite number. A series that has not met its tolerance after that many
#'   terms signals an error. The default is 100000.
#'
#' @return An S7 object of class `ExpectIntegrate`, which inherits from
#'   [ExpectMethod], with the properties `identity`, `quad_atol`, `quad_rtol`,
#'   `max_panels`, `series_atol`, `series_rtol` and `max_terms`.
#'
#' @examples
#' Dens <- S7::new_class("DensIntegrate", parent = continuous_distrib)
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
#' # The expected information of sigma is 2 / sigma^2 = 0.5 in absolute value,
#' # by either identity.
#' distrib_expected_hessian(d, 0, th, method = expect_integrate())$sigma_sigma
#' distrib_expected_hessian(d, 0, th,
#'                          method = expect_integrate(identity = "score"))$sigma_sigma
#'
#' # Tighter tolerances for the quadrature are set in the object.
#' a <- expect_integrate(quad_atol = 1e-12, quad_rtol = 1e-10)
#' a
#' distrib_expected_hessian(d, 0, th, method = a)$sigma_sigma
#'
#' @seealso [expect_opg()] and [expect_mc()] for the other methods;
#'   [expected_derivative_methods] for the comparison;
#'   [distrib_deriv3()], [distrib_deriv4()] and [distrib_expected_hessian()]
#'   for the functions that take the object.
#' @export
expect_integrate <- function(identity = c("hessian", "score"),
                             quad_atol = 1e-10, quad_rtol = 1e-8,
                             max_panels = 4096L,
                             series_atol = 1e-12, series_rtol = 1e-10,
                             max_terms = 100000L) {
  identity <- match.arg(identity)
  args <- list(quad_atol = quad_atol, quad_rtol = quad_rtol,
               max_panels = max_panels, series_atol = series_atol,
               series_rtol = series_rtol, max_terms = max_terms)
  for (nm in names(args)) {
    v <- args[[nm]]
    if (!is.numeric(v) || length(v) != 1L || !is.finite(v) || v <= 0) {
      stop(sprintf("'%s' must be a single positive finite number.", nm),
           call. = FALSE)
    }
  }
  do.call(ExpectIntegrate, c(list(identity = identity), args))
}

#' @title Monte Carlo as the Expected Derivative
#' @name expect_mc
#'
#' @description
#' Selects Monte Carlo as the method that obtains an expected derivative. The
#' function draws `nsim` observations from the distribution at each parameter
#' configuration, with [distrib_rng()], and averages over them the quantity
#' that `identity` names:
#' \deqn{\widehat{\mathbb{E}}[\partial^k \ell] = \frac{1}{n_{sim}}
#' \sum_{m=1}^{n_{sim}} \partial^k \ell(y_m;\theta)}
#' for `identity = "hessian"`, the observed derivative of the requested order,
#' and the same average of
#' \eqn{-\sum_{\pi} \prod_{B \in \pi} \ell_B(y_m;\theta)} for
#' `identity = "score"`, the Bartlett partition identity (see
#' [expect_integrate()] and [expected_derivative_methods]). Both estimate the
#' same quantity as [expect_integrate()] with the same `identity`. The error
#' decreases as \eqn{1/\sqrt{n_{sim}}}, and the result is reproducible only if
#' the seed is fixed. The cost is that of [distrib_rng()] times `nsim`, per
#' distinct parameter configuration. The choice is adequate where quadrature
#' performs poorly, as with heavy tails or an awkward support.
#'
#' For a family with a non-smooth parameter (see [param_smoothness()]) the
#' average of the observed second derivative omits the point mass at the kink,
#' so that `identity = "hessian"` at order 2 signals an error that names
#' `identity = "score"`.
#'
#' @param nsim The number of draws per parameter configuration, a single
#'   positive finite number. The default is 10000. A vector, a missing or
#'   non-finite value, or a number below 1 signals an error.
#' @param identity One of `"hessian"` (the default) and `"score"`. See the
#'   description and [expected_derivative_methods].
#'
#' @return An S7 object of class `ExpectMC`, which inherits from
#'   [ExpectMethod], with the properties `nsim` and `identity`.
#'
#' @examples
#' a <- expect_mc(nsim = 2000)
#' a
#' expect_mc(nsim = 2000, identity = "score")
#'
#' # A sample is drawn at every call, so the seed fixes the result.
#' d <- gaussian1_distrib()
#' th <- list(mu = 0, sigma = 2)
#' set.seed(1)
#' distrib_deriv3(d, 0, th, expected = TRUE, method = a)$mu_mu_sigma
#'
#' # The step of Fisher scoring can use the same object.
#' fisher_scoring(method = expect_mc(nsim = 1000))
#'
#' @seealso [expect_opg()] and [expect_integrate()] for the other methods;
#'   [expected_derivative_methods] for the comparison; [distrib_deriv3()],
#'   [distrib_deriv4()], [distrib_expected_hessian()] and [fisher_scoring()] for
#'   the functions that take the object.
#' @export
expect_mc <- function(nsim = 10000, identity = c("hessian", "score")) {
  identity <- match.arg(identity)
  if (!is.numeric(nsim) || length(nsim) != 1L || !is.finite(nsim) || nsim < 1) {
    stop("'nsim' must be a single positive number.", call. = FALSE)
  }
  ExpectMC(nsim = nsim, identity = identity)
}

#' Monte Carlo Sample Size of an Expectation Method
#'
#' @description
#' Returns the number of draws carried by an [expect_mc()] object, and 10000,
#' the default of [expect_mc()], for any other method. The multivariate third
#' and fourth derivatives take their expectation by sampling whatever the class
#' of `method` and read the sample size here.
#'
#' @param method An object inheriting from [ExpectMethod], or `NULL`.
#'
#' @return A single positive number.
#'
#' @seealso [expect_mc()], [expected_derivative_methods]
#' @keywords internal
method_nsim <- function(method) {
  if (S7::S7_inherits(method, ExpectMC)) method@nsim else 10000
}

#' Check That an Argument Is an Expectation Method
#'
#' @description
#' Returns `method` invisibly when it is `NULL` (if `null_ok`) or inherits from
#' [ExpectMethod], and signals an error otherwise.
#'
#' @param method The value to check.
#' @param arg The argument name used in the message.
#' @param null_ok Whether `NULL` is accepted.
#'
#' @return `method`, invisibly.
#'
#' @seealso [expected_derivative_methods]
#' @keywords internal
check_method <- function(method, arg = "method", null_ok = TRUE) {
  if (S7::S7_inherits(method, ExpectMethod)) return(invisible(method))
  if (null_ok && is.null(method)) return(invisible(method))
  stop(sprintf(paste0(
    "'%s' must be %san ExpectMethod object such as expect_opg(), ",
    "expect_integrate() or expect_mc()."), arg, if (null_ok) "NULL or " else ""),
    call. = FALSE)
}

#' Whether an Expectation Method Replaces the Family's Own Route
#'
#' @description
#' Returns `TRUE` when `method` is an object. The generics
#' [distrib_expected_hessian()], [distrib_deriv3()] and [distrib_deriv4()] then
#' compute the expectation with [route_expected()], bypassing the method
#' registered for the family, so that an explicit method is used also by a
#' family that writes the expectation in closed form.
#'
#' @param distrib An object inheriting from `distrib`.
#' @param method An object inheriting from [ExpectMethod], or `NULL`.
#'
#' @return A single logical.
#'
#' @seealso [expected_derivative_methods]
#' @keywords internal
explicit_method <- function(distrib, method) {
  !is.null(method)
}

#' Route an Explicit Expectation Method
#'
#' @description
#' Computes the expected derivative of an order by the method that the object
#' selects. A univariate family goes to [expected_derivative()]. A
#' multivariate family has no quadrature over the response, so it goes to
#' [mv_expected_sampled()], which takes the expectation by sampling.
#'
#' @param distrib An object inheriting from `distrib`.
#' @param y A numeric vector or matrix of observations.
#' @param theta A named list of parameters.
#' @param order The derivative order, 2 to 4.
#' @param method An object inheriting from [ExpectMethod].
#'
#' @return A named list of expected derivative component vectors.
#'
#' @seealso [explicit_method()], [expected_derivative_methods]
#' @keywords internal
route_expected <- function(distrib, y, theta, order, method) {
  if (S7::S7_inherits(distrib, multivariate_distrib)) {
    mv_expected_sampled(distrib, y, theta, order, method)
  } else {
    expected_derivative(distrib, y, theta, order = order, method = method)
  }
}

#' Resolve the Default Method of an Expected Derivative
#'
#' @description
#' Returns `method` when it is an object and, when it is `NULL`, the default
#' of the generics [distrib_expected_hessian()], [distrib_deriv3()] and
#' [distrib_deriv4()] for the order: [expect_opg()] at order 2 and
#' [expect_integrate()] with `identity = "hessian"` above it.
#'
#' @param method An object inheriting from [ExpectMethod], or `NULL`.
#' @param order The derivative order, 2 to 4.
#'
#' @return An object inheriting from [ExpectMethod].
#'
#' @seealso [expected_derivative_methods]
#' @keywords internal
resolve_method <- function(method, order) {
  if (!is.null(method)) return(method)
  if (order == 2L) expect_opg() else expect_integrate()
}

#' @title Print an Expectation Method
#' @name print.ExpectMethod
#'
#' @description
#' Writes one line naming the method and, for [expect_mc()], the number of
#' draws, and for [expect_integrate()] and [expect_mc()] the identity when it is
#' not the default.
#'
#' @param x An object inheriting from [ExpectMethod].
#' @param ... Unused, accepted for compatibility with [base::print()].
#'
#' @return `x`, invisibly. The function is called for the output that it writes.
#'
#' @examples
#' expect_opg()
#' expect_mc(nsim = 500)
#' expect_integrate(identity = "score")
#'
#' @seealso [expected_derivative_methods] for what each method computes.
#' @keywords internal
S7::method(print, ExpectMethod) <- function(x, ...) {
  cat(method_label(x), "\n", sep = "")
  invisible(x)
}

#' Description of an Expectation Method
#'
#' @description
#' Returns the sentence that the print method of [ExpectMethod] writes for a
#' method, including the sample size of [expect_mc()] and the settings of
#' [expect_integrate()] that differ from the defaults.
#'
#' @param x An object inheriting from [ExpectMethod].
#'
#' @return A single character string.
#'
#' @seealso [print.ExpectMethod]
#' @keywords internal
method_label <- function(x) {
  if (S7::S7_inherits(x, ExpectMC)) {
    sc <- if (identical(x@identity, "hessian")) "" else
      sprintf(", identity = %s", x@identity)
    sprintf("Expectation by Monte Carlo (nsim = %s%s)",
            format(x@nsim, scientific = FALSE), sc)
  } else if (S7::S7_inherits(x, ExpectOPG)) {
    "Expectation by the outer product of the scores"
  } else if (S7::S7_inherits(x, ExpectIntegrate)) {
    base <- "Expectation by integration against the density"
    ref <- ExpectIntegrate()
    set <- c("identity", "quad_atol", "quad_rtol", "max_panels", "series_atol",
             "series_rtol", "max_terms")
    set <- set[!vapply(set, function(nm) {
      isTRUE(S7::prop(x, nm) == S7::prop(ref, nm))
    }, logical(1))]
    if (length(set)) {
      val <- vapply(set, function(nm) {
        format(S7::prop(x, nm), scientific = FALSE, digits = 15)
      }, character(1))
      sprintf("%s (%s)", base, paste(set, val, sep = " = ", collapse = ", "))
    } else {
      base
    }
  } else {
    paste("Expectation method of class", class(x)[[1L]])
  }
}

#' Expected Derivative by the Route of a Method
#'
#' @description
#' Internal generic that dispatches on the class of `method` and computes the
#' expected derivative of a given order. The method for [ExpectOPG] calls
#' [expected_by_opg()]. The method for [ExpectIntegrate] calls
#' [expected_by_integrate()] for `identity = "hessian"` and
#' [expected_by_bartlett()] for `identity = "score"`. The method for
#' [ExpectMC] calls [expected_by_mc()] with the identity of the object. A
#' subclass of [ExpectMethod] defined elsewhere registers its own method.
#'
#' @param method An object inheriting from [ExpectMethod].
#' @param distrib An object inheriting from `distrib`.
#' @param y A numeric vector of observations; only its length is used by the
#'   methods that take an expectation.
#' @param theta A named list of parameters.
#' @param order The derivative order, 2 to 4.
#'
#' @return A named list of expected derivative component vectors.
#'
#' @seealso [expected_derivative()], [expected_derivative_methods]
#' @keywords internal
expected_by_method <- S7::new_generic("expected_by_method", "method",
  function(method, distrib, y, theta, order) S7::S7_dispatch())

S7::method(expected_by_method, ExpectOPG) <- function(method, distrib, y, theta, order) {
  expected_by_opg(distrib, y, theta, order)
}

S7::method(expected_by_method, ExpectIntegrate) <- function(method, distrib, y, theta, order) {
  if (identical(method@identity, "score")) {
    old <- set_integration_settings(method)
    on.exit(set_integration_settings(old))
    expected_by_bartlett(distrib, y, theta, order)
  } else {
    expected_by_integrate(distrib, y, theta, order, method = method)
  }
}

S7::method(expected_by_method, ExpectMC) <- function(method, distrib, y, theta, order) {
  expected_by_mc(distrib, y, theta, order, method@nsim, identity = method@identity)
}
