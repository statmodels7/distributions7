#' @include distrib.R
NULL

#' @title S7 Classes for the Approximation of an Expectation
#' @name ExpectApprox
#' @aliases ExpectOPG ExpectBartlett ExpectIntegrate ExpectMC
#'
#' @description
#' `ExpectApprox` is the abstract parent of the objects that select how an
#' expected derivative is computed when a family has no closed form for it.
#' Its four subclasses are `ExpectOPG`, `ExpectBartlett`, `ExpectIntegrate`
#' and `ExpectMC`. They are built by [expect_opg()], [expect_bartlett()],
#' [expect_integrate()] and [expect_mc()], which are the functions to call;
#' the classes are exported so that `S7::S7_inherits()` can test an object.
#' [expected_derivative_methods] compares the four.
#'
#' `ExpectOPG` and `ExpectBartlett` have no properties. `ExpectMC` has one,
#' `nsim`, the number of draws. `ExpectIntegrate` has six, the tolerances and
#' limits of the quadrature and of the series summation (see
#' [expect_integrate()]). The raw constructors `ExpectMC()` and
#' `ExpectIntegrate()` do not validate their properties, whereas [expect_mc()]
#' and [expect_integrate()] do.
#'
#' @param nsim The number of Monte Carlo draws, a single positive number, for
#'   `ExpectMC`.
#' @param quad_atol,quad_rtol,max_panels,series_atol,series_rtol,max_terms
#'   The properties of `ExpectIntegrate`, documented in [expect_integrate()].
#'
#' @return An S7 object of the class named, which inherits from `ExpectApprox`.
#'
#' @examples
#' a <- expect_mc(nsim = 500)
#' S7::S7_inherits(a, ExpectApprox)
#' S7::S7_inherits(a, ExpectMC)
#' a@nsim
#'
#' @seealso [expected_derivative_methods] for the comparison of the four
#'   approximations; [expect_opg()], [expect_bartlett()], [expect_integrate()]
#'   and [expect_mc()] for the constructors.
#' @export
ExpectApprox <- S7::new_class("ExpectApprox", abstract = TRUE)

#' @rdname ExpectApprox
#' @export
ExpectOPG <- S7::new_class("ExpectOPG", parent = ExpectApprox)

#' @rdname ExpectApprox
#' @export
ExpectBartlett <- S7::new_class("ExpectBartlett", parent = ExpectApprox)

#' @rdname ExpectApprox
#' @export
ExpectIntegrate <- S7::new_class("ExpectIntegrate", parent = ExpectApprox,
  properties = list(
    quad_atol = S7::new_property(S7::class_numeric, default = 1e-10),
    quad_rtol = S7::new_property(S7::class_numeric, default = 1e-8),
    max_panels = S7::new_property(S7::class_numeric, default = 4096L),
    series_atol = S7::new_property(S7::class_numeric, default = 1e-12),
    series_rtol = S7::new_property(S7::class_numeric, default = 1e-10),
    max_terms = S7::new_property(S7::class_numeric, default = 100000L)
  )
)

#' @rdname ExpectApprox
#' @export
ExpectMC <- S7::new_class("ExpectMC", parent = ExpectApprox,
  properties = list(nsim = S7::class_numeric)
)

#' @title The Outer Product of the Scores as the Expected Information
#' @name expect_opg
#'
#' @description
#' Selects the outer product of the observed scores (OPG, the estimator of
#' Berndt, Hall, Hall and Hausman) as the approximation of an expected
#' derivative. At order 2 the expected Hessian is replaced, observation by
#' observation, by the quantity
#' \deqn{-\,\ell_i(y)\,\ell_j(y),}
#' where \eqn{\ell_i} is the derivative of the log-density with respect to the
#' \eqn{i}th parameter. No expectation is taken, so the cost is one evaluation
#' of the score. The matrix is positive semidefinite and has the same fixed
#' point as the expected information, which makes it adequate for the step of
#' Fisher scoring; as an estimate of the information at a single observation it
#' is noisy, and a standard error should be read from the observed information
#' (see [fit_distrib()]).
#'
#' The outer product of the scores has no counterpart above order 2. For the
#' third and fourth derivatives, `expect_opg()` is routed to the Bartlett
#' computation of [expect_bartlett()].
#'
#' This is the default of [distrib_expected_hessian()] and of
#' [fisher_scoring()]. [expected_derivative_methods] compares it with the other
#' approximations.
#'
#' @return An S7 object of class `ExpectOPG`, which inherits from
#'   [ExpectApprox], with no properties.
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
#' distrib_expected_hessian(d, c(-1, 0.5, 2), th, approx = expect_opg())$mu_mu
#'
#' # The same object selects the step of Fisher scoring.
#' fisher_scoring(approx = expect_opg())
#'
#' @seealso [expect_bartlett()], [expect_integrate()] and [expect_mc()] for the
#'   other approximations; [expected_derivative_methods] for the comparison;
#'   [distrib_expected_hessian()] and [fisher_scoring()] for the functions that
#'   take the object.
#' @export
expect_opg <- function() ExpectOPG()

#' @title The Bartlett Identity as the Expected Derivative
#' @name expect_bartlett
#'
#' @description
#' Selects the Bartlett identity as the approximation of an expected
#' derivative. The expected derivative of order \eqn{k} is obtained from
#' expectations of products of lower-order derivatives, through the identity
#' \deqn{\sum_{\pi} \mathbb{E}\left[\prod_{B \in \pi} \ell_B\right] = 0,}
#' where \eqn{\pi} runs over the partitions of the index set and \eqn{\ell_B}
#' is the derivative of the log-density with respect to the parameters indexed
#' by the block \eqn{B}. The single-block partition is the target, so that it
#' equals minus the sum of the others. At order 2 the identity reads
#' \eqn{\mathbb{E}[\ell_{ij}] = -\mathbb{E}[\ell_i\ell_j]}, and the expectation
#' is evaluated as written, by a sum over the support for a discrete family and
#' by quadrature for a continuous one.
#'
#' The computation is deterministic and needs only derivatives of order lower
#' than the target, which is the reason to choose it for a family whose
#' top-order derivative is not available in closed form. It also remains valid
#' when the log-likelihood has a kink in a parameter, where the pointwise
#' average of the observed second derivative omits a point mass. It costs one
#' integral or series per partition, and it is dearer than [expect_opg()] at
#' order 2 and usually dearer than [expect_integrate()] at higher orders.
#'
#' @return An S7 object of class `ExpectBartlett`, which inherits from
#'   [ExpectApprox], with no properties.
#'
#' @examples
#' d <- pseudohuber_distrib()
#' th <- list(mu = 0, sigma = 1, nu = 1)
#' # A family with a closed form ignores the object.
#' identical(distrib_expected_hessian(d, 0, th, approx = expect_bartlett()),
#'           distrib_expected_hessian(d, 0, th))
#'
#' # The third expected derivative of a family without a closed form.
#' Dens <- S7::new_class("DensBartlett", parent = continuous_distrib)
#' S7::method(distrib_pdf, Dens) <- function(distrib, y, theta, log = FALSE, ...) {
#'   stats::dnorm(y, theta[[1]], theta[[2]], log = log)
#' }
#' dd <- Dens(distrib_name = "dens", dimension = "univariate",
#'            bounds = c(-Inf, Inf), params = c("mu", "sigma"),
#'            params_interpretation = c(mu = "location", sigma = "scale"),
#'            n_params = 2L,
#'            params_bounds = list(mu = c(-Inf, Inf), sigma = c(0, Inf)),
#'            link_params = list(mu = linkfunctions7::identity_link(),
#'                               sigma = linkfunctions7::log_link()),
#'            params_smooth = c(mu = TRUE, sigma = TRUE))
#' distrib_deriv3(dd, 0, list(mu = 0, sigma = 2), expected = TRUE,
#'                approx = expect_bartlett())$mu_mu_sigma
#'
#' @seealso [expect_opg()], [expect_integrate()] and [expect_mc()] for the
#'   other approximations; [expected_derivative_methods] for the comparison;
#'   [distrib_deriv3()], [distrib_deriv4()] and [distrib_expected_hessian()]
#'   for the functions that take the object.
#' @export
expect_bartlett <- function() ExpectBartlett()

#' @title Integration Against the Density as the Expected Derivative
#' @name expect_integrate
#'
#' @description
#' Selects numerical integration over the response as the approximation of an
#' expected derivative. The observed derivative of the requested order is
#' integrated against the density,
#' \deqn{\mathbb{E}[\partial^k \ell] = \int \partial^k \ell(y;\theta)\,
#' f(y;\theta)\,dy,}
#' by batched adaptive quadrature for a continuous family and by series
#' summation for a discrete one. The computation is deterministic and is
#' normally the most accurate of the approximations when the observed
#' derivative is available in closed form. The accuracy of the quadrature and
#' of the series is set by the arguments of the constructor, whose defaults are
#' those of [numericals7::quad_vec()] and [numericals7::series_vec()].
#'
#' At order 2, for a family that has a non-smooth parameter (see
#' [param_smoothness()]), the observed second derivative lacks the point mass
#' at the kink, and the expectation is taken as in [expect_bartlett()].
#'
#' This is the default of [distrib_deriv3()] and [distrib_deriv4()].
#'
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
#'   [ExpectApprox], with the properties `quad_atol`, `quad_rtol`,
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
#' # The expected information of sigma is 2 / sigma^2 = 0.5 in absolute value.
#' distrib_expected_hessian(d, 0, th, approx = expect_integrate())$sigma_sigma
#'
#' # Tighter tolerances for the quadrature are set in the object.
#' a <- expect_integrate(quad_atol = 1e-12, quad_rtol = 1e-10)
#' a
#' distrib_expected_hessian(d, 0, th, approx = a)$sigma_sigma
#'
#' @seealso [expect_opg()], [expect_bartlett()] and [expect_mc()] for the other
#'   approximations; [expected_derivative_methods] for the comparison;
#'   [distrib_deriv3()], [distrib_deriv4()] and [distrib_expected_hessian()]
#'   for the functions that take the object.
#' @export
expect_integrate <- function(quad_atol = 1e-10, quad_rtol = 1e-8,
                             max_panels = 4096L,
                             series_atol = 1e-12, series_rtol = 1e-10,
                             max_terms = 100000L) {
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
  do.call(ExpectIntegrate, args)
}

#' @title Monte Carlo as the Expected Derivative
#' @name expect_mc
#'
#' @description
#' Selects Monte Carlo as the approximation of an expected derivative. The
#' function draws `nsim` observations from the distribution at each parameter
#' configuration, with [distrib_rng()], and averages the observed derivative of
#' the requested order over them:
#' \deqn{\widehat{\mathbb{E}}[\partial^k \ell] = \frac{1}{n_{sim}}
#' \sum_{m=1}^{n_{sim}} \partial^k \ell(y_m;\theta).}
#' It estimates the same quantity as [expect_integrate()]. The error decreases
#' as \eqn{1/\sqrt{n_{sim}}}, and the result is reproducible only if the seed is
#' fixed. The cost is that of [distrib_rng()] times `nsim`, per distinct
#' parameter configuration. The choice is adequate where quadrature performs
#' poorly, as with heavy tails or an awkward support. For a family with a
#' non-smooth parameter the average of the observed second derivative omits the
#' point mass at the kink, so that on the Laplace location it returns 0 where
#' the information is \eqn{1/\sigma^2}.
#'
#' @param nsim The number of draws per parameter configuration, a single
#'   positive finite number. The default is 10000. A vector, a missing or
#'   non-finite value, or a number below 1 signals an error.
#'
#' @return An S7 object of class `ExpectMC`, which inherits from
#'   [ExpectApprox], with the property `nsim`.
#'
#' @examples
#' a <- expect_mc(nsim = 2000)
#' a
#'
#' # A sample is drawn at every call, so the seed fixes the result.
#' d <- gaussian1_distrib()
#' th <- list(mu = 0, sigma = 2)
#' set.seed(1)
#' distrib_deriv3(d, 0, th, expected = TRUE, approx = a)$mu_mu_sigma
#'
#' # The step of Fisher scoring can use the same object.
#' fisher_scoring(approx = expect_mc(nsim = 1000))
#'
#' @seealso [expect_opg()], [expect_bartlett()] and [expect_integrate()] for
#'   the other approximations; [expected_derivative_methods] for the
#'   comparison; [distrib_deriv3()], [distrib_deriv4()],
#'   [distrib_expected_hessian()] and [fisher_scoring()] for the functions that
#'   take the object.
#' @export
expect_mc <- function(nsim = 10000) {
  if (!is.numeric(nsim) || length(nsim) != 1L || !is.finite(nsim) || nsim < 1) {
    stop("'nsim' must be a single positive number.", call. = FALSE)
  }
  ExpectMC(nsim = nsim)
}

#' Monte Carlo Sample Size of an Approximation
#'
#' @description
#' Returns the number of draws carried by an [expect_mc()] object, and 10000,
#' the default of [expect_mc()], for any other approximation. The multivariate
#' third and fourth derivatives take their expectation by sampling whatever
#' the class of `approx` and read the sample size here.
#'
#' @param approx An object inheriting from [ExpectApprox].
#'
#' @return A single positive number.
#'
#' @seealso [expect_mc()], [expected_derivative_methods]
#' @keywords internal
approx_nsim <- function(approx) {
  if (S7::S7_inherits(approx, ExpectMC)) approx@nsim else 10000
}

#' Check That an Argument Is an Expectation Approximation
#'
#' @description
#' Returns `approx` invisibly when it inherits from [ExpectApprox] and
#' signals an error otherwise.
#'
#' @param approx The value to check.
#' @param arg The argument name used in the message.
#'
#' @return `approx`, invisibly.
#'
#' @seealso [expected_derivative_methods]
#' @keywords internal
check_approx <- function(approx, arg = "approx") {
  if (S7::S7_inherits(approx, ExpectApprox)) return(invisible(approx))
  stop(sprintf(paste0(
    "'%s' must be an ExpectApprox object such as expect_opg(), ",
    "expect_bartlett(), expect_integrate() or expect_mc()."), arg),
    call. = FALSE)
}

#' @title Print an Expectation Approximation
#' @name print.ExpectApprox
#'
#' @description
#' Writes one line naming the approximation and, for [expect_mc()], the number
#' of draws.
#'
#' @param x An object inheriting from [ExpectApprox].
#' @param ... Unused, accepted for compatibility with [base::print()].
#'
#' @return `x`, invisibly. The function is called for the output that it writes.
#'
#' @examples
#' expect_opg()
#' expect_mc(nsim = 500)
#'
#' @seealso [expected_derivative_methods] for what each approximation computes.
#' @keywords internal
S7::method(print, ExpectApprox) <- function(x, ...) {
  cat(approx_label(x), "\n", sep = "")
  invisible(x)
}

#' Description of an Expectation Approximation
#'
#' @description
#' Returns the sentence that the print method of [ExpectApprox] writes for an
#' approximation, including the sample size of [expect_mc()].
#'
#' @param x An object inheriting from [ExpectApprox].
#'
#' @return A single character string.
#'
#' @seealso [print.ExpectApprox]
#' @keywords internal
approx_label <- function(x) {
  if (S7::S7_inherits(x, ExpectMC)) {
    sprintf("Expectation by Monte Carlo (nsim = %s)", format(x@nsim, scientific = FALSE))
  } else if (S7::S7_inherits(x, ExpectOPG)) {
    "Expectation by the outer product of the scores"
  } else if (S7::S7_inherits(x, ExpectBartlett)) {
    "Expectation by the Bartlett identity"
  } else if (S7::S7_inherits(x, ExpectIntegrate)) {
    base <- "Expectation by integration against the density"
    ref <- ExpectIntegrate()
    set <- c("quad_atol", "quad_rtol", "max_panels", "series_atol",
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
    paste("Expectation approximation of class", class(x)[[1L]])
  }
}

#' Expected Derivative by the Route of an Approximation
#'
#' @description
#' Internal generic that dispatches on the class of `approx` and computes the
#' expected derivative of a given order. The method for [ExpectOPG] calls
#' [expected_by_opg()] at order 2 and [expected_by_bartlett()] above it; those
#' for [ExpectBartlett], [ExpectIntegrate] and [ExpectMC] call
#' [expected_by_bartlett()], [expected_by_integrate()] and
#' [expected_by_mc()]. A subclass of [ExpectApprox] defined elsewhere
#' registers its own method.
#'
#' @param approx An object inheriting from [ExpectApprox].
#' @param distrib An object inheriting from `distrib`.
#' @param y A numeric vector of observations; only its length is used.
#' @param theta A named list of parameters.
#' @param order The derivative order, 2 to 4.
#'
#' @return A named list of expected derivative component vectors.
#'
#' @seealso [expected_derivative()], [expected_derivative_methods]
#' @keywords internal
expected_by_approx <- S7::new_generic("expected_by_approx", "approx",
  function(approx, distrib, y, theta, order) S7::S7_dispatch())

S7::method(expected_by_approx, ExpectOPG) <- function(approx, distrib, y, theta, order) {
  # the outer product of scores is the SECOND-order identity and has no
  # counterpart above it
  if (order == 2L) expected_by_opg(distrib, y, theta)
  else expected_by_bartlett(distrib, y, theta, order)
}

S7::method(expected_by_approx, ExpectBartlett) <- function(approx, distrib, y, theta, order) {
  expected_by_bartlett(distrib, y, theta, order)
}

S7::method(expected_by_approx, ExpectIntegrate) <- function(approx, distrib, y, theta, order) {
  expected_by_integrate(distrib, y, theta, order, approx = approx)
}

S7::method(expected_by_approx, ExpectMC) <- function(approx, distrib, y, theta, order) {
  expected_by_mc(distrib, y, theta, order, approx@nsim)
}
