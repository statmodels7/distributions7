#' @include multivariate.R mvgaussian_distrib.R mv_plot.R
NULL

#' @title Multivariate Student t Distribution Class
#' @name MvStudentTDistrib
#'
#' @description
#' The S7 class of the elliptical Student t family on \eqn{\mathbb{R}^p}, with
#' density
#' \deqn{f(y) \propto \lvert\Sigma\rvert^{-1/2}
#'   \left(1 + \frac{(y-\mu)^\top \Sigma^{-1} (y-\mu)}{\nu}\right)^{-(\nu+p)/2}.}
#' The location \eqn{\mu} contributes \eqn{p} scalar parameters, the scale
#' matrix \eqn{\Sigma} is carried by a \pkg{parameters7} parametrization whose
#' free values become scalar parameters, and the degrees of freedom \eqn{\nu}
#' is added last. It inherits from `multivariate_distrib`, so the response is an
#' \eqn{n \times p} matrix and the distribution function and the quantile
#' function signal an error.
#'
#' Build one with [mvstudent_t1_distrib()], which fills the properties in and
#' checks the rank of the parametrization. This page documents the raw S7
#' constructor, which runs only the validators of `multivariate_distrib` and
#' `distrib`; nothing checks `param`.
#'
#' @param param A \pkg{parameters7} parametrization of the scale matrix or of
#'   its inverse, according to `inverted`, inheriting from
#'   `parameters7::parameter`. Its `n_free` free values are flattened into
#'   scalar parameters of the distribution.
#' @param inverted Logical of length 1. `TRUE` when `param` carries the inverse
#'   scale matrix \eqn{\Sigma^{-1}} and `FALSE` when it carries the scale matrix
#'   \eqn{\Sigma}. It decides which derivative arrays are the parametrization's
#'   own and which are obtained from them (see [mv_matrix_pieces()]), the sign
#'   of the log-determinant, the prefix of the free names and the quantities
#'   that [mv_summary()] reports. The same \eqn{\Sigma} gives the same law
#'   either way, and neither matrix is a moment of the response.
#' @param n_dim The dimension \eqn{p} of one observation. A single positive
#'   integer.
#' @inheritParams distrib
#'
#' @return An S7 object of class `MvStudentTDistrib`, inheriting from
#'   `multivariate_distrib` and from `distrib`. Beyond the parent's
#'   `distrib_name`, `dimension`, `bounds`, `params`,
#'   `params_interpretation`, `n_params`, `params_bounds`, `link_params`,
#'   `params_smooth` and `n_dim`, it carries `param` and `inverted` as
#'   supplied.
#'
#' @seealso [mvstudent_t1_distrib()] to build one, [mvgaussian1_distrib()] for
#'   the limit \eqn{\nu \to \infty}, [mv_sigma.MvStudentTDistrib()] for the
#'   scale matrix and [variance.MvStudentTDistrib()] for the covariance, which
#'   are different matrices here.
#'
#' @section Methods:
#' Registered on this class, each linked to the page that documents it:
#' - the distribution: [`distrib_pdf()`][distrib_pdf.MvStudentTDistrib] and
#'   [`distrib_rng()`][distrib_rng.MvStudentTDistrib];
#' - the derivatives in the parameters:
#'   [`distrib_gradient()`][distrib_gradient.MvStudentTDistrib],
#'   [`distrib_hessian()`][distrib_hessian.MvStudentTDistrib],
#'   [`distrib_expected_hessian()`][distrib_expected_hessian.MvStudentTDistrib],
#'   [`distrib_deriv3()`][distrib_deriv3.MvStudentTDistrib] and
#'   [`distrib_deriv4()`][distrib_deriv4.MvStudentTDistrib];
#' - the derivatives in the response:
#'   [`distrib_grad_y()`][distrib_grad_y.MvStudentTDistrib],
#'   [`distrib_hess_y()`][distrib_hess_y.MvStudentTDistrib],
#'   [`distrib_deriv3_y()`][distrib_deriv3_y.MvStudentTDistrib] and
#'   [`distrib_deriv4_y()`][distrib_deriv4_y.MvStudentTDistrib];
#' - the mixed derivatives in the response and the parameters:
#'   [`distrib_cross_y()`][distrib_cross_y.MvStudentTDistrib],
#'   [`distrib_cross2_y()`][distrib_cross2_y.MvStudentTDistrib],
#'   [`distrib_cross3_y()`][distrib_cross3_y.MvStudentTDistrib],
#'   [`distrib_grad_y_hess()`][distrib_grad_y_hess.MvStudentTDistrib] and
#'   [`distrib_hess_y_hess()`][distrib_hess_y_hess.MvStudentTDistrib];
#' - the moments: [`mean()`][mean.MvStudentTDistrib] and
#'   [`variance()`][variance.MvStudentTDistrib];
#' - the multivariate summaries:
#'   [`mv_location()`][mv_location.MvStudentTDistrib],
#'   [`mv_sigma()`][mv_sigma.MvStudentTDistrib],
#'   [`mv_marginal()`][mv_marginal.MvStudentTDistrib] and
#'   [`mv_derived()`][mv_derived.MvStudentTDistrib];
#' - the validation and fitting helpers:
#'   [`mv_reference_draw()`][mv_reference_draw.MvStudentTDistrib],
#'   [`generate_random_theta()`][generate_random_theta.MvStudentTDistrib] and
#'   [`distrib_start()`][distrib_start.MvStudentTDistrib].
#'
#' Every other generic is inherited from [multivariate_distrib()], whose methods
#' signal an error for the distribution function, the quantile function, the
#' skewness and the kurtosis.
#'
#' @examples
#' d <- mvstudent_t1_distrib(2)
#' S7::S7_inherits(d, multivariate_distrib)
#'
#' # Six parameters: two locations, three free values of the scale matrix, and
#' # the degrees of freedom last.
#' d@params
#' d@params_interpretation
#'
#' # nu is the one parameter with a bound, and the one with a link that is not
#' # the identity, so this family's link scale is not its parameter scale.
#' d@params_bounds$nu
#' vapply(d@link_params, function(l) l@link_name, character(1))
#'
#' @export
MvStudentTDistrib <- S7::new_class("MvStudentTDistrib",
  parent = multivariate_distrib,
  properties = list(
    param = parameters7::parameter,
    inverted = S7::class_logical
  )
)


#' @title S7 Classes for the Two Multivariate Student t Parametrizations
#'
#' @description
#' `MvStudentT1Distrib` is the family whose matrix parametrization is the
#' **scale matrix** \eqn{\Sigma} and `MvStudentT2Distrib` the family whose
#' parametrization is its **inverse** \eqn{\Sigma^{-1}}. Both inherit from
#' [MvStudentTDistrib], which carries every method; the two differ in what
#' their free values describe, and the computations are shared.
#'
#' @details
#' Neither matrix is a covariance or a precision of the response; they differ
#' from those by a factor in \eqn{\nu} that [mvstudent_t1_distrib()] writes
#' out. The numbering is the package's convention, the one that separates
#' [mvgaussian1_distrib()] from [mvgaussian2_distrib()] and
#' [gaussian1_distrib()] from [gaussian3_distrib()]: a parametrization decides
#' what a linear predictor acts on, so it has a name of its own.
#'
#' @inheritParams MvStudentTDistrib
#'
#' @return An object of class `MvStudentT1Distrib` or `MvStudentT2Distrib`, a
#'   subclass of [MvStudentTDistrib] adding no properties of its own.
#'
#' @seealso [mvstudent_t1_distrib()] and [mvstudent_t2_distrib()], the
#'   constructors, and [MvStudentTDistrib] for the properties and the methods.
#'
#' @examples
#' c(scale = S7::S7_inherits(mvstudent_t1_distrib(2), MvStudentT1Distrib),
#'   inverse = S7::S7_inherits(mvstudent_t2_distrib(2), MvStudentT2Distrib))
#'
#' @export
MvStudentT1Distrib <- S7::new_class("MvStudentT1Distrib",
  parent = MvStudentTDistrib
)

#' @rdname MvStudentT1Distrib
#' @export
MvStudentT2Distrib <- S7::new_class("MvStudentT2Distrib",
  parent = MvStudentTDistrib
)


#' @title Construct a Multivariate Student t Distribution
#'
#' @description
#' Two families on \eqn{\mathbb{R}^p}, differing in which matrix their
#' parametrization describes: `mvstudent_t1_distrib()` parametrizes the
#' **scale matrix** \eqn{\Sigma} and `mvstudent_t2_distrib()` its **inverse**
#' \eqn{\Sigma^{-1}}. Each carries a location vector of \eqn{p} free
#' parameters, the matrix through a \pkg{parameters7} parametrization, and a
#' positive degrees-of-freedom parameter. The free values become ordinary
#' scalar parameters, so the generics of the `distrib` class accept the object,
#' with `theta` a flat named list of numbers; the distribution function and the
#' quantile function signal an error. Both default to an unstructured matrix in
#' the log-Cholesky parametrization, which has `p * (p + 1) / 2` free values,
#' and to a log link on `nu`.
#'
#' @details
#' # The two matrices
#'
#' Neither matrix is a covariance, and neither is a precision. \eqn{\Sigma} is
#' the **scale matrix**, the matrix in the quadratic form of the density, and
#' the moments are obtained from it by scaling in \eqn{\nu}:
#'
#' \deqn{\operatorname{Var}(Y) = \frac{\nu}{\nu-2}\,\Sigma \quad (\nu > 2),
#'   \qquad
#'   \operatorname{Var}(Y)^{-1} = \frac{\nu-2}{\nu}\,\Sigma^{-1}.}
#'
#' So the matrix of `mvstudent_t2_distrib()` equals the precision of the
#' response only up to the factor \eqn{(\nu-2)/\nu}, which is 2/3 at
#' \eqn{\nu = 6}; a coordinate named `omega_log_L1` is a free value of
#' \eqn{\Sigma^{-1}} and not of the precision. The prefixes follow the
#' convention of [mvgaussian1_distrib()] and name which of the two matrices a
#' free value builds.
#'
#' For \eqn{\nu \le 2} the covariance does not exist while the density is well
#' defined, and for \eqn{\nu \le 1} the mean does not exist either. For that
#' reason [mv_sigma()] returns the scale matrix, which the parametrization
#' carries and which exists at every admissible \eqn{\nu}, and [variance()]
#' returns the covariance, a moment, which is `Inf` in every entry for
#' \eqn{\nu \le 2}; [mean()] returns `NaN` for \eqn{\nu \le 1}. A fit
#' therefore runs at \eqn{\nu = 1.5}.
#'
#' Two readings carry over from the matrix to the response without a factor,
#' because a positive multiple cancels out of a ratio:
#'
#' * the **correlations** of \eqn{\Sigma} are the correlations of \eqn{Y};
#' * the **partial correlations** read off \eqn{\Sigma^{-1}}, as
#'   \eqn{-\Omega_{jk}/\sqrt{\Omega_{jj}\Omega_{kk}}}, are the partial
#'   correlations of \eqn{Y}.
#'
#' The conditional variance does not carry over. For a gaussian
#' \eqn{1/\Omega_{jj}} is the conditional variance of \eqn{Y_j} given the rest;
#' here it is the Schur complement of the scale matrix, so it is a conditional
#' **scale**, and the conditional law of a \eqn{t} is itself a \eqn{t} whose
#' scale carries a further factor depending on the conditioning values.
#' [mv_summary()] names that quantity `cscale_vj` and not `cvar_vj` for this
#' reason, as it names the diagonal quantities `scale_sd_vj`.
#'
#' # Which side to write a structure on
#'
#' Where the matrix parametrization is closed under inversion the two families
#' describe the same laws in different coordinates. Where it is not they are
#' different models: `mvstudent_t2_distrib(p, parameters7::ar1(p))` is the
#' family whose inverse scale matrix has the entries \eqn{\sigma^2\rho^{|i-j|}},
#' with \eqn{\sigma^2} the scale of `ar1()`, and it is a valid law. The AR(1)
#' pattern on the scale matrix with the other side parametrized is
#' `mvstudent_t2_distrib(p, parameters7::ar1_inv(p))`, whose parametrization is
#' that inverse.
#'
#' # Parameters
#'
#' The location contributes `mu1`, ..., `mup`, the parametrization contributes
#' its free values under the `sigma_` prefix for the scale matrix and the
#' `omega_` prefix for its inverse, and `nu` comes last. The location and the
#' matrix parameters are unconstrained and carry identity links; `nu` is
#' positive and carries `link_nu`. Unlike [mvgaussian1_distrib()], this family
#' therefore has a link scale that differs from its parameter scale: under the
#' default log link, `scale = "link"` applies the chain rule of the log link to
#' every derivative in `nu`.
#'
#' # The weight of an observation
#'
#' A gaussian fitted to data with a few outlying rows inflates its covariance
#' to cover them. A \eqn{t} with \eqn{\nu} estimated does not: an observation at
#' squared Mahalanobis distance \eqn{q} enters every derivative through the
#' weight \eqn{c = (\nu+p)/(\nu+q)}, which falls as \eqn{1/q}. The gaussian is
#' the limit \eqn{\nu \to \infty}, where \eqn{c \equiv 1}.
#'
#' # The expected information
#'
#' Closed form, from the family's own scale mixture; see
#' [distrib_expected_hessian.MvStudentTDistrib()] for the four blocks.
#' [fit_distrib()] therefore rejects an `approx` argument on
#' [fisher_scoring()] for this family, as it does for any family that computes
#' its expected information exactly.
#'
#' # Reading a fit
#'
#' For the scale parametrization, [mv_summary()] reports the square roots of
#' the diagonal of the scale matrix, named `scale_sd_v1`, ..., `scale_sd_vp`
#' because they are not standard deviations of the response, and the
#' correlations, which are the response's as well. For the inverse
#' parametrization it reports instead the conditional scales `cscale_v1`, ...,
#' `cscale_vp` and the partial correlations `pcor_vj_vk`, which are the
#' response's exactly.
#'
#' # Degrees of freedom far from the gaussian limit and close to it
#'
#' Every quantity of \eqn{\nu} alone (the constant of the density and its
#' derivatives, and the expected information in \eqn{\nu}) is computed by
#' [mvt_cp_d()] and [mvt_info_nunu()] without the cancellation of the direct
#' difference of gamma and polygamma functions, and the part that carries the
#' data by [mvt_h_d()]; all of them vanish as \eqn{\nu} grows, while the terms
#' of the direct forms do not.
#'
#' # The response
#'
#' An \eqn{n \times p} matrix, one row per observation. A plain vector of
#' length \eqn{p} is read as a single observation. A parameter may not vary
#' from observation to observation.
#'
#' @section The distribution:
#' \deqn{f(y) = \frac{\Gamma\!\left(\frac{\nu+p}{2}\right)}{\Gamma\!\left(\frac{\nu}{2}\right)(\nu\pi)^{p/2}\lvert \Sigma \rvert^{1/2}}\left(1 + \frac{q}{\nu}\right)^{-(\nu+p)/2}, \qquad q = (y-\mu)^\top\Sigma^{-1}(y-\mu)}
#' on \eqn{y \in \mathbb{R}^{p}}, with
#' \deqn{\mathbb{E}[Y] = \mu \; (\nu > 1), \qquad \operatorname{Var}(Y) = \frac{\nu}{\nu-2}\Sigma \; (\nu > 2).}
#'
#' @section Notation:
#' \eqn{\mu} is the location, \eqn{\Sigma} the scale matrix, \eqn{\nu} the
#' degrees of freedom, \eqn{p} the dimension of one observation and \eqn{q} the
#' squared Mahalanobis distance of an observation from the location. \eqn{\eta}
#' is the free vector of the matrix parametrization.
#'
#' @param n_dim The dimension \eqn{p} of one observation. A single positive
#'   whole number, finite and at least 1. Any other value signals an error.
#' @param sigma A \pkg{parameters7} parametrization of the scale matrix, of
#'   dimension `n_dim` and of full rank. Defaults to
#'   `parameters7::log_cholesky(n_dim)`. A rank-deficient parametrization is
#'   rejected, the density not normalizing.
#' @param omega A \pkg{parameters7} parametrization of the inverse scale
#'   matrix, of dimension `n_dim` and of full rank. Defaults to
#'   `parameters7::log_cholesky(n_dim)`. It is the precision of the response
#'   only up to the factor \eqn{(\nu-2)/\nu}; see **Details**.
#' @param link_nu The link carrying the degrees of freedom to the
#'   unconstrained scale, a `linkfunctions7::link` object. Defaults to
#'   `linkfunctions7::log_link()`, which keeps \eqn{\nu} strictly positive at
#'   every point of the free scale.
#'
#' @return An S7 object of class [MvStudentT1Distrib] or [MvStudentT2Distrib],
#'   with `param` the parametrization supplied and `inverted` recording which
#'   matrix it carries. Its `params` are `mu1`, ..., `mup`, the prefixed
#'   free names, then `nu`; `n_params` is `p + param@n_free + 1`; every
#'   `params_bounds` entry is \eqn{(-\infty, \infty)} except that of `nu`,
#'   which is \eqn{(0, \infty)}; and every link is the identity except that of
#'   `nu`.
#'
#' @seealso [mvgaussian1_distrib()] for the limiting family,
#'   [student_t1_distrib()] for the one-dimensional case,
#'   [mv_marginal.MvStudentTDistrib()] for a marginal, which keeps the same
#'   \eqn{\nu}, [parameters7::inverse_of()] for writing a structure on the
#'   other side, and [fit_distrib()] to estimate one.
#'
#' @examples
#' d <- mvstudent_t1_distrib(2)
#' d
#'
#' theta <- list(mu1 = 0.5, mu2 = -0.3, sigma_log_L1 = 0.1,
#'               sigma_log_L2 = -0.2, sigma_L2.1 = 0.4, nu = 6)
#'
#' # The scale matrix is what the parametrization carries; the covariance is a
#' # moment and is larger by nu / (nu - 2).
#' mv_sigma(d, theta)
#' all.equal(variance(d, theta), (6 / 4) * mv_sigma(d, theta))
#'
#' # For nu <= 2 the covariance is infinite and the density is finite; for
#' # nu <= 1 the mean does not exist.
#' theta$nu <- 1.5
#' c(density = distrib_pdf(d, c(0, 0), theta), variance = variance(d, theta)[1, 1])
#' theta$nu <- 0.8
#' c(density = distrib_pdf(d, c(0, 0), theta), mean = mean(d, theta)[1])
#'
#' # At four units of scale from the location, a t at nu = 6 has a far larger
#' # density than a gaussian with the same matrix.
#' theta$nu <- 6
#' g <- mvgaussian1_distrib(2)
#' c(t = distrib_pdf(d, c(4.5, -4.3), theta),
#'   gaussian = distrib_pdf(g, c(4.5, -4.3), theta[1:5]))
#'
#' # The other parametrization carries the inverse scale matrix, and the prefix
#' # of its free values names it.
#' o <- mvstudent_t2_distrib(2)
#' o@params
#'
#' # Its matrix is the precision of the response only up to (nu - 2) / nu.
#' th_o <- list(mu1 = 0, mu2 = 0, omega_log_L1 = 0.1,
#'              omega_log_L2 = -0.2, omega_L2.1 = 0.4, nu = 6)
#' Om <- parameters7::param_value(o@param, unlist(th_o)[3:5])
#' all.equal(unname(solve(variance(o, th_o))),
#'           unname((6 - 2) / 6 * unclass(Om)))
#'
#' # The correlations are the response's exactly: a positive multiple cancels
#' # out of a ratio.
#' all.equal(cov2cor(mv_sigma(o, th_o)), cov2cor(variance(o, th_o)))
#'
#' @export
mvstudent_t1_distrib <- function(n_dim, sigma = NULL,
                                 link_nu = linkfunctions7::log_link()) {
  mvstudent_t_build(n_dim, sigma, link_nu, inverted = FALSE,
                    cls = MvStudentT1Distrib)
}

#' @rdname mvstudent_t1_distrib
#' @export
mvstudent_t2_distrib <- function(n_dim, omega = NULL,
                                 link_nu = linkfunctions7::log_link()) {
  mvstudent_t_build(n_dim, omega, link_nu, inverted = TRUE,
                    cls = MvStudentT2Distrib)
}


#' @title Build One of the Two Multivariate Student t Families
#'
#' @description
#' The body that [mvstudent_t1_distrib()] and [mvstudent_t2_distrib()] share.
#' The two differ in the class that they instantiate, in the `inverted`
#' property, in the prefix of their free names, in the tag of `distrib_name`
#' (`sigma=` or `omega=`) and in what `params_interpretation` calls the matrix.
#' Every method is registered on [MvStudentTDistrib], the parent of both.
#'
#' @param n_dim The dimension of one observation.
#' @param s The matrix parametrization, or `NULL` for
#'   `parameters7::log_cholesky(n_dim)`.
#' @param link_nu The link on the degrees of freedom.
#' @param inverted `TRUE` where `s` carries the inverse scale matrix.
#' @param cls The S7 class to instantiate.
#'
#' @return An object of class `cls`.
#'
#' @seealso [mvstudent_t1_distrib()], the constructors that it serves.
#'
#' @keywords internal
mvstudent_t_build <- function(n_dim, s, link_nu, inverted, cls) {
  if (!is.numeric(n_dim) || length(n_dim) != 1L || !is.finite(n_dim) ||
    n_dim < 1 || n_dim != round(n_dim)) {
    stop("'n_dim' must be a single positive integer.", call. = FALSE)
  }
  p <- as.integer(n_dim)
  if (is.null(s)) s <- parameters7::log_cholesky(p)

  if (!S7::S7_inherits(s, parameters7::parameter)) {
    stop("The matrix parameter must be a parameters7 'parameter' object.", call. = FALSE)
  }
  if (s@dimension != p) {
    stop(sprintf(
      "The matrix parameter has dimension %d but the distribution has dimension %d.",
      s@dimension, p
    ), call. = FALSE)
  }
  if (s@rank < s@dimension) {
    stop(sprintf(paste0(
      "The matrix parameter is rank deficient (%d of %d), so it does not describe a\n",
      "  density: the law would be supported on a subspace. Such a structure is\n",
      "  a penalty, not a distribution."
    ), s@rank, s@dimension), call. = FALSE)
  }
  if (!S7::S7_inherits(link_nu, linkfunctions7::link)) {
    stop("'link_nu' must be a linkfunctions7 link object.", call. = FALSE)
  }

  mu_names <- paste0("mu", seq_len(p))
  # The scale matrix of a t is written Sigma, so it takes the same prefix as a
  # covariance, and its inverse takes the precision's. Neither is a moment of
  # the response: what the prefix says is which of the two matrices a free
  # value builds, and the constructor's page says what those matrices are.
  free_names <- mv_prefixed_names(s@free_names, inverted)
  clash <- intersect(c(mu_names, "nu"), free_names)
  if (length(clash)) {
    stop(sprintf(paste0(
      "The matrix parameter's free value '%s' collides with a parameter name.\n",
      "  Every derivative component is keyed by these, so they must be unique."
    ), clash[1L]), call. = FALSE)
  }

  params <- c(mu_names, free_names, "nu")
  n_par <- length(params)
  ident <- linkfunctions7::identity_link()

  cls(
    distrib_name = sprintf(
      "multivariate student t [%dd, %s=%s]", p,
      if (inverted) "omega" else "sigma", s@param_name
    ),
    dimension = "multivariate",
    n_dim = p,
    bounds = c(-Inf, Inf),
    params = params,
    params_interpretation = stats::setNames(
      c(rep("location", p),
        rep(if (inverted) "precision" else "scale", s@n_free),
        "degrees of freedom"),
      params
    ),
    n_params = n_par,
    params_bounds = stats::setNames(
      c(rep(list(c(-Inf, Inf)), n_par - 1L), list(c(0, Inf))), params
    ),
    link_params = stats::setNames(
      c(rep(list(ident), n_par - 1L), list(link_nu)), params
    ),
    param = s,
    inverted = inverted
  )
}


#' @title The Pieces From Which a Multivariate t Evaluates
#'
#' @description
#' Assembles, once per call, the location, the scale matrix, its inverse, the
#' log-determinant and the degrees of freedom from a flat parameter vector,
#' together with the derivative arrays of \eqn{\Sigma} and of \eqn{\Sigma^{-1}}
#' when they are requested. Every method of [MvStudentTDistrib] that evaluates
#' the density or its derivatives calls this first and works from the result,
#' so a parameter vector is unpacked and a matrix factorized once per call and
#' not once per component.
#'
#' @details
#' The matrix quantities come from [mv_matrix_pieces()], which the gaussian
#' reads as well: where `distrib@inverted` is `TRUE` the parametrization
#' carries \eqn{\Sigma^{-1}}, and the helper returns \eqn{\Sigma} and the
#' derivatives of both matrices whichever side it carries.
#'
#' @param distrib An [MvStudentTDistrib] object.
#' @param theta A named list of parameters, already aligned by the generic, or
#'   any list whose components are in `distrib@params` order.
#' @param derivs Logical of length 1. When `TRUE` the first derivative arrays
#'   are computed. Defaults to `FALSE`.
#' @param derivs2 Logical of length 1. When `TRUE` the second derivative arrays
#'   are computed, and the first ones with them. Defaults to `FALSE`.
#'
#' @return A named list with `mu` (numeric of length \eqn{p}), `eta` (the
#'   matrix parametrization's free vector), `nu` (a single number), `p`, `s`
#'   (the parametrization itself) and the elements of [mv_matrix_pieces()]:
#'   `sigma` and `sigma_inv` (\eqn{p \times p} matrices), `logdet` (the
#'   log-determinant of \eqn{\Sigma}) and `inverted`; with derivatives
#'   requested also `a` and `ai` (the first derivative arrays of \eqn{\Sigma}
#'   and of \eqn{\Sigma^{-1}}) and `dlogdet`; with `derivs2` also `a2`, `ai2`
#'   and `d2logdet`.
#'
#' @section Notation:
#' \eqn{\mu} is the location, \eqn{\Sigma} the scale matrix, \eqn{\nu} the
#' degrees of freedom and \eqn{\eta} the free vector of the matrix
#' parametrization.
#'
#' @seealso [mvstudent_t1_distrib()] for the family, [mvt_weights()] for the
#'   quantities computed from these, and [mvg_pieces()] for the gaussian's
#'   twin.
#'
#' @examples
#' d <- mvstudent_t1_distrib(2)
#' theta <- list(mu1 = 0.5, mu2 = -0.3, sigma_log_L1 = 0.1,
#'               sigma_log_L2 = -0.2, sigma_L2.1 = 0.4, nu = 6)
#' pc <- distributions7:::mvt_pieces(d, theta, derivs = TRUE)
#' names(pc)
#'
#' # sigma_inv is the inverse of the scale matrix, and logdet its determinant.
#' all.equal(pc$sigma %*% pc$sigma_inv, diag(2), check.attributes = FALSE)
#' all.equal(pc$logdet, log(det(pc$sigma)))
#'
#' # nu is carried apart from the matrix, as its own number.
#' c(nu = pc$nu, n_free = length(pc$eta), n_deriv = length(pc$a))
#'
#' @keywords internal
mvt_pieces <- function(distrib, theta, derivs = FALSE, derivs2 = FALSE) {
  p <- distrib@n_dim
  s <- distrib@param
  v <- mv_flat_theta(distrib, theta)
  mu <- unname(v[seq_len(p)])
  eta <- unname(v[p + seq_len(s@n_free)])
  nu <- unname(v[[length(v)]])

  # the scale matrix, its inverse and their arrays, transported once where the
  # parametrization carries the inverse; the same helper the gaussian reads,
  # so the two families share the arithmetic rather than a copy of it
  c(list(mu = mu, eta = eta, nu = nu, p = p, s = s),
    mv_matrix_pieces(s, eta, distrib@inverted, derivs, derivs2))
}


#' @title Multivariate Student t Density
#' @name distrib_pdf.MvStudentTDistrib
#'
#' @description
#' Computes the log-density
#' \deqn{\ell = \log\Gamma\!\left(\tfrac{\nu+p}{2}\right)
#'   - \log\Gamma\!\left(\tfrac{\nu}{2}\right)
#'   - \tfrac{p}{2}\log(\nu\pi) - \tfrac{1}{2}\log\lvert\Sigma\rvert
#'   - \tfrac{\nu+p}{2}\log\!\left(1 + \tfrac{q}{\nu}\right),}
#' with \eqn{q = (y-\mu)^\top \Sigma^{-1}(y-\mu)}, row by row, and exponentiates
#' it unless `log = TRUE`. The gamma functions and \eqn{-\tfrac{p}{2}\log\nu}
#' enter as the single quantity of [mvt_cp_d()], whose direct form cancels as
#' \eqn{\nu} grows. The last logarithm is evaluated with [base::log1p()], which
#' keeps its relative accuracy when \eqn{q/\nu} is small, and \eqn{q} is the
#' squared norm of the residual whitened by the parametrization's factor (see
#' [mv_quad()]).
#'
#' @param distrib An [MvStudentTDistrib] object, from
#'   [mvstudent_t1_distrib()].
#' @param y An \eqn{n \times p} numeric matrix of observations, one row each.
#'   A plain numeric vector of length \eqn{p} is read as a single observation.
#'   Every point of \eqn{\mathbb{R}^p} is in the support, so no row is
#'   rejected. A matrix with zero rows returns `numeric(0)`.
#' @param theta A named list of parameters, each component a single number:
#'   `mu1`, ..., `mup`, the matrix parametrization's prefixed free values, and
#'   `nu`, which must be strictly positive.
#' @param log Logical of length 1. When `TRUE` the log-density is returned,
#'   which stays finite far into a tail. Defaults to `FALSE`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of length \eqn{n}, one density per row of `y`.
#'
#' @section Notation:
#' \eqn{\mu} is the location, \eqn{\Sigma} the scale matrix, \eqn{\nu} the
#' degrees of freedom, \eqn{p} the dimension and \eqn{q} the squared
#' Mahalanobis distance of an observation from the location.
#'
#' @seealso [distrib_gradient.MvStudentTDistrib()] for the derivatives,
#'   [distrib_rng.MvStudentTDistrib()] to draw from it,
#'   [distrib_pdf.MvGaussianDistrib()] for the limiting family, and
#'   [distrib_pdf()] for the generic.
#'
#' @examples
#' d <- mvstudent_t1_distrib(2)
#' theta <- list(mu1 = 0.5, mu2 = -0.3, sigma_log_L1 = 0.1,
#'               sigma_log_L2 = -0.2, sigma_L2.1 = 0.4, nu = 6)
#' y <- rbind(c(0, 0), c(1, -1), c(0.5, -0.3))
#'
#' distrib_pdf(d, y, theta)
#'
#' # Against the formula written out.
#' S <- mv_sigma(d, theta)
#' q <- mahalanobis(y, c(0.5, -0.3), S)
#' ref <- lgamma(4) - lgamma(3) - log(6 * pi) - 0.5 * log(det(S)) -
#'   4 * log1p(q / 6)
#' all.equal(distrib_pdf(d, y, theta, log = TRUE), ref)
#'
#' # As nu grows the density approaches the gaussian's at the same matrix.
#' g <- mvgaussian1_distrib(2)
#' vapply(c(6, 60, 6000), function(nu) {
#'   t2 <- theta; t2$nu <- nu
#'   distrib_pdf(d, c(1, -1), t2, log = TRUE)
#' }, numeric(1))
#' distrib_pdf(g, c(1, -1), theta[1:5], log = TRUE)
#'
#' @keywords internal
S7::method(distrib_pdf, MvStudentTDistrib) <- function(distrib, y, theta,
                                                       log = FALSE, ...) {
  y <- as_mv_matrix(distrib, y)
  pc <- mvt_pieces(distrib, theta)
  if (!nrow(y)) return(numeric(0))
  r <- sweep(y, 2L, pc$mu, "-")
  q <- mv_quad(r, pc$s, pc$eta, pc$inverted)
  nu <- pc$nu
  p <- pc$p
  # lgamma((nu + p)/2) - lgamma(nu/2) - (p/2) log(nu) through mvt_cp_d(): the
  # direct difference of two lgamma values cancels as nu grows
  ld <- mvt_cp_d(nu, p, 0L) - (p / 2) * base::log(pi) - 0.5 * pc$logdet -
    ((nu + p) / 2) * base::log1p(q / nu)
  if (log) ld else exp(ld)
}


#' @title Multivariate Student t Generator
#' @name distrib_rng.MvStudentTDistrib
#'
#' @description
#' Draws from the family through its scale-mixture representation, \eqn{y = \mu
#' + L z\sqrt{\nu/g}} with \eqn{z} a vector of independent standard normals,
#' \eqn{g \sim \chi^2_\nu} independent of them, and \eqn{LL^\top = \Sigma}. The
#' factor \eqn{L} is read from the parametrization (see [mv_sampling_factor()]),
#' so a scale matrix that is nearly singular is sampled without a decomposition.
#' A \eqn{t} is a gaussian whose precision has been multiplied by a gamma
#' variate, and the same representation gives the closed-form expected
#' information.
#'
#' The draws consume `n * p` normal variates and `n` chi-squared variates from
#' R's own generators, in that order, so the stream is reproducible under
#' [base::set.seed()].
#'
#' @param distrib An [MvStudentTDistrib] object, from
#'   [mvstudent_t1_distrib()].
#' @param n The number of observations to draw. A single non-negative whole
#'   number.
#' @param theta A named list of parameters, each component a single number.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return An \eqn{n \times p} numeric matrix, one draw per row, with column
#'   names `v1`, ..., `vp`.
#'
#' @section Notation:
#' \eqn{\mu} is the location, \eqn{\Sigma} the scale matrix, \eqn{\nu} the
#' degrees of freedom and \eqn{L} a factor of \eqn{\Sigma}, \eqn{LL^\top =
#' \Sigma}.
#'
#' @seealso [distrib_pdf.MvStudentTDistrib()] for the density this draws from,
#'   [variance.MvStudentTDistrib()] for the covariance the sample approaches
#'   when it exists, and [distrib_rng()] for the generic.
#'
#' @examples
#' d <- mvstudent_t1_distrib(2)
#' theta <- list(mu1 = 0.5, mu2 = -0.3, sigma_log_L1 = 0.1,
#'               sigma_log_L2 = -0.2, sigma_L2.1 = 0.4, nu = 6)
#'
#' set.seed(1)
#' distrib_rng(d, 4, theta)
#'
#' # At nu = 6 the covariance exists and the sample approaches it.
#' set.seed(2)
#' big <- distrib_rng(d, 40000, theta)
#' round(var(big), 2)
#' round(variance(d, theta), 2)
#'
#' # At nu = 1.5 the covariance is infinite, and the sample covariance grows
#' # with the sample.
#' t2 <- theta; t2$nu <- 1.5
#' set.seed(3)
#' vapply(c(2000, 20000, 200000),
#'        function(m) var(distrib_rng(d, m, t2))[1, 1], numeric(1))
#'
#' @keywords internal
S7::method(distrib_rng, MvStudentTDistrib) <- function(distrib, n, theta, ...) {
  pc <- mvt_pieces(distrib, theta)
  p <- pc$p
  l <- mv_sampling_factor(pc$s, pc$eta, distrib@inverted)
  z <- matrix(stats::rnorm(n * p), n, p) %*% t(l)
  w <- sqrt(pc$nu / stats::rchisq(n, df = pc$nu))
  out <- sweep(z * w, 2L, pc$mu, "+")
  colnames(out) <- paste0("v", seq_len(p))
  out
}


#' @title The Weight a Multivariate t Gives Each Observation
#'
#' @description
#' Computes the centered response \eqn{r_i = y_i - \mu}, its image
#' \eqn{w_i = \Sigma^{-1}r_i} under the inverse scale matrix, the squared
#' Mahalanobis distance \eqn{q_i = r_i^\top w_i}, and the weight
#' \eqn{c_i = (\nu + p)/(\nu + q_i)}. Every derivative of the family is written
#' in those four, so they are formed once per call.
#'
#' @details
#' The weight limits the influence of a distant observation. At \eqn{q = 0} it
#' is \eqn{(\nu+p)/\nu}, and it decays like \eqn{1/q}, so an observation far
#' from the location contributes less to every derivative. Letting \eqn{\nu \to
#' \infty} sends it to one and recovers the gaussian.
#'
#' @param y An \eqn{n \times p} numeric matrix of observations, already
#'   coerced by [as_mv_matrix()].
#' @param pc The result of [mvt_pieces()], from which `mu`, `sigma_inv`, `nu`
#'   and `p` are read.
#'
#' @return A named list with `r` and `w`, each an \eqn{n \times p} numeric
#'   matrix, and `q` and `cw`, each a numeric vector of length \eqn{n}.
#'
#' @section Notation:
#' \eqn{\mu} is the location, \eqn{\Sigma} the scale matrix, \eqn{\nu} the
#' degrees of freedom, \eqn{p} the dimension, \eqn{q} the squared Mahalanobis
#' distance and \eqn{c} the weight.
#'
#' @seealso [mvt_pieces()] for the argument and
#'   [distrib_gradient.MvStudentTDistrib()] for the first consumer.
#'
#' @examples
#' d <- mvstudent_t1_distrib(2)
#' theta <- list(mu1 = 0, mu2 = 0, sigma_log_L1 = 0,
#'               sigma_log_L2 = 0, sigma_L2.1 = 0, nu = 6)
#' pc <- distributions7:::mvt_pieces(d, theta)
#'
#' # At the identity scale matrix q is the squared distance from the origin.
#' y <- rbind(c(0, 0), c(1, 0), c(3, 0), c(10, 0))
#' z <- distributions7:::mvt_weights(y, pc)
#' rbind(q = z$q, weight = z$cw)
#'
#' # The weight starts at (nu + p) / nu and decays like 1 / q.
#' c(at_zero = (6 + 2) / 6, computed = z$cw[1])
#'
#' # As nu grows it flattens to one, and the family stops downweighting.
#' vapply(c(6, 60, 6000), function(nu) {
#'   t2 <- theta; t2$nu <- nu
#'   distributions7:::mvt_weights(y, distributions7:::mvt_pieces(d, t2))$cw[4]
#' }, numeric(1))
#'
#' @keywords internal
mvt_weights <- function(y, pc) {
  r <- sweep(y, 2L, pc$mu, "-")
  w <- r %*% pc$sigma_inv
  q <- mv_quad(r, pc$s, pc$eta, pc$inverted)
  list(r = r, w = w, q = q, cw = (pc$nu + pc$p) / (pc$nu + q))
}


#' @title Multivariate Student t Score
#' @name distrib_gradient.MvStudentTDistrib
#'
#' @description
#' Computes the first derivatives of the log-density in closed form. With \eqn{w
#' = \Sigma^{-1}(y-\mu)}, \eqn{q = (y-\mu)^\top w}, \eqn{c = (\nu+p)/(\nu+q)}
#' and \eqn{P_k} the derivative of \eqn{\Sigma^{-1}} in the \eqn{k}-th free
#' value, \deqn{\partial_\mu \ell = c\,w, \qquad \partial_{\eta_k}\ell =
#' -\tfrac{1}{2}\partial_{\eta_k}\log\lvert\Sigma\rvert - \tfrac{c}{2}\,
#' (y-\mu)^\top P_k (y-\mu),} \deqn{\partial_\nu \ell = \tfrac{1}{2}\left[
#' \psi\!\left(\tfrac{\nu+p}{2}\right) - \psi\!\left(\tfrac{\nu}{2}\right) -
#' \tfrac{p}{\nu} - \log\!\left(1+\tfrac{q}{\nu}\right) +
#' \tfrac{(\nu+p)q}{\nu(\nu+q)}\right].} The gaussian's score is the limit
#' \eqn{c \to 1}. Each observation enters the location and matrix components
#' through the weight \eqn{c}, which decreases as \eqn{q} grows. \eqn{P_k} is
#' read from [mv_matrix_pieces()] and \eqn{q} is computed by [mv_quad()], both
#' without the cancellation of an explicit inverse where the scale matrix is
#' nearly singular.
#'
#' @details
#' The \eqn{\nu} component is not evaluated in the form printed above. Both of
#' its bracketed pairs cancel to leading order as \eqn{\nu} grows, so it is
#' assembled as \eqn{A_p(\nu) + D(q/\nu) + (p/\nu)u/(1+u)}, with [mvt_A()] and
#' [mvt_D()] carrying out each cancellation exactly.
#'
#' @param distrib An [MvStudentTDistrib] object, from
#'   [mvstudent_t1_distrib()].
#' @param y An \eqn{n \times p} numeric matrix of observations. A vector of
#'   length \eqn{p} is read as a single observation.
#' @param theta A named list of parameters, each component a single number.
#' @param scale One of `"parameter"` (the default) or `"link"`. The method
#'   returns the parameter scale and the generic applies the link scale to its
#'   result. The two differ here: under the default log link the `nu` component
#'   is multiplied by \eqn{\nu}, and the rest are unchanged.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A named list with one numeric vector of length \eqn{n} per
#'   parameter, in `distrib@params` order: \eqn{p} location components, one per
#'   free value of the matrix parametrization, then `nu`.
#'
#' @section Notation:
#' \eqn{\mu} is the location, \eqn{\Sigma} the scale matrix, \eqn{\nu} the
#' degrees of freedom, \eqn{\eta} the free vector of the matrix
#' parametrization, \eqn{P_k = \partial\Sigma^{-1}/\partial\eta_k}, \eqn{q} the
#' squared Mahalanobis distance, \eqn{c} the weight, \eqn{\psi} the digamma
#' function and \eqn{\ell} the log-density of one observation.
#'
#' @seealso [distrib_hessian.MvStudentTDistrib()] for the observed curvature,
#'   [distrib_expected_hessian.MvStudentTDistrib()] for the information,
#'   [mvt_weights()] for the weight, [mvt_A()] for the cancellation, and
#'   [distrib_gradient()] for the generic.
#'
#' @examples
#' d <- mvstudent_t1_distrib(2)
#' theta <- list(mu1 = 0.5, mu2 = -0.3, sigma_log_L1 = 0.1,
#'               sigma_log_L2 = -0.2, sigma_L2.1 = 0.4, nu = 6)
#' set.seed(1)
#' y <- distrib_rng(d, 25, theta)
#'
#' g <- distrib_gradient(d, y, theta)
#' vapply(g, sum, numeric(1))
#'
#' # Against a numerical derivative of the log-likelihood.
#' ll <- function(v) {
#'   t2 <- as.list(v); names(t2) <- d@params
#'   sum(distrib_pdf(d, y, t2, log = TRUE))
#' }
#' max(abs(vapply(g, sum, numeric(1)) - numDeriv::grad(ll, unlist(theta))))
#'
#' # The location component is c w, the weight times the residual mapped
#' # through the inverse scale matrix.
#' z <- distributions7:::mvt_weights(y, distributions7:::mvt_pieces(d, theta))
#' all.equal(cbind(g$mu1, g$mu2), z$cw * z$w, check.attributes = FALSE)
#'
#' # The nu component decays as nu^-2.
#' vapply(c(1e3, 1e6, 1e10), function(nu) {
#'   t2 <- theta; t2$nu <- nu
#'   distrib_gradient(d, y[1, ], t2)$nu
#' }, numeric(1))
#'
#' # And the link scale differs, this family having one link that is not the
#' # identity.
#' c(parameter = sum(g$nu),
#'   link = sum(distrib_gradient(d, y, theta, scale = "link")$nu))
#'
#' @keywords internal
S7::method(distrib_gradient, MvStudentTDistrib) <- function(distrib, y, theta,
                                                            scale = c("parameter", "link"),
                                                            ...) {
  y <- as_mv_matrix(distrib, y)
  pc <- mvt_pieces(distrib, theta, derivs = TRUE)
  z <- mvt_weights(y, pc)
  nu <- pc$nu
  p <- pc$p

  out <- vector("list", distrib@n_params)
  names(out) <- distrib@params
  for (j in seq_len(p)) out[[j]] <- z$cw * z$w[, j]

  # w' A_k w = -r' d_k(Sigma^-1) r, read through the chart's derivative of the
  # inverse: the sandwich Sigma^-1 A_k Sigma^-1 cancels near a singular matrix
  dld <- pc$dlogdet
  for (k in seq_along(pc$ai)) {
    out[[p + k]] <- -0.5 * dld[[k]] -
      0.5 * z$cw * rowSums((z$r %*% pc$ai[[k]]) * z$r)
  }

  # A_p(nu) + D(u) + (p/nu) u/(1+u), u = q/nu: the two cancellations written
  # out.  See mvt_A() for what the direct form costs at large nu.
  u_nu <- z$q / nu
  out[[distrib@n_params]] <- 0.5 * (
    mvt_A(nu, p) + mvt_D(u_nu) + (p / nu) * u_nu / (1 + u_nu)
  )
  out
}


#' @title The Digamma Difference of a Multivariate t, Without the Cancellation
#'
#' @description
#' Computes
#' \eqn{A_p(\nu) = \psi\!\left(\tfrac{\nu+p}{2}\right)
#' - \psi\!\left(\tfrac{\nu}{2}\right) - \tfrac{p}{\nu}}, the part of the score
#' in \eqn{\nu} that does not involve the data, as twice the first derivative
#' of [mvt_cp_d()].
#'
#' @details
#' As \eqn{\nu} grows the multivariate t tends to the multivariate gaussian and
#' every derivative in \eqn{\nu} vanishes, so \eqn{A_p} is a difference of
#' terms agreeing to leading order:
#' \eqn{\psi(\tfrac{\nu+p}{2}) - \psi(\tfrac{\nu}{2})} is close to \eqn{p/\nu},
#' the quantity from which it is subtracted. In floating point the direct form
#' loses its significant digits as \eqn{\nu} grows and eventually has the wrong
#' sign. For even \eqn{p} the recurrence \eqn{\psi(x+1) = \psi(x) + 1/x} gives
#' a sum whose terms carry one sign,
#' \deqn{A_p(\nu) = -\sum_{j=0}^{p/2-1} \frac{4j}{\nu(\nu+2j)},}
#' exactly zero at \eqn{p = 2}; for odd \eqn{p} [mvt_cp_d()] sums an
#' asymptotic series above \eqn{\nu = \max(20, 4p)}.
#'
#' @param nu The degrees of freedom, a numeric vector of positive values. Not
#'   validated.
#' @param p The dimension, a single positive whole number.
#'
#' @return A numeric vector as long as `nu`: positive at \eqn{p = 1}, exactly
#'   zero at \eqn{p = 2}, and negative for \eqn{p \ge 3}.
#'
#' @section Notation:
#' \eqn{\nu} is the degrees of freedom, \eqn{p} the dimension and \eqn{\psi}
#' the digamma function.
#'
#' @seealso [mvt_cp_d()] for the computation, [mvt_D()] for the part that
#'   carries the data, and [distrib_gradient.MvStudentTDistrib()] for the
#'   consumer.
#'
#' @examples
#' # At p = 2 the quantity is exactly zero and the direct form is noise.
#' direct <- function(nu, p) digamma((nu + p) / 2) - digamma(nu / 2) - p / nu
#' rbind(exact = distributions7:::mvt_A(c(3, 10, 1e6), 2),
#'       direct = direct(c(3, 10, 1e6), 2))
#'
#' # At p = 4 the two agree while the direct form has digits left, and part
#' # company at large nu.
#' nu <- c(3, 1e3, 1e6, 1e8, 1e9)
#' signif(rbind(exact = distributions7:::mvt_A(nu, 4),
#'              direct = direct(nu, 4)), 5)
#'
#' # The exact form decays as -p(p - 2) / (2 nu^2), which is -4 / nu^2 at
#' # p = 4.
#' cbind(nu = nu, scaled = distributions7:::mvt_A(nu, 4) * nu^2)
#'
#' @keywords internal
mvt_A <- function(nu, p) {
  2 * mvt_cp_d(nu, p, 1L)
}


# Commented out on 2026-10-12 (T4): mvt_A() is twice mvt_cp_d(nu, p, 1), whose series is not cut at three
# terms; the direct branch of this one lost up to 1e-10 below its crossover
# (3.5e-11 at nu = 150).
# Kept rather than deleted.
# @title The Univariate Digamma Difference at One Degree of Freedom
#
# @description
# Computes
# \eqn{A_1(\nu) = \psi\!\left(\tfrac{\nu+1}{2}\right)
# - \psi\!\left(\tfrac{\nu}{2}\right) - \tfrac{1}{\nu}}, the odd-dimension
# base case of [mvt_A()]. The half-integer shift here leaves no recurrence to
# telescope, so above \eqn{\nu = 200} the expansion
# \eqn{1/(2\nu^2) - 1/(4\nu^4) + 1/(2\nu^6)} is used, and below it the direct
# difference, which still has digits there.
#
# @details
# The crossover of 200 is the one measured for the same expansion in this
# package's `src/student_t.cpp`, and this is the one place in the package
# where that series exists twice. The two copies are pinned against each other
# in the tests. Across the switch the two branches agree to eleven figures:
# at \eqn{\nu = 199, 200, 201} the values run
# \eqn{1.2626\times10^{-5}}, \eqn{1.2500\times10^{-5}},
# \eqn{1.2376\times10^{-5}} with no step.
#
# @param nu The degrees of freedom, a numeric vector of any length, strictly
#   positive. The branch is taken elementwise.
#
# @return A numeric vector as long as `nu`, positive and decaying as
#   \eqn{1/(2\nu^2)}.
#
# @section Notation:
# \eqn{\nu} is the degrees of freedom and \eqn{\psi} the digamma function.
#
# @seealso [mvt_A()], which calls this on an odd dimension, and [mvt_T1()] for
#   the trigamma twin.
#
# @examples
# # No step across the crossover at nu = 200.
# signif(distributions7:::mvt_A1(c(199, 200, 201)), 10)
#
# # Below it the branch is the direct difference and agrees exactly.
# direct <- function(nu) digamma((nu + 1) / 2) - digamma(nu / 2) - 1 / nu
# identical(distributions7:::mvt_A1(c(3, 50)), direct(c(3, 50)))
#
# # Above it the series is the accurate one: scaled by 2 nu^2 both approach 1,
# # and the direct form starts to wander.
# nu <- c(1e3, 1e4, 1e5)
# rbind(series = distributions7:::mvt_A1(nu) * 2 * nu^2,
#       direct = direct(nu) * 2 * nu^2)
#
# @keywords internal
# mvt_A1 <- function(nu) {
#   # psi((nu+1)/2) - psi(nu/2) - 1/nu = 1/(2 nu^2) - 1/(4 nu^4) + 1/(2 nu^6)
#   # above the crossover measured for student_t.cpp's own copy of it
#   out <- numeric(length(nu))
#   big <- nu >= 200
#   if (any(big)) {
#     u2 <- 1 / nu[big]^2
#     out[big] <- u2 * (0.5 + u2 * (-0.25 + u2 * 0.5))
#   }
#   if (any(!big)) {
#     v <- nu[!big]
#     out[!big] <- digamma((v + 1) / 2) - digamma(v / 2) - 1 / v
#   }
#   out
# }


# Commented out on 2026-10-12 (T4): the observed Hessian takes D^2 c_p from mvt_cp_d(nu, p, 2) directly, so
# this has no caller.
# Kept rather than deleted.
# @title The Trigamma Difference of a Multivariate t, Without the Cancellation
#
# @description
# Computes
# \eqn{T_p(\nu) = \tfrac12\left[\psi'\!\left(\tfrac{\nu+p}{2}\right)
# - \psi'\!\left(\tfrac{\nu}{2}\right)\right] + \tfrac{p}{\nu^2}}, the part of
# the curvature in \eqn{\nu} that does not involve the data, written so that
# its terms carry one sign. It is the second-order twin of [mvt_A()] and
# cancels for the same reason: every derivative in \eqn{\nu} vanishes as the
# family approaches the gaussian.
#
# @details
# For even \eqn{p} the recurrence \eqn{\psi'(x+1) = \psi'(x) - 1/x^2}
# telescopes into
# \deqn{T_p(\nu) = \sum_{j=0}^{p/2-1} \frac{8j(\nu+j)}{\nu^2(\nu+2j)^2},}
# exactly zero at \eqn{p = 2} where the direct form returns noise. For odd
# \eqn{p} the half-integer shift carries the quantity onto the univariate
# [mvt_T1()]. Measured at \eqn{p = 4}, the direct form is wrong by a relative
# \eqn{4.4\times10^{-5}} at \eqn{\nu = 10^6} and by a factor of five at
# \eqn{10^8}.
#
# @param nu The degrees of freedom, a numeric vector of any length, strictly
#   positive. Nothing is validated.
# @param p The dimension, a single positive whole number. Even and odd take
#   different branches.
#
# @return A numeric vector as long as `nu`, non-negative at every \eqn{p} and
#   exactly zero at \eqn{p = 2}.
#
# @section Notation:
# \eqn{\nu} is the degrees of freedom, \eqn{p} the dimension and \eqn{\psi'}
# the trigamma function.
#
# @seealso [mvt_A()] for the first-order twin, [mvt_T1()] for the univariate
#   base case, and [distrib_hessian.MvStudentTDistrib()] for the consumer.
#
# @examples
# # Exactly zero at p = 2, where the direct form is noise.
# direct <- function(nu, p)
#   0.5 * (trigamma((nu + p) / 2) - trigamma(nu / 2)) + p / nu^2
# rbind(exact = distributions7:::mvt_T(c(3, 10, 1e6), 2),
#       direct = direct(c(3, 10, 1e6), 2))
#
# # At p = 4 the two agree until the cancellation bites.
# nu <- c(3, 1e4, 1e6, 1e8)
# signif(rbind(exact = distributions7:::mvt_T(nu, 4),
#              direct = direct(nu, 4)), 5)
#
# @keywords internal
# mvt_T <- function(nu, p) {
#   2 * mvt_cp_d(nu, p, 2L)
# }


# Commented out on 2026-10-12 (T4): mvt_cp_d() replaces it; its three-term series was off by 1.7e-11 at
# nu = 100.
# Kept rather than deleted.
# @title The Univariate Trigamma Difference at One Degree of Freedom
#
# @description
# Computes
# \eqn{T_1(\nu) = \tfrac12\left[\psi'\!\left(\tfrac{\nu+1}{2}\right)
# - \psi'\!\left(\tfrac{\nu}{2}\right)\right] + \tfrac{1}{\nu^2}}, the
# odd-dimension base case of [mvt_T()]. Above \eqn{\nu = 100} it uses the
# expansion \eqn{-1/\nu^3 + 1/\nu^5 - 3/\nu^7}, and below it the direct
# difference.
#
# @details
# The series comes from the same duplication formula this package's
# `src/student_t.cpp` uses, halved and with the \eqn{1/\nu^2} folded in.
# Across the switch the two branches agree: at \eqn{\nu = 99, 100, 101} the
# values run \eqn{-1.0305\times10^{-6}}, \eqn{-9.9990\times10^{-7}},
# \eqn{-9.7049\times10^{-7}} with no step.
#
# Note the sign. \eqn{T_1} is NEGATIVE, where \eqn{T_p} for even \eqn{p} is
# non-negative; an odd dimension assembles from this base case plus a
# non-negative sum.
#
# @param nu The degrees of freedom, a numeric vector of any length, strictly
#   positive. The branch is taken elementwise.
#
# @return A numeric vector as long as `nu`, negative and decaying as
#   \eqn{-1/\nu^3}.
#
# @section Notation:
# \eqn{\nu} is the degrees of freedom and \eqn{\psi'} the trigamma function.
#
# @seealso [mvt_T()], which calls this on an odd dimension, and [mvt_A1()] for
#   the digamma twin.
#
# @examples
# # No step across the crossover at nu = 100.
# signif(distributions7:::mvt_T1(c(99, 100, 101)), 8)
#
# # Scaled by nu^3 the series approaches -1.
# nu <- c(1e2, 1e3, 1e4)
# distributions7:::mvt_T1(nu) * nu^3
#
# @keywords internal
# mvt_T1 <- function(nu) {
#   # [psi'((nu+1)/2) - psi'(nu/2)]/2 + 1/nu^2 = -1/nu^3 + 1/nu^5 - 3/nu^7 + ...
#   # from the same duplication student_t.cpp's t_S() uses, halved and with the
#   # 1/nu^2 folded in
#   out <- numeric(length(nu))
#   big <- nu >= 100
#   if (any(big)) {
#     u <- 1 / nu[big]; u2 <- u * u
#     out[big] <- u2 * u * (-1 + u2 * (1 - u2 * 3))
#   }
#   if (any(!big)) {
#     v <- nu[!big]
#     out[!big] <- 0.5 * (trigamma((v + 1) / 2) - trigamma(v / 2)) + 1 / v^2
#   }
#   out
# }


#' @title The Log-Ratio Cancellation of a Multivariate t
#'
#' @description
#' Computes \eqn{D(u) = \tfrac{u}{1+u} - \log(1+u)}, the data-carrying part of
#' the score in \eqn{\nu}, taken at \eqn{u = q/\nu}. The two terms agree to
#' first order, so below \eqn{\lvert u\rvert = 10^{-3}} the expansion
#' \eqn{-u^2/2 + 2u^3/3 - 3u^4/4 + 4u^5/5} is used instead of the difference.
#'
#' @details
#' The argument \eqn{u} is the squared Mahalanobis distance divided by the
#' degrees of freedom. It is small when the family is close to the gaussian or
#' when the observation is close to the location, so small values of \eqn{u} are
#' the ordinary case. Evaluated directly at \eqn{u = 10^{-8}}, the difference
#' keeps about eight significant figures.
#'
#' @param u A numeric vector of any length, taken as \eqn{q/\nu}. Values at or
#'   below \eqn{-1} are outside the domain of the logarithm and give `NaN` from
#'   [base::log1p()]; the intended argument is non-negative.
#'
#' @return A numeric vector as long as `u`, non-positive, of order
#'   \eqn{-u^2/2} near zero.
#'
#' @section Notation:
#' \eqn{q} is the squared Mahalanobis distance of an observation from the
#' location and \eqn{\nu} the degrees of freedom.
#'
#' @seealso [mvt_A()] for the part of the score without the data, [mvt_h_d()]
#'   for the higher derivatives of the data term, and
#'   [distrib_gradient.MvStudentTDistrib()] for the consumer.
#'
#' @examples
#' u <- c(1e-8, 1e-5, 1e-3, 0.1, 1)
#'
#' # Near zero the value is -u^2 / 2, and it approaches that ratio smoothly.
#' cbind(u = u, scaled = distributions7:::mvt_D(u) / (-u^2 / 2))
#'
#' # Against the direct difference, relatively. The two agree where the
#' # difference still has digits and part company as u shrinks.
#' cbind(u = u,
#'       rel_gap = abs(distributions7:::mvt_D(u) - (u / (1 + u) - log1p(u))) /
#'         abs(distributions7:::mvt_D(u)))
#'
#' @keywords internal
mvt_D <- function(u) {
  # u/(1+u) - log1p(u) = -u^2/2 + 2u^3/3 - 3u^4/4 + 4u^5/5 - ...
  out <- numeric(length(u))
  sm <- abs(u) < 1e-3
  if (any(sm)) {
    v <- u[sm]
    out[sm] <- v * v * (-0.5 + v * (2 / 3 + v * (-0.75 + v * 0.8)))
  }
  if (any(!sm)) {
    v <- u[!sm]
    out[!sm] <- v / (1 + v) - log1p(v)
  }
  out
}

#' @title The Taylor Coefficients Behind a Multivariate t's Degrees of Freedom
#'
#' @description
#' Returns the Taylor coefficients at \eqn{t = 0} of
#' \eqn{g(t) = (1 - e^{-pt})/(1 - e^{-2t})}, from which every quantity of
#' \eqn{\nu} alone of a multivariate t follows by a Laplace transform. With
#' \eqn{c_p(\nu) = \log\Gamma\{(\nu+p)/2\} - \log\Gamma(\nu/2) - (p/2)\log\nu},
#' \deqn{c_p'(\nu) = \int_0^\infty e^{-\nu t}\{g(t) - p/2\}\,dt,}
#' so that \eqn{c_p'(\nu) = \sum_{n \ge 1} g_n\, n!/\nu^{n+1}}, a series that
#' converges for even \eqn{p} and is asymptotic for odd \eqn{p}.
#'
#' @details
#' The coefficients come from the division of the two power series of
#' \eqn{(1 - e^{-pt})/t} and \eqn{(1 - e^{-2t})/t}, and are kept by dimension
#' in a cache, so each dimension is computed once per session.
#'
#' @param p The dimension, a single positive whole number.
#' @param K The highest power kept. Defaults to 40.
#'
#' @return A numeric vector of length `K + 1`, the coefficients of
#'   \eqn{t^0, \ldots, t^K}; the first is \eqn{p/2}.
#'
#' @seealso [mvt_cp_d()] and [mvt_info_nunu()], the two consumers.
#'
#' @examples
#' g <- distributions7:::mvt_nu_coef(2)
#' # At p = 2 the ratio is identically 1.
#' g[1:4]
#' distributions7:::mvt_nu_coef(3)[1:4]
#'
#' @keywords internal
mvt_nu_coef <- function(p, K = 40L) {
  key <- paste(p, K)
  got <- .mvt_coef_cache[[key]]
  if (!is.null(got)) return(got)
  n <- 0:K
  a <- -(-p)^(n + 1) / factorial(n + 1)
  b <- -(-2)^(n + 1) / factorial(n + 1)
  q <- numeric(K + 1L)
  for (i in seq_len(K + 1L)) {
    acc <- a[i]
    if (i > 1L) acc <- acc - sum(b[2:i] * q[(i - 1L):1L])
    q[i] <- acc / b[1L]
  }
  assign(key, q, envir = .mvt_coef_cache)
  q
}

.mvt_coef_cache <- new.env(parent = emptyenv())


#' @title The Constant of a Multivariate t and Its Derivatives in nu
#'
#' @description
#' Computes the \eqn{k}-th derivative of
#' \eqn{c_p(\nu) = \log\Gamma\{(\nu+p)/2\} - \log\Gamma(\nu/2) -
#' (p/2)\log\nu}, for \eqn{k} from 0 to 4, without the cancellation of the
#' direct form: every derivative vanishes as \eqn{\nu} grows, while the terms
#' of the direct form do not.
#'
#' @details
#' For even \eqn{p} the recurrence of the gamma function telescopes into
#' \deqn{c_p(\nu) = -\tfrac{p}{2}\log 2 + \sum_{j=1}^{p/2-1}
#'   \log\!\left(1 + \tfrac{2j}{\nu}\right),}
#' whose derivatives are
#' \eqn{(-1)^{k-1}(k-1)!\,\nu^{-k}\sum_j \{(1 + 2j/\nu)^{-k} - 1\}}, each term
#' formed with [base::expm1()] and [base::log1p()]. For odd \eqn{p} and
#' \eqn{\nu \ge \max(20, 4p)} the series of [mvt_nu_coef()] is summed,
#' \deqn{c_p^{(k)}(\nu) = (-1)^{k-1}\sum_{n \ge 1} g_n
#'   \frac{(n+k-1)!}{\nu^{n+k}}, \qquad
#'   c_p(\nu) = -\tfrac{p}{2}\log 2 - \sum_{n \ge 1} g_n \frac{(n-1)!}{\nu^n};}
#' below that bound the direct form in [base::lgamma()] and
#' [base::psigamma()] is used, where it keeps its digits.
#'
#' @param nu The degrees of freedom, a numeric vector of positive values.
#' @param p The dimension, a single positive whole number.
#' @param k The order of the derivative, 0 to 4.
#'
#' @return A numeric vector as long as `nu`.
#'
#' @seealso [mvt_h_d()] for the part of the log-density that carries the data,
#'   [mvt_A()], which is twice its first derivative, and [mvt_nu_coef()] for
#'   the series.
#'
#' @examples
#' nu <- c(5, 50, 1e6)
#' # Against the direct form, which loses its digits as nu grows.
#' rbind(exact = distributions7:::mvt_cp_d(nu, 3, 0L),
#'       direct = lgamma((nu + 3) / 2) - lgamma(nu / 2) - 1.5 * log(nu))
#'
#' @keywords internal
mvt_cp_d <- function(nu, p, k) {
  if (p %% 2L == 0L) {
    j <- seq_len(p %/% 2L - 1L)
    if (k == 0L) {
      acc <- rep(-(p / 2) * log(2), length(nu))
      for (jj in j) acc <- acc + log1p(2 * jj / nu)
      return(acc)
    }
    acc <- numeric(length(nu))
    for (jj in j) acc <- acc + expm1(-k * log1p(2 * jj / nu))
    return((-1)^(k - 1L) * factorial(k - 1L) * acc / nu^k)
  }
  out <- numeric(length(nu))
  big <- nu >= max(20, 4 * p)
  if (any(big)) {
    x <- 1 / nu[big]
    g <- mvt_nu_coef(p)
    n <- seq_len(length(g) - 1L)
    cf <- if (k == 0L) {
      -g[n + 1L] * factorial(n - 1L)
    } else {
      (-1)^(k - 1L) * g[n + 1L] * factorial(n + k - 1L)
    }
    s <- numeric(length(x))
    for (i in rev(n)) s <- (s + cf[i]) * x
    out[big] <- if (k == 0L) -(p / 2) * log(2) + s else s * x^k
  }
  if (any(!big)) {
    v <- nu[!big]
    out[!big] <- if (k == 0L) {
      lgamma((v + p) / 2) - lgamma(v / 2) - (p / 2) * log(v)
    } else {
      0.5^k * (psigamma((v + p) / 2, k - 1L) - psigamma(v / 2, k - 1L)) -
        (p / 2) * (-1)^(k - 1L) * factorial(k - 1L) / v^k
    }
  }
  out
}


#' @title The Data Term of a Multivariate t and Its Derivatives in nu
#'
#' @description
#' Computes the \eqn{k}-th derivative in \eqn{\nu} of
#' \eqn{h = \tfrac{\nu+p}{2}\log(1 + q/\nu)}, the part of the log-density that
#' carries the observation, at fixed \eqn{q}. With \eqn{u = q/\nu} and
#' \eqn{a_k(u) = (1+u)^{-k} - 1}, for \eqn{k \ge 2}
#' \deqn{h^{(k)} = \tfrac{1}{2}(-1)^{k-1}(k-1)!\,\nu^{-(k-1)}
#'   \left\{B_k(u) + \tfrac{p}{\nu}\,a_k(u)\right\}, \qquad
#'   B_k(u) = a_k(u) - \tfrac{k}{k-1}\,a_{k-1}(u),}
#' and \eqn{h' = -\tfrac{1}{2}\{D(u) + (p/\nu)\,u/(1+u)\}} with \eqn{D} from
#' [mvt_D()].
#'
#' @details
#' The two terms of \eqn{B_k} agree to first order in \eqn{u}, so below
#' \eqn{u = 0.1} it is summed as the power series
#' \eqn{\sum_{m \ge 2} (-1)^m\{\binom{k+m-1}{m} - \tfrac{k}{k-1}
#' \binom{k+m-2}{m}\}u^m}, to \eqn{m = 22}; above it the direct form loses at
#' most one digit.
#'
#' @param u The ratio \eqn{q/\nu}, a numeric vector.
#' @param nu The degrees of freedom, a single positive number.
#' @param p The dimension.
#' @param k The order of the derivative, 1 to 4.
#'
#' @return A numeric vector as long as `u`.
#'
#' @seealso [mvt_cp_d()] for the part of \eqn{\nu} alone; the derivative of the
#'   log-density in \eqn{\nu} alone is \eqn{c_p^{(k)} - h^{(k)}}.
#'
#' @examples
#' # Against the direct form at nu = 1e9, where its terms cancel.
#' nu <- 1e9; q <- 1.25; p <- 2
#' h <- function(v) (v + p) / 2 * log1p(q / v)
#' distributions7:::mvt_h_d(q / nu, nu, p, 2L)
#' (h(nu * (1 + 1e-4)) - 2 * h(nu) + h(nu * (1 - 1e-4))) / (nu * 1e-4)^2
#'
#' @keywords internal
mvt_h_d <- function(u, nu, p, k) {
  if (k == 1L) return(-0.5 * (mvt_D(u) + (p / nu) * u / (1 + u)))
  ak <- expm1(-k * log1p(u))
  bk <- numeric(length(u))
  sm <- abs(u) < 0.1
  if (any(sm)) {
    m <- 2:22
    cf <- (-1)^m * (choose(k + m - 1, m) - k / (k - 1) * choose(k + m - 2, m))
    v <- u[sm]
    s <- numeric(length(v))
    for (i in rev(seq_along(m))) s <- (s + cf[i]) * v
    bk[sm] <- s * v
  }
  if (any(!sm)) {
    v <- u[!sm]
    bk[!sm] <- expm1(-k * log1p(v)) - k / (k - 1) * expm1(-(k - 1) * log1p(v))
  }
  0.5 * (-1)^(k - 1L) * factorial(k - 1L) / nu^(k - 1L) * (bk + (p / nu) * ak)
}


#' @title The Expected Information in nu of a Multivariate t
#'
#' @description
#' Computes
#' \deqn{I_{\nu\nu} = \tfrac{1}{4}\left\{\psi'\!\left(\tfrac{\nu}{2}\right) -
#'   \psi'\!\left(\tfrac{\nu+p}{2}\right)\right\} -
#'   \frac{p(\nu+p+4)}{2\nu(\nu+p)(\nu+p+2)}}
#' without the cancellation of the direct form, whose terms are of order
#' \eqn{\nu^{-2}} while the result is of order \eqn{\nu^{-4}}.
#'
#' @details
#' Both parts are Laplace transforms:
#' \eqn{I_{\nu\nu} = \int_0^\infty e^{-\nu t}\phi(t)\,dt} with
#' \deqn{\phi(t) = t\,g(t) - \alpha + e^{-pt} - \gamma\,e^{-(p+2)t}, \qquad
#'   \alpha = \frac{p+4}{2(p+2)}, \quad \gamma = \frac{p}{2(p+2)},}
#' \eqn{g} as in [mvt_nu_coef()]. The first three derivatives of \eqn{\phi} at
#' zero vanish, so \eqn{I_{\nu\nu} = \sum_{k \ge 3} \phi^{(k)}(0)/\nu^{k+1}};
#' this series is summed for \eqn{\nu \ge \max(20, 4(p+2))} and the direct form
#' is used below.
#'
#' @param nu The degrees of freedom, a numeric vector of positive values.
#' @param p The dimension, a single positive whole number.
#'
#' @return A numeric vector as long as `nu`, the information, which is
#'   positive; the expected second derivative is its negative.
#'
#' @seealso [distrib_expected_hessian.MvStudentTDistrib()], the consumer, and
#'   [mvt_nu_coef()] for the series.
#'
#' @examples
#' nu <- c(10, 1e4, 1e8)
#' direct <- (trigamma(nu / 2) - trigamma((nu + 3) / 2)) / 4 -
#'   3 * (nu + 7) / (2 * nu * (nu + 3) * (nu + 5))
#' rbind(exact = distributions7:::mvt_info_nunu(nu, 3), direct = direct)
#'
#' @keywords internal
mvt_info_nunu <- function(nu, p) {
  out <- numeric(length(nu))
  big <- nu >= max(20, 4 * (p + 2))
  if (any(big)) {
    g <- mvt_nu_coef(p)
    K <- length(g) - 1L
    ga <- p / (2 * (p + 2))
    k <- 3:K
    # phi^(k)(0) = k g_(k-1) + (-p)^k - gamma (-(p+2))^k
    M <- factorial(k) * g[k] + (-p)^k - ga * (-(p + 2))^k
    x <- 1 / nu[big]
    s <- numeric(length(x))
    for (i in rev(seq_along(k))) s <- (s + M[i]) * x
    out[big] <- s * x^3
  }
  if (any(!big)) {
    v <- nu[!big]
    out[!big] <- (trigamma(v / 2) - trigamma((v + p) / 2)) / 4 -
      p * (v + p + 4) / (2 * v * (v + p) * (v + p + 2))
  }
  out
}


#' @title The Weights of a Multivariate t, Scaled for a Distant Observation
#'
#' @description
#' Computes \eqn{w/s}, \eqn{c = (\nu+p)/s} and \eqn{1/s}, with
#' \eqn{w = \Sigma^{-1}(y-\mu)} and \eqn{s = \nu + q}, without forming
#' \eqn{q}, which overflows for \eqn{\lvert y\rvert} above about
#' \eqn{10^{154}}. Each residual is divided by \eqn{m = \max(1,
#' \max_j\lvert r_j\rvert)} first, and the powers of \eqn{m} are carried
#' separately.
#'
#' @param y An \eqn{n \times p} numeric matrix.
#' @param pc The result of [mvt_pieces()].
#'
#' @return A named list with `ws` (an \eqn{n \times p} matrix), `cw` and `is`
#'   (numeric vectors of length \eqn{n}).
#'
#' @seealso [distrib_hess_y.MvStudentTDistrib()] and [mvt_y_tensor()], the
#'   consumers, and [mvt_weights()] for the unscaled quantities.
#'
#' @examples
#' d <- mvstudent_t1_distrib(1)
#' pc <- distributions7:::mvt_pieces(d, list(mu1 = 0, sigma_log_L1 = 0, nu = 5))
#' distributions7:::mvt_yscaled(matrix(c(1, 1e200)), pc)$cw
#'
#' @keywords internal
mvt_yscaled <- function(y, pc) {
  r <- sweep(y, 2L, pc$mu, "-")
  m <- if (nrow(r)) pmax(1, apply(abs(r), 1L, max)) else numeric(0)
  rt <- r / m
  wt <- rt %*% pc$sigma_inv
  qt <- mv_quad(rt, pc$s, pc$eta, pc$inverted)
  im2 <- (1 / m)^2
  st <- pc$nu * im2 + qt
  is <- im2 / st
  list(ws = wt / (m * st), cw = (pc$nu + pc$p) * is, is = is)
}


#' @title A Multivariate t Proposal for Integrating the Density
#' @name mv_reference_draw.MvStudentTDistrib
#'
#' @description
#' Draws from a multivariate t with the family's location, twice its scale
#' matrix and half its degrees of freedom, and returns the log-density of the
#' draws under that law. Its tails are heavier than the family's at every
#' \eqn{\nu}, so the importance ratio is bounded and its second moment finite;
#' the default gaussian proposal needs a finite covariance, which a t does not
#' have at \eqn{\nu \le 2}.
#'
#' @param distrib An [MvStudentTDistrib] object.
#' @param theta A named list of parameters, already aligned by the generic.
#' @param n The number of draws.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A named list with `y`, an \eqn{n \times p} numeric matrix, and
#'   `logd`, a numeric vector of length \eqn{n}.
#'
#' @seealso [mv_reference_draw()] for the generic and [check_distrib()] for the
#'   consumer.
#'
#' @examples
#' d <- mvstudent_t1_distrib(2)
#' theta <- list(mu1 = 0, mu2 = 0, sigma_log_L1 = 0, sigma_log_L2 = 0,
#'               sigma_L2.1 = 0.5, nu = 1.5)
#' set.seed(1)
#' r <- mv_reference_draw(d, theta, 20000)
#' mean(exp(distrib_pdf(d, r$y, theta, log = TRUE) - r$logd))
#'
#' @keywords internal
S7::method(mv_reference_draw, MvStudentTDistrib) <- function(distrib, theta, n, ...) {
  pc <- mvt_pieces(distrib, theta)
  p <- pc$p
  nu <- pc$nu / 2
  l <- sqrt(2) * mv_sampling_factor(pc$s, pc$eta, distrib@inverted)
  z <- matrix(stats::rnorm(n * p), n, p) %*% t(l)
  y <- sweep(z * sqrt(nu / stats::rchisq(n, df = nu)), 2L, pc$mu, "+")
  r <- sweep(y, 2L, pc$mu, "-")
  q <- mv_quad(r, pc$s, pc$eta, pc$inverted) / 2
  logd <- mvt_cp_d(nu, p, 0L) - (p / 2) * log(pi) -
    0.5 * (pc$logdet + p * log(2)) - ((nu + p) / 2) * log1p(q / nu)
  list(y = y, logd = logd)
}

#' @title Multivariate Student t Observed Hessian
#' @name distrib_hessian.MvStudentTDistrib
#'
#' @description
#' Computes the second derivatives of the log-density in closed form, by
#' differentiating the score of [distrib_gradient.MvStudentTDistrib()] once
#' more. Every block has a term in a derivative of the weight \eqn{c =
#' (\nu+p)/(\nu+q)}, which depends on the observation through \eqn{q}; for the
#' gaussian \eqn{c} is one and those terms are absent. Writing \eqn{w =
#' \Sigma^{-1}(y-\mu)} and \eqn{d = \nu + q}, the location block is
#' \deqn{\ell^{(\mu_a\mu_b)} = \frac{2c\,w_a w_b}{d} - c\,(\Sigma^{-1})_{ab},}
#' and \eqn{\partial c/\partial\nu = (q-p)/d^2} carries \eqn{\nu} into
#' every mixed block.
#'
#' @details
#' The pure \eqn{\nu} component is \eqn{c_p''(\nu) - h''(\nu)}, with
#' [mvt_cp_d()] for the part of \eqn{\nu} alone and [mvt_h_d()] for the part
#' that carries the data, each written without the cancellation of its direct
#' form; the result stays finite and accurate as \eqn{\nu} grows. The blocks in
#' the matrix read the derivatives \eqn{P_k} and \eqn{P_{kl}} of
#' \eqn{\Sigma^{-1}} from [mv_matrix_pieces()]: with \eqn{q_k = r^\top P_k r},
#' the matrix block is \eqn{-\tfrac{1}{2}\partial_{kl}\log\lvert\Sigma\rvert +
#' c\,q_kq_l/(2d) - \tfrac{c}{2}\, r^\top P_{kl} r}.
#'
#' @param distrib An [MvStudentTDistrib] object, from
#'   [mvstudent_t1_distrib()].
#' @param y An \eqn{n \times p} numeric matrix of observations. A vector of
#'   length \eqn{p} is read as a single observation.
#' @param theta A named list of parameters, each component a single number.
#' @param scale One of `"parameter"` (the default) or `"link"`. The method
#'   returns the parameter scale and the generic applies the link scale to its
#'   result. The two differ here, `nu` carrying a log link by default.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A named list of numeric vectors of length \eqn{n}, one per unordered
#'   pair of parameters, keyed as
#'   [`hess_names(distrib@params)`][hess_names]. Unlike the gaussian's, no
#'   block here is constant across rows.
#'
#' @section Notation:
#' \eqn{\mu} is the location, \eqn{\Sigma} the scale matrix, \eqn{\nu} the
#' degrees of freedom, \eqn{\eta} the free vector of the matrix
#' parametrization, \eqn{P_k} and \eqn{P_{kl}} the first and second derivatives
#' of \eqn{\Sigma^{-1}}, \eqn{r = y - \mu}, \eqn{q} the squared Mahalanobis
#' distance, \eqn{c} the weight and \eqn{\ell^{(ij)}} the second derivative of
#' the log-density in parameters \eqn{i} and \eqn{j}.
#'
#' @seealso [distrib_expected_hessian.MvStudentTDistrib()] for the expectation
#'   of this matrix, which is closed form,
#'   [distrib_gradient.MvStudentTDistrib()] for the score, [mvt_cp_d()] and
#'   [mvt_h_d()] for the \eqn{\nu} component, and [distrib_hessian()] for the
#'   generic.
#'
#' @examples
#' d <- mvstudent_t1_distrib(2)
#' theta <- list(mu1 = 0.5, mu2 = -0.3, sigma_log_L1 = 0.1,
#'               sigma_log_L2 = -0.2, sigma_L2.1 = 0.4, nu = 6)
#' set.seed(1)
#' y <- distrib_rng(d, 25, theta)
#'
#' H <- distrib_hessian(d, y, theta)
#' vapply(H, sum, numeric(1))
#'
#' # Against a numerical Hessian of the log-likelihood.
#' ll <- function(v) {
#'   t2 <- as.list(v); names(t2) <- d@params
#'   sum(distrib_pdf(d, y, t2, log = TRUE))
#' }
#' Hn <- numDeriv::hessian(ll, unlist(theta))
#' ref <- vapply(distributions7:::hess_pairs(d@params),
#'               function(q) Hn[match(q[1], d@params), match(q[2], d@params)],
#'               numeric(1))
#' max(abs(vapply(H, sum, numeric(1)) - ref))
#'
#' # The location block moves with the observation, where a gaussian's does
#' # not: the curvature at a far row is smaller.
#' z <- distributions7:::mvt_weights(y, distributions7:::mvt_pieces(d, theta))
#' near <- which.min(z$q); far <- which.max(z$q)
#' c(q_near = z$q[near], curv_near = H$mu1_mu1[near],
#'   q_far = z$q[far], curv_far = H$mu1_mu1[far])
#'
#' @keywords internal
S7::method(distrib_hessian, MvStudentTDistrib) <- function(distrib, y, theta,
                                                           scale = c("parameter", "link"),
                                                           ...) {
  y <- as_mv_matrix(distrib, y)
  pc <- mvt_pieces(distrib, theta, derivs2 = TRUE)
  z <- mvt_weights(y, pc)
  n <- nrow(y)
  p <- pc$p
  nu <- pc$nu
  si <- pc$sigma_inv
  cw <- z$cw
  den <- nu + z$q

  d2ld <- pc$d2logdet
  spair <- param_pair_lookup(pc$s)
  # P_k r and q_k = r' P_k r, with P_k = d_k(Sigma^-1) read from the chart:
  # the sandwich Sigma^-1 A_k Sigma^-1 cancels near a singular scale matrix.
  # q_k is the derivative of q in the matrix directions and appears in three
  # of the blocks.
  rp <- lapply(pc$ai, function(ak) z$r %*% ak)
  qk <- lapply(rp, function(m) rowSums(m * z$r))
  # dc/dnu, which is where the degrees of freedom enter every mixed block.
  dc_dnu <- (z$q - p) / den^2

  n_par <- distrib@n_params
  pairs <- mv_hess_indices(distrib)
  out <- vector("list", length(pairs))
  names(out) <- hess_names(distrib@params)

  for (i in seq_along(pairs)) {
    a <- pairs[[i]][1L]
    b <- pairs[[i]][2L]
    out[[i]] <- if (a <= p && b <= p) {
      # d/dmu_b [c w_a]: dc/dmu_b = 2 c w_b / (nu + q), dw_a/dmu_b = -Sigma^-1
      2 * cw * z$w[, a] * z$w[, b] / den - cw * si[a, b]
    } else if (a <= p && b < n_par) {
      k <- b - p
      # d/deta_k [c w_a]: dc = -c q_k / (nu + q), dw_a = (P_k r)_a
      -cw * qk[[k]] * z$w[, a] / den + cw * rp[[k]][, a]
    } else if (a <= p && b == n_par) {
      z$w[, a] * dc_dnu
    } else if (a < n_par && b < n_par) {
      k <- a - p
      l <- b - p
      idx <- spair[[paste(min(k, l), max(k, l), sep = ":")]]
      -0.5 * d2ld[[idx]] +
        cw * qk[[l]] * qk[[k]] / (2 * den) -
        0.5 * cw * rowSums((z$r %*% pc$ai2[[idx]]) * z$r)
    } else if (a < n_par && b == n_par) {
      -0.5 * qk[[a - p]] * dc_dnu
    } else {
      # d2l/dnu2 = D^2 c_p - D^2 h, each written without its cancellation:
      # see mvt_cp_d() and mvt_h_d()
      mvt_cp_d(nu, p, 2L) - mvt_h_d(z$q / nu, nu, p, 2L)
    }
    if (length(out[[i]]) == 1L) out[[i]] <- rep(out[[i]], n)
  }
  out
}


#' @title Multivariate Student t Expected Information
#' @name distrib_expected_hessian.MvStudentTDistrib
#'
#' @description
#' Computes the expectation of the observed Hessian in closed form, taken from
#' the family's own scale mixture. Every component is a Beta moment or a
#' polygamma, so nothing is sampled and no quadrature runs, and two calls at the
#' same parameters return the same numbers. The components do not depend on the
#' data, so every returned vector is constant across rows and `y` is used only
#' for its row count.
#'
#' @details
#' # Why the expectations separate
#'
#' Writing \eqn{q = z^\top\Sigma^{-1}z} and \eqn{c = (\nu+p)/(\nu+q)}, the
#' variable \eqn{v = q/(q+\nu)} has exactly a \eqn{\mathrm{Beta}(p/2, \nu/2)}
#' distribution and is independent of the direction \eqn{e =
#' \Sigma^{-1/2}z/\lVert\Sigma^{-1/2}z\rVert} of the whitened vector, which is
#' uniform on the sphere. Each expectation therefore factors into a radial part
#' and an
#' angular part: the radial one gives \eqn{cq = (\nu+p)v} and
#' \eqn{c^2q^2 = (\nu+p)^2v^2}, and the angular one gives
#' \eqn{\mathbb{E}[ee^\top] = I/p} and
#' \eqn{\mathbb{E}[(e^\top Be)(e^\top Ce)] =
#' \{\operatorname{tr}B\operatorname{tr}C + 2\operatorname{tr}(BC)\}/
#' \{p(p+2)\}}.
#'
#' # The blocks
#'
#' With \eqn{t_k = \operatorname{tr}(\Sigma^{-1}A_k)} and
#' \eqn{T_{kl} = \operatorname{tr}(\Sigma^{-1}A_k\Sigma^{-1}A_l)},
#' \deqn{I_{\mu\mu} = \frac{\nu+p}{\nu+p+2}\,\Sigma^{-1}, \qquad
#'   I_{kl} = \frac{(\nu+p)T_{kl} - t_kt_l}{2(\nu+p+2)}, \qquad
#'   I_{k\nu} = \frac{-t_k}{(\nu+p)(\nu+p+2)},}
#' \deqn{I_{\nu\nu} = \frac{\psi'(\nu/2) - \psi'\{(\nu+p)/2\}}{4}
#'   - \frac{p}{\nu(\nu+p)} + \frac{p}{2\nu(\nu+p+2)}.}
#' The location is orthogonal to everything else, every cross-expectation with
#' it being odd in \eqn{z}.
#'
#' # The gaussian limit
#'
#' Each block reduces to the gaussian's as \eqn{\nu \to \infty}: \eqn{I_{\mu\mu}
#' \to \Sigma^{-1}}, \eqn{I_{kl} \to T_{kl}/2}, and \eqn{I_{\nu\nu} \to 0}, the
#' degrees of freedom ceasing to be identified in the limit; \eqn{I_{\nu\nu}} is
#' of order \eqn{\nu^{-4}}, so the standard error of \eqn{\nu} is large unless
#' the tail is heavy. Its terms as written above are of order \eqn{\nu^{-2}} and
#' cancel, so the block is computed by [mvt_info_nunu()].
#'
#' @param distrib An [MvStudentTDistrib] object, from
#'   [mvstudent_t1_distrib()].
#' @param y An \eqn{n \times p} numeric matrix of observations. Only its row
#'   count is used: the expectation is taken over the law, so no observation
#'   enters any component.
#' @param theta A named list of parameters, each component a single number.
#' @param scale One of `"parameter"` (the default) or `"link"`. The method
#'   returns the parameter scale and the generic applies the link scale to its
#'   result. The two differ here, `nu` carrying a log link by default.
#' @param approx Ignored: the expectation is exact and no approximation
#'   strategy is consulted. [fit_distrib()] rejects a `fisher_scoring(approx =)`
#'   for this family for that reason.
#' @param nsim Ignored, for the same reason.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A named list of numeric vectors of length \eqn{n}, keyed as
#'   [`hess_names(distrib@params)`][hess_names], each vector constant. The
#'   returned quantities are the expected second derivatives, so the
#'   information is their negative.
#'
#' @section Notation:
#' \eqn{\mu} is the location, \eqn{\Sigma} the scale matrix, \eqn{\nu} the
#' degrees of freedom, \eqn{\eta} the free vector of the matrix
#' parametrization, \eqn{A_k = \partial\Sigma/\partial\eta_k}, \eqn{\psi'} the
#' trigamma function, \eqn{e} a uniform direction on the unit sphere and
#' \eqn{I} an expected information block.
#'
#' @seealso [distrib_hessian.MvStudentTDistrib()] for the observed matrix,
#'   [distrib_expected_hessian.MvGaussianDistrib()] for the limiting family,
#'   [fit_distrib()], whose Fisher scoring inverts this matrix, and
#'   [distrib_expected_hessian()] for the generic.
#'
#' @examples
#' d <- mvstudent_t1_distrib(2)
#' theta <- list(mu1 = 0.5, mu2 = -0.3, sigma_log_L1 = 0.1,
#'               sigma_log_L2 = -0.2, sigma_L2.1 = 0.4, nu = 6)
#'
#' EH <- distrib_expected_hessian(d, matrix(0, 3, 2), theta)
#' round(vapply(EH, function(z) z[1], numeric(1)), 5)
#'
#' # Closed form, so two calls agree to the bit and nothing is sampled.
#' identical(EH, distrib_expected_hessian(d, matrix(0, 3, 2), theta))
#'
#' # Every location-matrix and location-nu component is exactly zero.
#' orth <- grep("^mu[0-9]+_(sigma|nu)", names(EH), value = TRUE)
#' vapply(EH[orth], function(z) z[1], numeric(1))
#'
#' # And the closed form is what averaging the observed Hessian converges to.
#' set.seed(3)
#' big <- distrib_rng(d, 50000, theta)
#' round(rbind(sampled = vapply(distrib_hessian(d, big, theta), mean, numeric(1)),
#'             closed = vapply(EH, function(z) z[1], numeric(1))), 3)
#'
#' # The nu block vanishes as the tail lightens.
#' vapply(c(6, 30, 200, 5000), function(nu) {
#'   t2 <- theta; t2$nu <- nu
#'   -distrib_expected_hessian(d, matrix(0, 1, 2), t2)$nu_nu[1]
#' }, numeric(1))
#'
#' @keywords internal
S7::method(distrib_expected_hessian, MvStudentTDistrib) <- function(
    distrib, y, theta, scale = c("parameter", "link"),
    approx = c("opg", "bartlett", "integrate", "mc"), nsim = 10000, ...) {
  y <- as_mv_matrix(distrib, y)
  pc <- mvt_pieces(distrib, theta, derivs = TRUE)
  n <- nrow(y)
  p <- pc$p
  nu <- pc$nu
  si <- pc$sigma_inv
  sa <- lapply(pc$a, function(ak) si %*% ak)
  tk <- vapply(sa, function(m) sum(diag(m)), numeric(1))
  npr <- nu + p

  # (psi'(nu/2) - psi'((nu + p)/2))/4 - p (nu + p + 4)/(2 nu (nu + p)(nu + p + 2)),
  # whose terms cancel to order nu^-4: see mvt_info_nunu()
  i_nunu <- mvt_info_nunu(nu, p)

  pairs <- mv_hess_indices(distrib)
  nq <- distrib@n_params
  out <- vector("list", length(pairs))
  names(out) <- hess_names(distrib@params)

  for (i in seq_along(pairs)) {
    a <- pairs[[i]][1L]
    b <- pairs[[i]][2L]
    a_mu <- a <= p
    b_mu <- b <= p
    a_nu <- a == nq
    b_nu <- b == nq
    v <- if (a_mu && b_mu) {
      -npr / (npr + 2) * si[a, b]
    } else if (a_mu || b_mu) {
      0
    } else if (a_nu && b_nu) {
      -i_nunu
    } else if (a_nu || b_nu) {
      k <- if (a_nu) b - p else a - p
      tk[[k]] / (npr * (npr + 2))
    } else {
      ka <- a - p
      kb <- b - p
      -(npr * sum(sa[[ka]] * t(sa[[kb]])) - tk[[ka]] * tk[[kb]]) /
        (2 * (npr + 2))
    }
    out[[i]] <- rep(v, n)
  }
  out
}


#' @title Multivariate Student t Response Gradient
#' @name distrib_grad_y.MvStudentTDistrib
#'
#' @description
#' Computes the derivative of the log-density in the response,
#' \deqn{\frac{\partial \ell}{\partial y} = -c\,\Sigma^{-1}(y-\mu), \qquad c =
#' \frac{\nu+p}{\nu+q},} one row per observation: the gaussian's expression with
#' the family's weight in front of it. A gaussian's response gradient grows
#' without limit as the observation moves away; this one rises, turns and decays
#' like \eqn{1/\lVert y\rVert}, the redescending score of a heavy-tailed family.
#' It is minus the score in the location, this being a location family. It is
#' computed as \eqn{-(\nu+p)\,w/s} with \eqn{s = \nu + q} (see [mvt_yscaled()]),
#' which stays finite for a distant observation.
#'
#' @param distrib An [MvStudentTDistrib] object, from
#'   [mvstudent_t1_distrib()].
#' @param y An \eqn{n \times p} numeric matrix of observations. A vector of
#'   length \eqn{p} is read as a single observation.
#' @param theta A named list of parameters, each component a single number.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return An \eqn{n \times p} numeric matrix, row \eqn{i} holding
#'   \eqn{\partial\ell_i/\partial y_i}.
#'
#' @section Notation:
#' \eqn{\mu} is the location, \eqn{\Sigma} the scale matrix, \eqn{\nu} the
#' degrees of freedom, \eqn{p} the dimension, \eqn{q} the squared Mahalanobis
#' distance and \eqn{c} the weight.
#'
#' @seealso [distrib_hess_y.MvStudentTDistrib()] for the second derivative,
#'   [distrib_cross_y.MvStudentTDistrib()] for the mixed block,
#'   [mvt_weights()] for the weight, and [distrib_grad_y()] for the generic.
#'
#' @examples
#' d <- mvstudent_t1_distrib(2)
#' theta <- list(mu1 = 0.5, mu2 = -0.3, sigma_log_L1 = 0.1,
#'               sigma_log_L2 = -0.2, sigma_L2.1 = 0.4, nu = 6)
#' set.seed(1)
#' y <- distrib_rng(d, 4, theta)
#'
#' distrib_grad_y(d, y, theta)
#'
#' # Against a numerical derivative taken row by row.
#' num <- t(apply(y, 1, function(r)
#'   numDeriv::grad(function(z) distrib_pdf(d, z, theta, log = TRUE), r)))
#' max(abs(distrib_grad_y(d, y, theta) - num))
#'
#' # It is minus the score in the location.
#' g <- distrib_gradient(d, y, theta)
#' all.equal(distrib_grad_y(d, y, theta), -cbind(g$mu1, g$mu2),
#'           check.attributes = FALSE)
#'
#' # It redescends: the gaussian's grows with the distance and the t's does
#' # not, along the first coordinate from the location.
#' far <- cbind(0.5 + c(1, 3, 10, 40), -0.3)
#' g0 <- mvgaussian1_distrib(2)
#' rbind(t = distrib_grad_y(d, far, theta)[, 1],
#'       gaussian = distrib_grad_y(g0, far, theta[1:5])[, 1])
#'
#' @keywords internal
S7::method(distrib_grad_y, MvStudentTDistrib) <- function(distrib, y, theta, ...) {
  y <- as_mv_matrix(distrib, y)
  pc <- mvt_pieces(distrib, theta)
  z <- mvt_yscaled(y, pc)
  -(pc$nu + pc$p) * z$ws
}


#' @title The Pieces Every Mixed Response Derivative of the t Is Written In
#'
#' @description
#' Assembles, once per call, the first and second derivatives in the parameters
#' of the quantities from which every mixed response derivative of this family
#' is built: \eqn{s = \nu + q}, the residual mapped through the inverse scale
#' matrix \eqn{w = \Sigma^{-1}(y-\mu)}, the inverse scale matrix
#' \eqn{\Sigma^{-1}}, and the two scalars \eqn{c = (\nu+p)/s} and \eqn{d =
#' (\nu+p)/s^2}. The response derivatives are \eqn{\ell^{(y)} = -c\,w} and
#' \eqn{\ell^{(yy)} = -c\,\Sigma^{-1} + 2d\,ww^\top}, so differentiating either
#' in the parameters is a product rule over these pieces.
#'
#' @details
#' # First derivatives
#'
#' Writing \eqn{A_k} for the matrix parametrization's first derivatives,
#' \eqn{u_j} for the \eqn{j}th column of \eqn{\Sigma^{-1}} and \eqn{P_k} for
#' \eqn{\Sigma^{-1}A_k\Sigma^{-1}},
#' \deqn{\partial_{\mu_j} s = -2w_j, \quad \partial_{\eta_k} s = -w^\top A_kw,
#'   \quad \partial_\nu s = 1,}
#' \deqn{\partial_{\mu_j} w = -u_j, \quad
#'   \partial_{\eta_k} w = -\Sigma^{-1}A_kw, \quad \partial_\nu w = 0,
#'   \qquad \partial_{\eta_k}\Sigma^{-1} = -P_k.}
#'
#' # Second derivatives
#'
#' All vanish except
#' \deqn{\partial_{\mu_j\mu_l} s = 2(\Sigma^{-1})_{jl}, \quad
#'   \partial_{\mu_j\eta_k} s = 2(\Sigma^{-1}A_kw)_j, \quad
#'   \partial_{\mu_j\eta_k} w = P_ke_j,}
#' and, with \eqn{B_{kl} = A_k\Sigma^{-1}A_l + A_l\Sigma^{-1}A_k - A_{kl}},
#' \deqn{\partial_{\eta_k\eta_l} s = w^\top B_{kl}w, \quad
#'   \partial_{\eta_k\eta_l} w = \Sigma^{-1}B_{kl}w, \quad
#'   \partial_{\eta_k\eta_l}\Sigma^{-1} = \Sigma^{-1}B_{kl}\Sigma^{-1}.}
#' The arrays \eqn{P_k} and \eqn{\Sigma^{-1}B_{kl}\Sigma^{-1}} are the first and
#' second derivatives of \eqn{\Sigma^{-1}}, read from [mv_matrix_pieces()] and
#' not formed from the derivatives of \eqn{\Sigma}; the derivatives in \eqn{s}
#' and \eqn{w} are quadratic and linear forms in them. \eqn{c} and \eqn{d}
#' follow by the quotient rule, \eqn{\nu+p} being linear in \eqn{\nu} and
#' constant in everything else; the derivative of \eqn{c} in \eqn{\nu} is
#' written as \eqn{(q-p)/s^2} and its second derivative as \eqn{2(p-q)/s^3}, the
#' two terms of the quotient rule agreeing to leading order as \eqn{\nu} grows.
#'
#' @param distrib An [MvStudentTDistrib] object, from
#'   [mvstudent_t1_distrib()].
#' @param y An \eqn{n \times p} numeric matrix of observations, already coerced
#'   by [as_mv_matrix()].
#' @param theta A named list of parameters, each component a single number.
#'
#' @return A named list of the quantities above, indexed by parameter position:
#'   the \eqn{p} locations, then the matrix parametrization's free values, then
#'   \eqn{\nu}. Its components are `p`, `n` and `npar` (the dimension, the
#'   number of rows and the number of parameters), `w`, `si`, `s`, `cw` and `dw`
#'   (the values), `sa`, `wa`, `Sa`, `ca` and `da` (lists of first derivatives,
#'   one per parameter), and `pair(a, b)`, a function returning the second-order
#'   pieces `s`, `w`, `S`, `c` and `d` for one pair of positions.
#'
#' @section Notation:
#' \eqn{\mu} is the location, \eqn{\Sigma} the scale matrix, \eqn{\nu} the
#' degrees of freedom, \eqn{\eta} the free vector of the matrix
#' parametrization, \eqn{A_k} and \eqn{A_{kl}} its derivative arrays, \eqn{q}
#' the squared Mahalanobis distance, \eqn{w = \Sigma^{-1}(y-\mu)} and \eqn{e_j}
#' the \eqn{j}th standard basis vector.
#'
#' @seealso [distrib_cross_y.MvStudentTDistrib()],
#'   [distrib_cross2_y.MvStudentTDistrib()],
#'   [distrib_grad_y_hess.MvStudentTDistrib()] and
#'   [distrib_hess_y_hess.MvStudentTDistrib()], the four consumers, and
#'   [mvt_weights()] for the values these differentiate.
#'
#' @examples
#' d <- mvstudent_t1_distrib(2)
#' theta <- list(mu1 = 0.5, mu2 = -0.3, sigma_log_L1 = 0.1,
#'               sigma_log_L2 = -0.2, sigma_L2.1 = 0.4, nu = 6)
#' set.seed(1)
#' y <- distrib_rng(d, 4, theta)
#' z <- distributions7:::mvt_dpieces(d, y, theta)
#' names(z)
#'
#' # Six parameters, so six first-derivative entries per quantity.
#' c(npar = z$npar, n_sa = length(z$sa), n_wa = length(z$wa))
#'
#' # ds/dnu is exactly one, s being nu + q.
#' z$sa[[6]]
#'
#' # And ds/dmu1 is -2 w1, against a difference of q + nu.
#' h <- 1e-5
#' sfun <- function(m1) {
#'   t2 <- theta; t2$mu1 <- m1
#'   pc <- distributions7:::mvt_pieces(d, t2)
#'   distributions7:::mvt_weights(y, pc)$q + t2$nu
#' }
#' max(abs(z$sa[[1]] - (sfun(0.5 + h) - sfun(0.5 - h)) / (2 * h)))
#'
#' @keywords internal
mvt_dpieces <- function(distrib, y, theta) {
  p <- distrib@n_dim
  pc <- mvt_pieces(distrib, theta, derivs2 = TRUE)
  z <- mvt_weights(y, pc)
  si <- pc$sigma_inv
  w <- z$w
  n <- nrow(w)
  nf <- pc$s@n_free
  npar <- p + nf + 1L
  nu <- pc$nu
  N <- nu + p
  s <- nu + z$q
  zero_m <- matrix(0, p, p)

  # Sigma^-1 A_k w, as rows, and w' A_k w -- the two quantities every matrix
  # component is written in
  # Sigma^-1 A_k Sigma^-1 = -d_k(Sigma^-1), from the chart and not as a
  # sandwich, which cancels near a singular scale matrix
  v <- lapply(pc$ai, function(ak) -(z$r %*% ak))
  gq <- lapply(pc$ai, function(ak) -rowSums((z$r %*% ak) * z$r))
  P <- lapply(pc$ai, function(ak) -ak)

  sa <- wa <- Sa <- vector("list", npar)
  Na <- numeric(npar)
  for (j in seq_len(p)) {
    sa[[j]] <- -2 * w[, j]
    wa[[j]] <- matrix(rep(-si[j, ], each = n), n, p)
    Sa[[j]] <- zero_m
  }
  for (k in seq_len(nf)) {
    sa[[p + k]] <- -gq[[k]]
    wa[[p + k]] <- -v[[k]]
    Sa[[p + k]] <- -P[[k]]
  }
  sa[[npar]] <- rep(1, n)
  wa[[npar]] <- matrix(0, n, p)
  Sa[[npar]] <- zero_m
  Na[npar] <- 1

  # d/dnu of c = N/s is 1/s - N/s^2, written as (q - p)/s^2: the two terms
  # agree to leading order as nu grows
  ca <- lapply(seq_len(npar), function(a) {
    if (a == npar) (z$q - p) / s^2 else -N * sa[[a]] / s^2
  })
  da <- lapply(seq_len(npar), function(a) Na[a] / s^2 - 2 * N * sa[[a]] / s^3)

  pair <- function(a, b) {
    if (a > b) {
      t <- a
      a <- b
      b <- t
    }
    s_ab <- rep(0, n)
    w_ab <- matrix(0, n, p)
    S_ab <- zero_m
    if (a <= p && b <= p) {
      s_ab <- rep(2 * si[a, b], n)
    } else if (a <= p && b <= p + nf) {
      k <- b - p
      s_ab <- 2 * v[[k]][, a]
      w_ab <- matrix(rep(P[[k]][a, ], each = n), n, p)
    } else if (a > p && a <= p + nf && b <= p + nf) {
      ka <- a - p
      kb <- b - p
      # the second derivative of the inverse, from the chart: s_ab = r' P_kl r
      # and w_ab = P_kl r, neither formed as a sandwich of Sigma^-1
      S_ab <- mv_ai2(pc, ka, kb)
      w_ab <- z$r %*% S_ab
      s_ab <- rowSums(w_ab * z$r)
    }
    # nu is linear in the recursion: every pair naming it has vanishing
    # second derivatives of s, w and the inverse, and only c and d move
    c_ab <- if (a == npar && b == npar) {
      # -2/s^2 + 2N/s^3, whose terms agree to leading order as nu grows
      2 * (p - z$q) / s^3
    } else {
      -Na[a] * sa[[b]] / s^2 - Na[b] * sa[[a]] / s^2 -
        N * s_ab / s^2 + 2 * N * sa[[a]] * sa[[b]] / s^3
    }
    d_ab <- -2 * Na[a] * sa[[b]] / s^3 - 2 * Na[b] * sa[[a]] / s^3 -
      2 * N * s_ab / s^3 + 6 * N * sa[[a]] * sa[[b]] / s^4
    list(s = s_ab, w = w_ab, S = S_ab, c = c_ab, d = d_ab)
  }

  list(p = p, n = n, npar = npar, si = si, w = w, s = s,
       cw = N / s, dw = N / s^2, sa = sa, wa = wa, Sa = Sa,
       ca = ca, da = da, pair = pair)
}

# Commented out on 2026-10-12 (T4): the mixed response derivatives read the chart's second derivatives of the
# inverse, so this has no caller.
# Kept rather than deleted.
# The Second Derivative of the Scale Matrix, by Position
# @description As [mvg_a2()], for the Student t's pieces.
# @param pc The result of [mvt_pieces()] with `derivs2`.
# @param k,l Positions among the structure's free values.
# @return A \eqn{p \times p} numeric matrix.
# @keywords internal
# mvt_a2 <- function(pc, k, l) {
#   nm <- pc$s@free_names
#   ij <- sort(c(k, l))
#   pc$a2[[paste(nm[ij], collapse = ":")]]
# }

#' @title An Outer Product Per Observation
#'
#' @description
#' Builds the \eqn{p \times p \times n} array whose \eqn{i}th slice is
#' \eqn{a_ib_i^\top}, the outer product of row \eqn{i} of `a` with row \eqn{i}
#' of `b`. The whole array is assembled by two column replications and one
#' permutation, with no loop over the observations.
#'
#' @param a,b Numeric matrices with the same dimensions, \eqn{n} rows and
#'   \eqn{p} columns. Not validated: matrices of different dimensions signal the
#'   error of the arithmetic underneath.
#'
#' @return A \eqn{p \times p \times n} numeric array. It is symmetric slice by
#'   slice only when `a` and `b` are the same matrix.
#'
#' @seealso [.mvt_scale_slices()] for the other array helper and
#'   [distrib_cross2_y.MvStudentTDistrib()] for the consumer.
#'
#' @examples
#' a <- matrix(1:6, 3, 2)
#' b <- matrix(c(1, 0, -1, 2, 1, 0), 3, 2)
#' o <- distributions7:::mv_outer_rows(a, b)
#' dim(o)
#'
#' # Slice i is the outer product of the two i-th rows.
#' o[, , 1]
#' outer(a[1, ], b[1, ])
#' all.equal(o[, , 2], outer(a[2, ], b[2, ]))
#'
#' # With one matrix given twice, every slice is symmetric.
#' s <- distributions7:::mv_outer_rows(a, a)
#' all(vapply(seq_len(dim(s)[3]),
#'            function(i) isSymmetric(s[, , i]), TRUE))
#'
#' @keywords internal
mv_outer_rows <- function(a, b) {
  n <- nrow(a)
  p <- ncol(a)
  aperm(array(a[, rep(seq_len(p), times = p), drop = FALSE] *
                b[, rep(seq_len(p), each = p), drop = FALSE], c(n, p, p)),
        c(2L, 3L, 1L))
}

#' @title Multivariate Student t Mixed Response-Parameter Derivatives
#' @name distrib_cross_y.MvStudentTDistrib
#'
#' @description
#' Computes \eqn{\partial^2\ell/\partial y\,\partial\theta_k}, one
#' \eqn{n \times p} matrix per parameter. The response gradient is
#' \eqn{-c\,w}, so every component carries the derivative of the weight beside
#' the gaussian term it multiplies. With
#' \eqn{A_k = \partial\Sigma/\partial\eta_k},
#' \deqn{\frac{\partial^2\ell}{\partial y\,\partial\mu_j}
#'     = c\,\Sigma^{-1}e_j - \frac{\partial c}{\partial\mu_j}\,w,
#'   \qquad \frac{\partial c}{\partial\mu_j}
#'     = \frac{2(\nu+p)\,w_j}{(\nu+q)^2},}
#' \deqn{\frac{\partial^2\ell}{\partial y\,\partial\eta_k}
#'     = c\,\Sigma^{-1}A_k w - \frac{\partial c}{\partial\eta_k}\,w,
#'   \qquad \frac{\partial c}{\partial\eta_k}
#'     = \frac{(\nu+p)\,w^\top A_k w}{(\nu+q)^2},}
#' \deqn{\frac{\partial^2\ell}{\partial y\,\partial\nu}
#'     = -\frac{(q-p)\,w}{(\nu+q)^2}.}
#' For the gaussian the location block, \eqn{\Sigma^{-1}e_j}, is the same for
#' every row; here every component varies with the row, the location block
#' carrying the observation through \eqn{c}.
#'
#' @details
#' The log-density contains no distribution function, only `lgamma`, a logarithm
#' and a quadratic form, so every derivative in \eqn{\nu} is closed form. As
#' \eqn{\nu \to \infty} the weight tends to one and its derivatives to zero, so
#' every component that the two families share tends to the gaussian's.
#'
#' @param distrib An [MvStudentTDistrib] object, from
#'   [mvstudent_t1_distrib()].
#' @param y An \eqn{n \times p} numeric matrix of observations. A vector of
#'   length \eqn{p} is read as a single observation.
#' @param theta A named list of parameters, each component a single number.
#' @param scale One of `"parameter"` (the default) or `"link"`. The method
#'   returns the parameter scale and the generic applies the link scale to its
#'   result. The two differ here, `nu` carrying a log link by default.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A named list of \eqn{n \times p} numeric matrices, one per
#'   parameter, in `distrib@params` order.
#'
#' @section Notation:
#' \eqn{\mu} is the location, \eqn{\Sigma} the scale matrix, \eqn{\nu} the
#' degrees of freedom, \eqn{\eta} the free vector of the matrix
#' parametrization, \eqn{q} the squared Mahalanobis distance, \eqn{c} the
#' weight, \eqn{w = \Sigma^{-1}(y-\mu)} and \eqn{e_j} the \eqn{j}th standard
#' basis vector.
#'
#' @seealso [distrib_grad_y.MvStudentTDistrib()], whose derivative in the
#'   parameters this is, [mvt_dpieces()] for the assembly,
#'   [distrib_cross_y.MvGaussianDistrib()] for the limiting family, and
#'   [distrib_cross_y()] for the generic.
#'
#' @examples
#' d <- mvstudent_t1_distrib(2)
#' theta <- list(mu1 = 0.5, mu2 = -0.3, sigma_log_L1 = 0.1,
#'               sigma_log_L2 = -0.2, sigma_L2.1 = 0.4, nu = 6)
#' set.seed(1)
#' y <- distrib_rng(d, 4, theta)
#'
#' cy <- distrib_cross_y(d, y, theta)
#' names(cy)
#'
#' # Against a difference of the response gradient, parameter by parameter.
#' h <- 1e-5
#' vapply(seq_along(d@params), function(k) {
#'   tp <- theta; tp[[k]] <- tp[[k]] + h
#'   tm <- theta; tm[[k]] <- tm[[k]] - h
#'   max(abs(cy[[k]] -
#'           (distrib_grad_y(d, y, tp) - distrib_grad_y(d, y, tm)) / (2 * h)))
#' }, numeric(1))
#'
#' # As nu grows every shared component approaches the gaussian's.
#' t2 <- theta; t2$nu <- 1e7
#' g <- mvgaussian1_distrib(2)
#' cg <- distrib_cross_y(g, y, theta[1:5])
#' ct <- distrib_cross_y(d, y, t2)
#' max(vapply(1:5, function(k) max(abs(cg[[k]] - ct[[k]])), numeric(1)))
#'
#' @keywords internal
S7::method(distrib_cross_y, MvStudentTDistrib) <-
  function(distrib, y, theta, scale = c("parameter", "link"), ...) {
    y <- as_mv_matrix(distrib, y)
    # one source for every order: l^(y) = -c w, so the first derivative is
    # -c_a w - c w_a and the pieces supply both
    z <- mvt_dpieces(distrib, y, theta)
    stats::setNames(lapply(seq_len(z$npar), function(a) {
      -z$ca[[a]] * z$w - z$cw * z$wa[[a]]
    }), distrib@params)
  }


#' @title Multivariate Student t Third Derivative in Two Responses and One Parameter
#' @name distrib_cross2_y.MvStudentTDistrib
#'
#' @description
#' Computes \eqn{\partial^3\ell/\partial y\,\partial y^\top \partial\theta_a},
#' one \eqn{p \times p \times n} array per parameter. Writing
#' \eqn{M = \ell^{(yy)} = -c\,\Sigma^{-1} + 2d\,ww^\top} with
#' \eqn{c = (\nu+p)/s}, \eqn{d = (\nu+p)/s^2} and \eqn{s = \nu+q},
#' \deqn{\partial_a M = -c_a\Sigma^{-1} - c\,\partial_a\Sigma^{-1}
#'   + 2d_a ww^\top + 2d\left(w_aw^\top + ww_a^\top\right),}
#' every piece coming from [mvt_dpieces()].
#'
#' @details
#' For the gaussian the response Hessian is \eqn{-\Sigma^{-1}}, one matrix for
#' the whole sample, so its derivative in a parameter is one \eqn{p \times p}
#' matrix. For the t the response Hessian depends on the observation through the
#' weight, so this method and [distrib_hess_y_hess.MvStudentTDistrib()] return
#' one matrix per row. The location components of the gaussian are exactly zero,
#' while those of the t are not zero in general.
#'
#' This is one of the three derivatives that a marginal criterion reads when the
#' family is a prior over a coefficient block.
#'
#' @param distrib An [MvStudentTDistrib] object, from
#'   [mvstudent_t1_distrib()].
#' @param y An \eqn{n \times p} numeric matrix of observations. A vector of
#'   length \eqn{p} is read as a single observation.
#' @param theta A named list of parameters, each component a single number.
#' @param scale One of `"parameter"` (the default) or `"link"`. The method
#'   returns the parameter scale and the generic applies the link scale to its
#'   result. The two differ here.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A named list of \eqn{p \times p \times n} numeric arrays, one per
#'   parameter, in `distrib@params` order.
#'
#' @section Notation:
#' \eqn{\Sigma} is the scale matrix, \eqn{\nu} the degrees of freedom, \eqn{q}
#' the squared Mahalanobis distance, \eqn{s = \nu+q}, \eqn{c} and \eqn{d} the
#' two weights, \eqn{w = \Sigma^{-1}(y-\mu)}, and a subscript \eqn{a} denotes
#' a derivative in parameter \eqn{a}.
#'
#' @seealso [distrib_hess_y.MvStudentTDistrib()], whose derivative in the
#'   parameters this is, [distrib_grad_y_hess.MvStudentTDistrib()] and
#'   [distrib_hess_y_hess.MvStudentTDistrib()] for the other two a marginal
#'   criterion reads, [mvt_dpieces()] for the assembly, and
#'   [distrib_cross2_y()] for the generic.
#'
#' @examples
#' d <- mvstudent_t1_distrib(2)
#' theta <- list(mu1 = 0.5, mu2 = -0.3, sigma_log_L1 = 0.1,
#'               sigma_log_L2 = -0.2, sigma_L2.1 = 0.4, nu = 6)
#' set.seed(1)
#' y <- distrib_rng(d, 4, theta)
#'
#' c2 <- distrib_cross2_y(d, y, theta)
#' dim(c2$sigma_L2.1)
#'
#' # One matrix per row, so the location component is not zero as it is for a
#' # gaussian.
#' c2$mu1[, , 1]
#' distrib_cross2_y(mvgaussian1_distrib(2), y, theta[1:5])$mu1
#'
#' # Against a difference of the response Hessian.
#' h <- 1e-5
#' vapply(seq_along(d@params), function(k) {
#'   tp <- theta; tp[[k]] <- tp[[k]] + h
#'   tm <- theta; tm[[k]] <- tm[[k]] - h
#'   max(abs(c2[[k]] -
#'           (distrib_hess_y(d, y, tp) - distrib_hess_y(d, y, tm)) / (2 * h)))
#' }, numeric(1))
#'
#' @keywords internal
S7::method(distrib_cross2_y, MvStudentTDistrib) <-
  function(distrib, y, theta, scale = c("parameter", "link"), ...) {
    y <- as_mv_matrix(distrib, y)
    z <- mvt_dpieces(distrib, y, theta)
    ww <- mv_outer_rows(z$w, z$w)
    stats::setNames(lapply(seq_len(z$npar), function(a) {
      wa <- z$wa[[a]]
      .mvt_const_slices(z$si, -z$ca[[a]]) +
        .mvt_const_slices(z$Sa[[a]], -z$cw) +
        .mvt_scale_slices(ww, 2 * z$da[[a]]) +
        .mvt_scale_slices(mv_outer_rows(wa, z$w) + mv_outer_rows(z$w, wa),
                          2 * z$dw)
    }), distrib@params)
  }

#' @title Multivariate Student t Third Derivative in One Response and Two Parameters
#' @name distrib_grad_y_hess.MvStudentTDistrib
#'
#' @description
#' Computes \eqn{\partial^3\ell/\partial y\,\partial\theta_a\partial\theta_b},
#' one \eqn{n \times p} matrix per unordered pair of parameters. The response
#' gradient is \eqn{\ell^{(y)} = -c\,w}, so a second derivative in the
#' parameters is the ordinary product rule
#' \deqn{\partial_{ab}\ell^{(y)} = -c_{ab}w - c_aw_b - c_bw_a - c\,w_{ab},}
#' with every piece from [mvt_dpieces()]. For a gaussian the pair of two
#' location parameters is exactly zero, its response gradient being linear in
#' the location; the response gradient of the t is not, and these pairs are not
#' zero in general.
#'
#' @param distrib An [MvStudentTDistrib] object, from
#'   [mvstudent_t1_distrib()].
#' @param y An \eqn{n \times p} numeric matrix of observations. A vector of
#'   length \eqn{p} is read as a single observation.
#' @param theta A named list of parameters, each component a single number.
#' @param scale One of `"parameter"` (the default) or `"link"`. The method
#'   returns the parameter scale and the generic applies the link scale to its
#'   result. The two differ here.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A named list of \eqn{n \times p} numeric matrices, keyed as
#'   [`hess_names(distrib@params)`][hess_names].
#'
#' @section Notation:
#' \eqn{\nu} is the degrees of freedom, \eqn{c} the weight,
#' \eqn{w = \Sigma^{-1}(y-\mu)}, and a subscript \eqn{a} or \eqn{ab} denotes a
#' derivative in the parameters named.
#'
#' @seealso [distrib_cross_y.MvStudentTDistrib()] for the order below,
#'   [distrib_hess_y_hess.MvStudentTDistrib()] for the sibling with two
#'   response indices, [mvt_dpieces()] for the assembly, and
#'   [distrib_grad_y_hess()] for the generic.
#'
#' @examples
#' d <- mvstudent_t1_distrib(2)
#' theta <- list(mu1 = 0.5, mu2 = -0.3, sigma_log_L1 = 0.1,
#'               sigma_log_L2 = -0.2, sigma_L2.1 = 0.4, nu = 6)
#' set.seed(1)
#' y <- distrib_rng(d, 4, theta)
#'
#' gh <- distrib_grad_y_hess(d, y, theta)
#' dim(gh$mu1_mu2)
#'
#' # A gaussian's two-location pair is exactly zero and this one is not.
#' gh$mu1_mu2
#' distrib_grad_y_hess(mvgaussian1_distrib(2), y, theta[1:5])$mu1_mu2
#'
#' # Against a second difference of the response gradient. The reference
#' # amplifies rounding by h^-2, so 1e-6 on a quantity of order 1 is its own
#' # accuracy rather than a disagreement.
#' h <- 1e-5
#' f <- function(a, b) {
#'   t2 <- theta; t2$mu1 <- t2$mu1 + a; t2$mu2 <- t2$mu2 + b
#'   distrib_grad_y(d, y, t2)
#' }
#' c(gap = max(abs(gh$mu1_mu2 -
#'                 (f(h, h) - f(h, -h) - f(-h, h) + f(-h, -h)) / (4 * h * h))),
#'   scale = max(abs(gh$mu1_mu2)))
#'
#' @keywords internal
S7::method(distrib_grad_y_hess, MvStudentTDistrib) <-
  function(distrib, y, theta, scale = c("parameter", "link"), ...) {
    y <- as_mv_matrix(distrib, y)
    z <- mvt_dpieces(distrib, y, theta)
    nm <- distrib@params
    stats::setNames(lapply(hess_names(nm), function(k) {
      ij <- hess_pairs(nm)[[k]]
      a <- match(ij[1L], nm)
      b <- match(ij[2L], nm)
      ab <- z$pair(a, b)
      -ab$c * z$w - z$ca[[a]] * z$wa[[b]] - z$ca[[b]] * z$wa[[a]] -
        z$cw * ab$w
    }), hess_names(nm))
  }

#' @title Multivariate Student t Fourth Derivative in Two Responses and Two Parameters
#' @name distrib_hess_y_hess.MvStudentTDistrib
#'
#' @description
#' Computes \eqn{\partial^4\ell/\partial y\,\partial y^\top
#' \partial\theta_a\partial\theta_b}, one \eqn{p \times p \times n} array per
#' unordered pair of parameters. It is the expansion of
#' \eqn{M = \ell^{(yy)} = -c\,\Sigma^{-1} + 2d\,ww^\top} carried one order
#' further than [distrib_cross2_y.MvStudentTDistrib()]:
#' \deqn{\partial_{ab} M = -c_{ab}\Sigma^{-1} - c_a\partial_b\Sigma^{-1}
#'   - c_b\partial_a\Sigma^{-1} - c\,\partial_{ab}\Sigma^{-1}
#'   + 2d_{ab}\,ww^\top + 2d_a S(w_b, w) + 2d_b S(w_a, w)
#'   + 2d\left\{S(w_{ab}, w) + S(w_a, w_b)\right\},}
#' with \eqn{S(u, v) = uv^\top + vu^\top} and every piece from
#' [mvt_dpieces()]. For a gaussian every pair naming a location parameter is
#' exactly zero; for the t these pairs are not zero in general.
#'
#' @param distrib An [MvStudentTDistrib] object, from
#'   [mvstudent_t1_distrib()].
#' @param y An \eqn{n \times p} numeric matrix of observations. A vector of
#'   length \eqn{p} is read as a single observation.
#' @param theta A named list of parameters, each component a single number.
#' @param scale One of `"parameter"` (the default) or `"link"`. The method
#'   returns the parameter scale and the generic applies the link scale to its
#'   result. The two differ here.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A named list of \eqn{p \times p \times n} numeric arrays, keyed as
#'   [`hess_names(distrib@params)`][hess_names]. Every slice is symmetric, the
#'   derivative of a symmetric matrix being symmetric.
#'
#' @section Notation:
#' \eqn{\Sigma} is the scale matrix, \eqn{c} and \eqn{d} the two weights,
#' \eqn{w = \Sigma^{-1}(y-\mu)}, and a subscript \eqn{a} or \eqn{ab} denotes a
#' derivative in the parameters named.
#'
#' @seealso [distrib_cross2_y.MvStudentTDistrib()] for the same derivative one
#'   parameter down, [distrib_grad_y_hess.MvStudentTDistrib()] for the sibling
#'   with one response index, [mvt_dpieces()] for the assembly, and
#'   [distrib_hess_y_hess()] for the generic.
#'
#' @examples
#' d <- mvstudent_t1_distrib(2)
#' theta <- list(mu1 = 0.5, mu2 = -0.3, sigma_log_L1 = 0.1,
#'               sigma_log_L2 = -0.2, sigma_L2.1 = 0.4, nu = 6)
#' set.seed(1)
#' y <- distrib_rng(d, 4, theta)
#'
#' hh <- distrib_hess_y_hess(d, y, theta)
#' dim(hh$nu_nu)
#'
#' # Every slice is symmetric.
#' all(vapply(seq_len(4), function(i) isSymmetric(hh$nu_nu[, , i]), TRUE))
#'
#' # Against a second difference of the response Hessian, with the reference's
#' # own h^-2 amplification beside it.
#' h <- 1e-5
#' f <- function(a, b) {
#'   t2 <- theta
#'   t2$sigma_log_L1 <- t2$sigma_log_L1 + a
#'   t2$sigma_L2.1 <- t2$sigma_L2.1 + b
#'   distrib_hess_y(d, y, t2)
#' }
#' key <- "sigma_log_L1_sigma_L2.1"
#' c(gap = max(abs(hh[[key]] -
#'                 (f(h, h) - f(h, -h) - f(-h, h) + f(-h, -h)) / (4 * h * h))),
#'   scale = max(abs(hh[[key]])))
#'
#' @keywords internal
S7::method(distrib_hess_y_hess, MvStudentTDistrib) <-
  function(distrib, y, theta, scale = c("parameter", "link"), ...) {
    y <- as_mv_matrix(distrib, y)
    z <- mvt_dpieces(distrib, y, theta)
    nm <- distrib@params
    ww <- mv_outer_rows(z$w, z$w)
    stats::setNames(lapply(hess_names(nm), function(k) {
      ij <- hess_pairs(nm)[[k]]
      a <- match(ij[1L], nm)
      b <- match(ij[2L], nm)
      ab <- z$pair(a, b)
      wa <- z$wa[[a]]
      wb <- z$wa[[b]]
      sym <- function(u, v) mv_outer_rows(u, v) + mv_outer_rows(v, u)
      .mvt_const_slices(z$si, -ab$c) +
        .mvt_const_slices(z$Sa[[b]], -z$ca[[a]]) +
        .mvt_const_slices(z$Sa[[a]], -z$ca[[b]]) +
        .mvt_const_slices(ab$S, -z$cw) +
        .mvt_scale_slices(ww, 2 * ab$d) +
        .mvt_scale_slices(sym(wb, z$w), 2 * z$da[[a]]) +
        .mvt_scale_slices(sym(wa, z$w), 2 * z$da[[b]]) +
        .mvt_scale_slices(sym(ab$w, z$w) + sym(wa, wb), 2 * z$dw)
    }), hess_names(nm))
  }

#' @title Scale the Slices of an Array
#'
#' @description
#' Multiplies the \eqn{i}th slice of a \eqn{p \times p \times n} array by
#' `v[i]`, which is how a term whose matrix part is the same at every
#' observation but whose scalar part is not enters an assembled derivative.
#' The multiplication is one vectorized recycle over the flattened array, so
#' nothing loops over the observations.
#'
#' @param arr A \eqn{p \times p \times n} numeric array.
#' @param v A numeric vector of length \eqn{n}. A shorter vector recycles
#'   silently, as the arithmetic underneath would.
#'
#' @return A \eqn{p \times p \times n} numeric array with the same dimensions
#'   as `arr`.
#'
#' @seealso [.mvt_const_slices()] for the case where the matrix is what stays
#'   constant, [mv_outer_rows()] for the array these are combined with, and
#'   [distrib_cross2_y.MvStudentTDistrib()] for the consumer.
#'
#' @examples
#' arr <- array(1, c(2, 2, 3))
#' distributions7:::.mvt_scale_slices(arr, c(1, 10, 100))[1, 1, ]
#'
#' # Slice by slice it is the scalar times the slice.
#' a <- array(rnorm(12), c(2, 2, 3))
#' s <- distributions7:::.mvt_scale_slices(a, c(2, -1, 0.5))
#' all.equal(s[, , 2], -a[, , 2])
#'
#' @keywords internal
.mvt_scale_slices <- function(arr, v) {
  d <- dim(arr)
  arr * rep(v, each = d[1L] * d[2L])
}

#' @title Repeat a Constant Matrix Across Slices
#'
#' @description
#' Builds the \eqn{p \times p \times n} array whose \eqn{i}th slice is
#' `v[i] * m`, the contribution of a term with a fixed matrix and a
#' per-observation scalar. It is the mirror of
#' [.mvt_scale_slices()], where the array varies and the scalar is applied to
#' it; here the matrix is the same everywhere and only the scalar moves.
#'
#' @param m A \eqn{p \times p} numeric matrix.
#' @param v A numeric vector of length \eqn{n}, which sets the third dimension
#'   of the result.
#'
#' @return A \eqn{p \times p \times n} numeric array, with `n = length(v)`.
#'
#' @seealso [.mvt_scale_slices()] for the mirror case and
#'   [distrib_hess_y_hess.MvStudentTDistrib()] for the consumer.
#'
#' @examples
#' m <- matrix(c(1, 2, 2, 4), 2, 2)
#' cs <- distributions7:::.mvt_const_slices(m, c(1, -1))
#' dim(cs)
#' cs[, , 2]
#'
#' # Every slice is a multiple of the one matrix.
#' all.equal(cs[, , 1], m)
#' all.equal(cs[, , 2], -m)
#'
#' @keywords internal
.mvt_const_slices <- function(m, v) {
  p <- nrow(m)
  array(as.numeric(m), c(p, p, length(v))) * rep(v, each = p * p)
}


#' @title Multivariate Student t Response Hessian
#' @name distrib_hess_y.MvStudentTDistrib
#'
#' @description
#' Computes the second derivative of the log-density in the response,
#' \deqn{\frac{\partial^2 \ell}{\partial y\, \partial y^\top}
#'   = -c\,\Sigma^{-1} + 2d\,ww^\top, \qquad
#'   c = \frac{\nu+p}{\nu+q}, \quad d = \frac{\nu+p}{(\nu+q)^2},}
#' one \eqn{p \times p} matrix per observation. The gaussian's response Hessian
#' is one matrix for the whole sample; this one moves with the observation, and
#' far enough out its rank-one term makes the curvature positive along the
#' direction of the residual. The rank-one term is formed as
#' \eqn{2(\nu+p)(w/s)(w/s)^\top} with \eqn{s = \nu + q} (see [mvt_yscaled()]),
#' so it neither underflows nor overflows for a distant observation.
#'
#' @param distrib An [MvStudentTDistrib] object, from
#'   [mvstudent_t1_distrib()].
#' @param y An \eqn{n \times p} numeric matrix of observations. A vector of
#'   length \eqn{p} is read as a single observation.
#' @param theta A named list of parameters, each component a single number.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A \eqn{p \times p \times n} numeric array, slice \eqn{i} holding
#'   \eqn{\partial^2\ell_i/\partial y_i\partial y_i^\top}. Compare
#'   [distrib_hess_y.MvGaussianDistrib()], which returns a
#'   \eqn{p \times p} matrix and no third dimension at all.
#'
#' @section Notation:
#' \eqn{\Sigma} is the scale matrix, \eqn{\nu} the degrees of freedom, \eqn{p}
#' the dimension, \eqn{q} the squared Mahalanobis distance, \eqn{c} and
#' \eqn{d} the two weights and \eqn{w = \Sigma^{-1}(y-\mu)}.
#'
#' @seealso [distrib_grad_y.MvStudentTDistrib()] for the first derivative,
#'   [distrib_cross2_y.MvStudentTDistrib()] for its derivative in the
#'   parameters, and [distrib_hess_y()] for the generic.
#'
#' @examples
#' d <- mvstudent_t1_distrib(2)
#' theta <- list(mu1 = 0.5, mu2 = -0.3, sigma_log_L1 = 0.1,
#'               sigma_log_L2 = -0.2, sigma_L2.1 = 0.4, nu = 6)
#' set.seed(1)
#' y <- distrib_rng(d, 4, theta)
#'
#' hy <- distrib_hess_y(d, y, theta)
#' dim(hy)
#'
#' # Against a numerical Hessian at one observation.
#' max(abs(hy[, , 1] -
#'         numDeriv::hessian(function(z) distrib_pdf(d, z, theta, log = TRUE),
#'                           y[1, ])))
#'
#' # Far out along one coordinate the curvature turns positive, where a
#' # gaussian's never does.
#' far <- cbind(0.5 + c(1, 3, 12), -0.3)
#' vapply(1:3, function(i) distrib_hess_y(d, far, theta)[1, 1, i], numeric(1))
#' distrib_hess_y(mvgaussian1_distrib(2), far, theta[1:5])[1, 1]
#'
#' @keywords internal
S7::method(distrib_hess_y, MvStudentTDistrib) <- function(distrib, y, theta, ...) {
  y <- as_mv_matrix(distrib, y)
  pc <- mvt_pieces(distrib, theta)
  z <- mvt_yscaled(y, pc)
  n <- nrow(y)
  p <- pc$p
  N <- pc$nu + p
  out <- array(0, dim = c(p, p, n))
  for (i in seq_len(n)) {
    # -c Sigma^-1 + 2 N w w' / s^2, the second term formed from w/s
    out[, , i] <- -z$cw[i] * pc$sigma_inv + 2 * N * tcrossprod(z$ws[i, ])
  }
  out
}


#' @title Multivariate Student t Third and Fourth Response Derivatives
#' @name distrib_deriv3_y.MvStudentTDistrib
#'
#' @description
#' Compute the third and fourth derivatives of the log-density in the
#' response, one array per observation. With \eqn{w = \Sigma^{-1}(y-\mu)},
#' \eqn{q = (y-\mu)^\top w}, \eqn{s = \nu + q}, \eqn{N = \nu + p} and
#' \eqn{\Omega = \Sigma^{-1}}, the log-density is
#' \eqn{-\tfrac{N}{2}\log s} plus terms free of \eqn{y}, and
#' \deqn{\ell_{ijk} = \frac{2N}{s^2}\left(\Omega_{ij}w_k + \Omega_{ik}w_j +
#'   \Omega_{jk}w_i\right) - \frac{8N}{s^3}\,w_iw_jw_k,}
#' \deqn{\ell_{ijkl} = \frac{2N}{s^2}\left(\Omega_{ij}\Omega_{kl} +
#'   \Omega_{ik}\Omega_{jl} + \Omega_{il}\Omega_{jk}\right)
#'   - \frac{8N}{s^3}\sum_{(xy)(zt)}\Omega_{xy}w_zw_t
#'   + \frac{48N}{s^4}\,w_iw_jw_kw_l,}
#' where the middle sum runs over the six ways of choosing the pair
#' \eqn{(x, y)} from \eqn{\{i, j, k, l\}}, the other two indices taking the
#' \eqn{w}s. Each follows from the order below by \eqn{\partial s/\partial y =
#' 2w} and \eqn{\partial w/\partial y = \Omega}. Every term is evaluated in
#' \eqn{w/s} and \eqn{1/s} (see [mvt_yscaled()]), so no power of \eqn{s} is
#' formed and the arrays stay finite for an observation as far as \eqn{10^{300}}
#' from the location.
#'
#' @details
#' A multivariate Student t prior over a block of coefficients has a Hessian
#' that moves with the coefficients, so the derivative of a marginal
#' criterion in the hyperparameters reads these arrays through
#' [penalties7::penalty_dhessian_beta()] and
#' [penalties7::penalty_d2hessian_beta()].
#'
#' @param distrib An [MvStudentTDistrib] object, from
#'   [mvstudent_t1_distrib()] or [mvstudent_t2_distrib()].
#' @param y An \eqn{n \times p} numeric matrix of observations. A vector of
#'   length \eqn{p} is read as a single observation.
#' @param theta A named list of parameters, each component a single number.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return `distrib_deriv3_y()` a \eqn{p \times p \times p \times n} numeric
#'   array, `distrib_deriv4_y()` a \eqn{p \times p \times p \times p \times n}
#'   one; slice \eqn{i} in the last dimension belongs to row \eqn{i} of `y`.
#'
#' @seealso [distrib_hess_y.MvStudentTDistrib()] for the order below,
#'   [distrib_cross3_y.MvStudentTDistrib()] for the derivative in the
#'   parameters, and [distrib_deriv3_y()] for the generics.
#'
#' @examples
#' d <- mvstudent_t1_distrib(2)
#' theta <- list(mu1 = 0.5, mu2 = -0.3, sigma_log_L1 = 0.1,
#'               sigma_log_L2 = -0.2, sigma_L2.1 = 0.4, nu = 6)
#' set.seed(1)
#' y <- distrib_rng(d, 3, theta)
#' d3 <- distrib_deriv3_y(d, y, theta)
#' dim(d3)
#'
#' # Against a difference of the response Hessian along the first coordinate.
#' h <- 1e-5
#' yp <- y; yp[, 1] <- yp[, 1] + h
#' ym <- y; ym[, 1] <- ym[, 1] - h
#' max(abs(d3[, , 1, ] -
#'         (distrib_hess_y(d, yp, theta) - distrib_hess_y(d, ym, theta)) / (2 * h)))
#'
#' @keywords internal
S7::method(distrib_deriv3_y, MvStudentTDistrib) <- function(distrib, y, theta,
                                                            ...) {
  y <- as_mv_matrix(distrib, y)
  pc <- mvt_pieces(distrib, theta)
  mvt_y_tensor(mvt_yscaled(y, pc), pc, 3L)
}

#' @rdname distrib_deriv3_y.MvStudentTDistrib
#' @name distrib_deriv4_y.MvStudentTDistrib
#' @keywords internal
S7::method(distrib_deriv4_y, MvStudentTDistrib) <- function(distrib, y, theta,
                                                            ...) {
  y <- as_mv_matrix(distrib, y)
  pc <- mvt_pieces(distrib, theta)
  mvt_y_tensor(mvt_yscaled(y, pc), pc, 4L)
}

#' @title Multivariate Student t Third Response Derivative in a Parameter
#' @name distrib_cross3_y.MvStudentTDistrib
#'
#' @description
#' Computes \eqn{\partial^4\ell/\partial y_i\partial y_j\partial y_k\,
#' \partial\theta_a}, one \eqn{p \times p \times p \times n} array per
#' parameter. With the notation of
#' [distrib_deriv3_y.MvStudentTDistrib()], write
#' \eqn{T_{ijk} = \Omega_{ij}w_k + \Omega_{ik}w_j + \Omega_{jk}w_i}. The
#' derivative of
#' \eqn{\ell_{ijk} = 2Ns^{-2}T_{ijk} - 8Ns^{-3}w_iw_jw_k} in \eqn{\theta_a} is
#' \deqn{\left(\frac{2N_a}{s^2} - \frac{4Ns_a}{s^3}\right)T_{ijk}
#'   + \frac{2N}{s^2}\,\partial_aT_{ijk}
#'   - \left(\frac{8N_a}{s^3} - \frac{24Ns_a}{s^4}\right)w_iw_jw_k
#'   - \frac{8N}{s^3}\,\partial_a(w_iw_jw_k),}
#' with \eqn{N_a}, \eqn{s_a}, \eqn{\partial_aw} and \eqn{\partial_a\Omega}
#' from [mvt_dpieces()]: \eqn{N_a} is one for \eqn{\nu} and zero otherwise.
#'
#' @param distrib An [MvStudentTDistrib] object.
#' @param y An \eqn{n \times p} numeric matrix of observations.
#' @param theta A named list of parameters, each component a single number.
#' @param scale One of `"parameter"` (the default) or `"link"`. The method
#'   returns the parameter scale and the generic applies the link scale to its
#'   result.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A named list of \eqn{p \times p \times p \times n} numeric arrays,
#'   one per parameter, in `distrib@params` order.
#'
#' @seealso [distrib_deriv3_y.MvStudentTDistrib()], whose derivative in the
#'   parameters this is, and [distrib_cross3_y()] for the generic.
#'
#' @examples
#' d <- mvstudent_t1_distrib(2)
#' theta <- list(mu1 = 0.5, mu2 = -0.3, sigma_log_L1 = 0.1,
#'               sigma_log_L2 = -0.2, sigma_L2.1 = 0.4, nu = 6)
#' set.seed(1)
#' y <- distrib_rng(d, 3, theta)
#' c3 <- distrib_cross3_y(d, y, theta)
#'
#' # Against a difference of the third response derivative in nu.
#' h <- 1e-5
#' tp <- theta; tp$nu <- 6 + h
#' tm <- theta; tm$nu <- 6 - h
#' max(abs(c3$nu - (distrib_deriv3_y(d, y, tp) -
#'                  distrib_deriv3_y(d, y, tm)) / (2 * h)))
#'
#' @keywords internal
S7::method(distrib_cross3_y, MvStudentTDistrib) <-
  function(distrib, y, theta, scale = c("parameter", "link"), ...) {
    y <- as_mv_matrix(distrib, y)
    z <- mvt_dpieces(distrib, y, theta)
    N <- z$cw * z$s
    stats::setNames(lapply(seq_len(z$npar), function(a) {
      Na <- if (a == z$npar) 1 else 0
      mvt_y_tensor3_d(z, N[1L], Na, z$sa[[a]], z$wa[[a]], z$Sa[[a]])
    }), distrib@params)
  }

#' @title Response Tensors of the Multivariate Student t
#'
#' @description
#' `mvt_y_tensor()` assembles the third or fourth response derivative of
#' [distrib_deriv3_y.MvStudentTDistrib()] from the scaled quantities of
#' [mvt_yscaled()], one observation at a time: every term is written in
#' \eqn{w/s} and \eqn{1/s}, so no power of \eqn{s} is formed.
#' `mvt_y_tensor3_d()` assembles the derivative of the third one in a
#' parameter, given that parameter's derivatives of \eqn{N}, \eqn{s}, \eqn{w}
#' and \eqn{\Omega}.
#'
#' @param z For `mvt_y_tensor()` the result of [mvt_yscaled()]; for
#'   `mvt_y_tensor3_d()` the result of [mvt_dpieces()].
#' @param pc The result of [mvt_pieces()].
#' @param order `3L` or `4L`.
#' @param N The single number \eqn{\nu + p}.
#' @param Na The derivative of \eqn{N} in the parameter.
#' @param sa A numeric vector of length \eqn{n}, the derivative of \eqn{s}.
#' @param wa An \eqn{n \times p} matrix, the derivative of \eqn{w}.
#' @param Oa A \eqn{p \times p} matrix, the derivative of \eqn{\Omega}.
#'
#' @return A numeric array with \eqn{n} as its last dimension.
#'
#' @seealso [distrib_deriv3_y.MvStudentTDistrib()],
#'   [distrib_cross3_y.MvStudentTDistrib()].
#'
#' @examples
#' d <- mvstudent_t1_distrib(2)
#' theta <- list(mu1 = 0, mu2 = 0, sigma_log_L1 = 0, sigma_log_L2 = 0,
#'               sigma_L2.1 = 0, nu = 4)
#' pc <- distributions7:::mvt_pieces(d, theta)
#' z <- distributions7:::mvt_yscaled(matrix(c(1, -1), 1, 2), pc)
#' dim(distributions7:::mvt_y_tensor(z, pc, 3L))
#'
#' @keywords internal
mvt_y_tensor <- function(z, pc, order) {
  p <- pc$p
  n <- length(z$is)
  N <- pc$nu + p
  O <- pc$sigma_inv
  out <- array(0, dim = c(rep(p, order), n))
  # every term written in ws = w/s and is = 1/s, so that neither w nor s
  # appears raised to a power: 2N/s^2 O w = 2N is (O ws), 8N/s^3 w w w =
  # 8N ws ws ws, and so on
  for (i in seq_len(n)) {
    w <- z$ws[i, ]
    is <- z$is[i]
    if (order == 3L) {
      T3 <- outer(O, w) + aperm(outer(O, w), c(1L, 3L, 2L)) +
        aperm(outer(O, w), c(3L, 1L, 2L))
      out[, , , i] <- 2 * N * is * T3 - 8 * N * outer(outer(w, w), w)
    } else {
      OO <- outer(O, O)
      Q <- OO + aperm(OO, c(1L, 3L, 2L, 4L)) + aperm(OO, c(1L, 4L, 3L, 2L))
      Ow <- outer(O, outer(w, w))
      M <- Ow + aperm(Ow, c(1L, 3L, 2L, 4L)) + aperm(Ow, c(1L, 3L, 4L, 2L)) +
        aperm(Ow, c(3L, 1L, 2L, 4L)) + aperm(Ow, c(3L, 1L, 4L, 2L)) +
        aperm(Ow, c(3L, 4L, 1L, 2L))
      w4 <- outer(outer(w, w), outer(w, w))
      out[, , , , i] <- 2 * N * is^2 * Q - 8 * N * is * M + 48 * N * w4
    }
  }
  out
}

#' @rdname mvt_y_tensor
#' @keywords internal
mvt_y_tensor3_d <- function(z, N, Na, sa, wa, Oa) {
  p <- z$p
  n <- z$n
  O <- z$si
  sym3 <- function(A) A + aperm(A, c(1L, 3L, 2L)) + aperm(A, c(3L, 1L, 2L))
  out <- array(0, dim = c(p, p, p, n))
  for (i in seq_len(n)) {
    w <- z$w[i, ]
    wd <- wa[i, ]
    s <- z$s[i]
    sd <- sa[i]
    T3 <- sym3(outer(O, w))
    dT3 <- sym3(outer(Oa, w) + outer(O, wd))
    ww <- outer(outer(w, w), w)
    dww <- outer(outer(wd, w), w) + outer(outer(w, wd), w) +
      outer(outer(w, w), wd)
    out[, , , i] <- (2 * Na / s^2 - 4 * N * sd / s^3) * T3 +
      2 * N / s^2 * dT3 - (8 * Na / s^3 - 24 * N * sd / s^4) * ww -
      8 * N / s^3 * dww
  }
  out
}


#' @title Mean of a Multivariate Student t
#' @name mean.MvStudentTDistrib
#'
#' @description
#' Returns \eqn{\mathbb{E}[Y] = \mu} for \eqn{\nu > 1}, and a vector of `NaN`
#' otherwise. For one degree of freedom or fewer the defining integral does not
#' converge: the density decays like \eqn{\lVert y\rVert^{-(\nu+p)}} and the
#' first absolute moment integrates \eqn{\lVert y\rVert^{1-\nu-p}} over a shell
#' of surface \eqn{\lVert y\rVert^{p-1}}, which is finite only for \eqn{\nu >
#' 1}. The location \eqn{\mu} is a parameter and exists at every \eqn{\nu},
#' while the mean is a moment; for \eqn{\nu \le 1} this method returns `NaN` and
#' [mv_location.MvStudentTDistrib()] returns \eqn{\mu}.
#'
#' @param x An [MvStudentTDistrib] object, from [mvstudent_t1_distrib()].
#' @param theta A named list of parameters, each component a single number.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of length \eqn{p}, named `v1`, ..., `vp`, holding
#'   the location for \eqn{\nu > 1} and `NaN` throughout for \eqn{\nu \le 1}.
#'
#' @seealso [variance.MvStudentTDistrib()] for the second moment, which needs
#'   \eqn{\nu > 2}, [mv_location.MvStudentTDistrib()] for the location, which
#'   exists at every \eqn{\nu}, and [base::mean()] for the generic.
#'
#' @examples
#' d <- mvstudent_t1_distrib(2)
#' theta <- list(mu1 = 0.5, mu2 = -0.3, sigma_log_L1 = 0.1,
#'               sigma_log_L2 = -0.2, sigma_L2.1 = 0.4, nu = 6)
#'
#' mean(d, theta)
#'
#' # Below one degree of freedom the mean does not exist, while the location
#' # and the density do.
#' t2 <- theta; t2$nu <- 0.8
#' mean(d, t2)
#' mv_location(d, t2)
#' distrib_pdf(d, c(0, 0), t2)
#'
#' @keywords internal
S7::method(mean, MvStudentTDistrib) <- function(x, theta, ...) {
  mu <- mv_location(x, theta)
  nu <- mv_flat_theta(x, align_theta(x, theta))[[x@n_params]]
  if (nu <= 1) mu[] <- NaN
  mu
}


#' @title Covariance of a Multivariate Student t
#' @name variance.MvStudentTDistrib
#'
#' @description
#' Returns \eqn{\operatorname{Var}(Y) = \nu\Sigma/(\nu-2)} for \eqn{\nu > 2},
#' and a matrix of `Inf` otherwise. The scale matrix \eqn{\Sigma} is what the
#' parametrization carries and exists at every admissible \eqn{\nu}; the
#' covariance is a moment, and the two differ by that factor, which is 3 at
#' \eqn{\nu = 3} and approaches 1 as the tail lightens. [mv_sigma()] returns
#' the first and this returns the second.
#'
#' @param x An [MvStudentTDistrib] object, from [mvstudent_t1_distrib()].
#' @param theta A named list of parameters, each component a single number.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A \eqn{p \times p} numeric matrix with both dimnames `v1`, ...,
#'   `vp`: the covariance for \eqn{\nu > 2}, and `Inf` in every entry for
#'   \eqn{\nu \le 2}, the boundary \eqn{\nu = 2} included.
#'
#' @seealso [mv_sigma.MvStudentTDistrib()] for the scale matrix,
#'   [mean.MvStudentTDistrib()] for the first moment, which needs only
#'   \eqn{\nu > 1}, and [variance()] for the generic.
#'
#' @examples
#' d <- mvstudent_t1_distrib(2)
#' theta <- list(mu1 = 0.5, mu2 = -0.3, sigma_log_L1 = 0.1,
#'               sigma_log_L2 = -0.2, sigma_L2.1 = 0.4, nu = 6)
#'
#' variance(d, theta)
#' mv_sigma(d, theta)
#' all.equal(variance(d, theta), (6 / 4) * mv_sigma(d, theta))
#'
#' # The factor grows without bound as nu falls to 2.
#' vapply(c(10, 3, 2.5, 2.05, 2), function(nu) {
#'   t2 <- theta; t2$nu <- nu
#'   variance(d, t2)[1, 1]
#' }, numeric(1))
#'
#' # The correlations are the same: a positive multiple cancels out of them.
#' cv <- function(m) m[1, 2] / sqrt(m[1, 1] * m[2, 2])
#' c(scale = cv(mv_sigma(d, theta)), covariance = cv(variance(d, theta)))
#'
#' @keywords internal
S7::method(variance, MvStudentTDistrib) <- function(x, theta, ...) {
  sg <- mv_sigma(x, theta)
  nu <- mv_flat_theta(x, align_theta(x, theta))[[x@n_params]]
  if (nu <= 2) {
    sg[] <- Inf
    return(sg)
  }
  sg * nu / (nu - 2)
}


#' @title Marginal of a Multivariate Student t
#' @name mv_marginal.MvStudentTDistrib
#'
#' @description
#' Returns the marginal law of a subset of coordinates, which for this family is
#' again a Student t: the subvector of the location, the corresponding block of
#' the scale matrix, and the same degrees of freedom. The mixing variable of the
#' scale-mixture representation is shared by every coordinate, so it survives
#' integrating any of them out. The returned object is a new
#' [mvstudent_t1_distrib()] of the reduced dimension when `distrib` is a
#' [MvStudentT1Distrib], and a new [mvstudent_t2_distrib()] when it is a
#' [MvStudentT2Distrib], both on an unstructured log-Cholesky matrix. The free
#' values are computed by `parameters7::param_free()` from the block of the
#' scale matrix (from its inverse for the second family), so they are not a
#' subset of the original's.
#'
#' @param distrib An [MvStudentTDistrib] object, from [mvstudent_t1_distrib()]
#'   or [mvstudent_t2_distrib()].
#' @param theta A named list of parameters, each component a single number.
#' @param which An integer vector of distinct coordinates to keep, each between
#'   1 and \eqn{p}, already checked by the generic [mv_marginal()]. The marginal
#'   follows the order given.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A named list with `distrib`, an `MvStudentTDistrib` of dimension
#'   `length(which)` of the parent's subclass, and `theta`, its parameters as a
#'   named list.
#'
#' @section Notation:
#' \eqn{\mu} is the location, \eqn{\Sigma} the scale matrix, \eqn{\nu} the
#' degrees of freedom and \eqn{p} the dimension of the full law.
#'
#' @seealso [mv_sigma.MvStudentTDistrib()] for the matrix it takes a block of,
#'   [plot.multivariate_distrib()], whose diagonal panels are these marginals,
#'   and [mv_marginal()] for the generic.
#'
#' @examples
#' d <- mvstudent_t1_distrib(3)
#' theta <- as.list(stats::setNames(
#'   c(1, -2, 0.5, 0.1, -0.2, 0.3, 0.4, -0.1, 0.2, 6), d@params))
#'
#' m <- mv_marginal(d, theta, c(1, 3))
#' m$distrib@n_dim
#'
#' # The degrees of freedom are unchanged, and the scale block is the
#' # submatrix of the full one.
#' c(full = theta$nu, marginal = m$theta$nu)
#' all.equal(mv_sigma(m$distrib, m$theta), mv_sigma(d, theta)[c(1, 3), c(1, 3)],
#'           check.attributes = FALSE)
#'
#' # A single coordinate is the univariate t, against stats::dt scaled.
#' m1 <- mv_marginal(d, theta, 1)
#' s1 <- sqrt(mv_sigma(d, theta)[1, 1])
#' c(ours = distrib_pdf(m1$distrib, 2.4, m1$theta, log = TRUE),
#'   dt = dt((2.4 - 1) / s1, df = 6, log = TRUE) - log(s1))
#'
#' @keywords internal
S7::method(mv_marginal, MvStudentTDistrib) <- function(distrib, theta, which, ...) {
  mu <- as.numeric(mv_location(distrib, theta))[which]
  sg <- mv_sigma(distrib, theta)[which, which, drop = FALSE]
  nu <- mv_flat_theta(distrib, align_theta(distrib, theta))[[distrib@n_params]]
  # the marginal keeps the parametrization it came from: a family written in
  # the inverse scale matrix reports one, its free values read off the inverse
  # of the scale submatrix
  md <- if (isTRUE(distrib@inverted)) {
    mvstudent_t2_distrib(length(which))
  } else {
    mvstudent_t1_distrib(length(which))
  }
  eta <- parameters7::param_free(
    md@param, if (isTRUE(distrib@inverted)) solve(unname(sg)) else unname(sg)
  )
  list(
    distrib = md,
    theta = as.list(stats::setNames(c(mu, unname(eta), nu), md@params))
  )
}


#' @title The Scale Matrix of a Multivariate Student t
#' @name mv_sigma.MvStudentTDistrib
#'
#' @description
#' Returns the scale matrix \eqn{\Sigma}, assembled from the matrix
#' parametrization's free values. This is the matrix in the density, and it is
#' not the covariance: the covariance is \eqn{\nu\Sigma/(\nu-2)} and exists
#' only for \eqn{\nu > 2}, while the scale matrix exists at every admissible
#' \eqn{\nu}. [variance.MvStudentTDistrib()] returns the covariance. The
#' correlations read off either are the same, a positive multiple of a matrix
#' leaving them alone.
#'
#' @param distrib An [MvStudentTDistrib] object, from
#'   [mvstudent_t1_distrib()].
#' @param theta A named list of parameters, each component a single number. The
#'   generic requires every parameter and checks `nu` against its domain; the
#'   location components and `nu` do not enter the result.
#'
#' @return A \eqn{p \times p} symmetric positive definite numeric matrix, with
#'   both dimnames `v1`, ..., `vp`.
#'
#' @seealso [variance.MvStudentTDistrib()] for the covariance,
#'   [mv_sigma.MvGaussianDistrib()], where the two coincide,
#'   [mv_summary()], which reports the square roots of the diagonal of this matrix
#'   as `scale_sd_v1` and so on, and [mv_sigma()] for the generic.
#'
#' @examples
#' d <- mvstudent_t1_distrib(2)
#' theta <- list(mu1 = 0.5, mu2 = -0.3, sigma_log_L1 = 0.1,
#'               sigma_log_L2 = -0.2, sigma_L2.1 = 0.4, nu = 6)
#'
#' mv_sigma(d, theta)
#'
#' # Positive definite whatever the free values are.
#' eigen(mv_sigma(d, theta), only.values = TRUE)$values
#'
#' # It does not move with nu, where the covariance does.
#' vapply(c(3, 6, 60), function(nu) {
#'   t2 <- theta; t2$nu <- nu
#'   c(scale = mv_sigma(d, t2)[1, 1], covariance = variance(d, t2)[1, 1])
#' }, numeric(2))
#'
#' @keywords internal
S7::method(mv_sigma, MvStudentTDistrib) <- function(distrib, theta) {
  pc <- mvt_pieces(distrib, theta)
  nm <- paste0("v", seq_len(distrib@n_dim))
  dimnames(pc$sigma) <- list(nm, nm)
  pc$sigma
}


#' @title Random Parameters for a Multivariate Student t
#' @name generate_random_theta.MvStudentTDistrib
#'
#' @description
#' Draws a parameter vector for testing: each location uniform on \eqn{(-1, 1)},
#' each free value of the matrix parametrization uniform on \eqn{(-0.4, 0.4)},
#' and the degrees of freedom uniform on \eqn{(5, 12)}. The matrix band puts the
#' scale matrix near the identity, as it does for the gaussian. The \eqn{\nu}
#' band keeps the family heavy-tailed and above \eqn{\nu = 4}, so the covariance
#' and the fourth moment exist at every draw; the comparison of the generator
#' with the covariance in [check_distrib()] needs both.
#'
#' @param distrib An [MvStudentTDistrib] object, from
#'   [mvstudent_t1_distrib()].
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A named list of `distrib@n_params` single numbers, named and ordered
#'   as `distrib@params`.
#'
#' @seealso [generate_random_theta.MvGaussianDistrib()] for the gaussian's
#'   bands, [check_distrib()], which draws parameters this way, and
#'   [generate_random_theta()] for the generic.
#'
#' @examples
#' d <- mvstudent_t1_distrib(2)
#'
#' set.seed(11)
#' unlist(generate_random_theta(d))
#'
#' # The three bands, over 400 draws.
#' set.seed(12)
#' round(apply(replicate(400, unlist(generate_random_theta(d))), 1, range), 2)
#'
#' # Every draw keeps nu above 4, so the covariance and the fourth moment of
#' # every coordinate exist.
#' set.seed(13)
#' all(replicate(50, is.finite(variance(d, generate_random_theta(d))[1, 1])))
#'
#' @keywords internal
S7::method(generate_random_theta, MvStudentTDistrib) <- function(distrib, ...) {
  p <- distrib@n_dim
  s <- distrib@param
  as.list(stats::setNames(
    c(stats::runif(p, -1, 1), stats::runif(s@n_free, -0.4, 0.4),
      stats::runif(1, 5, 12)),
    distrib@params
  ))
}

#' @title Location of a Multivariate Student t
#' @name mv_location.MvStudentTDistrib
#'
#' @description
#' Returns the location \eqn{\mu}, the first \eqn{p} parameters read off
#' `theta` in order. The method is [mv_leading_location()], shared with the
#' gaussian. For this family the location is the center of symmetry at every
#' admissible \eqn{\nu} and the mean only for \eqn{\nu > 1}:
#' [mean.MvStudentTDistrib()] returns `NaN` for \eqn{\nu \le 1}, and this method
#' returns the location at every \eqn{\nu}.
#'
#' @param distrib An [MvStudentTDistrib] object, from
#'   [mvstudent_t1_distrib()].
#' @param theta A named list of parameters, one number each. Only the \eqn{p}
#'   location components are read.
#'
#' @return A numeric vector of length \eqn{p}, named `v1`, ..., `vp` after the
#'   coordinates of the response.
#'
#' @seealso [mean.MvStudentTDistrib()] for the moment,
#'   [mv_sigma.MvStudentTDistrib()] for the matrix, and [mv_location()] for the
#'   generic.
#'
#' @examples
#' d <- mvstudent_t1_distrib(2)
#' theta <- list(mu1 = 0.5, mu2 = -0.3, sigma_log_L1 = 0.1,
#'               sigma_log_L2 = -0.2, sigma_L2.1 = 0.4, nu = 6)
#'
#' mv_location(d, theta)
#'
#' # Above one degree of freedom the location is the mean; below it, only the
#' # location survives.
#' t2 <- theta; t2$nu <- 0.8
#' rbind(location = mv_location(d, t2), mean = mean(d, t2))
#'
#' @keywords internal
S7::method(mv_location, MvStudentTDistrib) <- mv_leading_location
