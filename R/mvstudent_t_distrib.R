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
#' is added last. It inherits from `multivariate_distrib`, so the response is
#' an \eqn{n \times p} matrix and the distribution function and the quantile
#' function are refused.
#'
#' Build one with [mvstudent_t_distrib()], which fills the properties in and
#' checks the rank of the parametrization. This page documents the raw S7
#' constructor, which validates nothing.
#'
#' @param param A \pkg{parameters7} parametrization of the SCALE matrix,
#'   inheriting from `parameters7::parameter`. Its `n_free` free values are
#'   flattened into scalar parameters of the distribution.
#' @param n_dim The dimension \eqn{p} of one observation. A single positive
#'   integer.
#' @inheritParams distrib
#'
#' @return An S7 object of class `MvStudentTDistrib`, inheriting from
#'   `multivariate_distrib` and from `distrib`. Beyond the parent's
#'   `distrib_name`, `dimension`, `bounds`, `params`,
#'   `params_interpretation`, `n_params`, `params_bounds`, `link_params`,
#'   `params_smooth` and `n_dim`, it carries `param` as supplied. Unlike
#'   [MvGaussianDistrib] it has no `inverted` property: the scale matrix is
#'   always parametrized as itself.
#'
#' @seealso [mvstudent_t_distrib()] to build one, [mvgaussian_distrib()] for
#'   the limit \eqn{\nu \to \infty}, [mv_sigma.MvStudentTDistrib()] for the
#'   scale matrix and [variance.MvStudentTDistrib()] for the covariance, which
#'   are different matrices here.
#'
#' @section Methods:
#' Registered on this class:
#'   [`distrib_cross2_y()`][distrib_cross2_y.MvStudentTDistrib],
#'   [`distrib_cross_y()`][distrib_cross_y.MvStudentTDistrib],
#'   [`distrib_expected_hessian()`][distrib_expected_hessian.MvStudentTDistrib],
#'   [`distrib_grad_y()`][distrib_grad_y.MvStudentTDistrib],
#'   [`distrib_gradient()`][distrib_gradient.MvStudentTDistrib],
#'   [`distrib_hess_y()`][distrib_hess_y.MvStudentTDistrib],
#'   [`distrib_hessian()`][distrib_hessian.MvStudentTDistrib],
#'   [`distrib_pdf()`][distrib_pdf.MvStudentTDistrib],
#'   [`distrib_rng()`][distrib_rng.MvStudentTDistrib],
#'   [`generate_random_theta()`][generate_random_theta.MvStudentTDistrib],
#'   [`mean()`][mean.MvStudentTDistrib],
#'   [`mv_location()`][mv_location.MvStudentTDistrib],
#'   [`mv_marginal()`][mv_marginal.MvStudentTDistrib],
#'   [`mv_sigma()`][mv_sigma.MvStudentTDistrib],
#'   [`variance()`][variance.MvStudentTDistrib]
#'
#' Everything else comes from [multivariate_distrib()]. The third and fourth
#' derivatives are registered in `mv_higher.R`, which the two multivariate
#' families share.
#'
#' @examples
#' d <- mvstudent_t_distrib(2)
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
  properties = list(param = parameters7::parameter)
)


#' @title Construct a Multivariate Student t Distribution
#'
#' @description
#' Builds the elliptical Student t family on \eqn{\mathbb{R}^p}: a location
#' vector of \eqn{p} free parameters, a scale matrix carried by a
#' \pkg{parameters7} parametrization, and a positive degrees-of-freedom
#' parameter. The free values of the parametrization become ordinary scalar
#' parameters, so the object answers every generic of the `distrib` contract
#' with `theta` a flat named list of numbers. The default scale matrix is
#' unstructured in the log-Cholesky parametrization, which is
#' `p * (p + 1) / 2` free values, and the default link for `nu` is the log.
#'
#' @details
#' # The scale matrix is not the covariance
#'
#' \eqn{\Sigma} is the SCALE matrix. The covariance is \eqn{\nu\Sigma/(\nu-2)}
#' where it exists, and for \eqn{\nu \le 2} it does not exist while the density
#' is perfectly well defined. [mv_sigma()] returns the scale matrix, which is
#' what the parametrization carries, and [variance()] returns the covariance,
#' which is a moment. The separation is what allows a fit to run at
#' \eqn{\nu = 1.5}, and below \eqn{\nu = 1} the mean does not exist either and
#' [base::mean()] returns `NaN`.
#'
#' # Parameters
#'
#' The location contributes `mu1`, ..., `mup`, the parametrization contributes
#' its free values under the `sigma_` prefix, and `nu` comes last. The location
#' and the matrix parameters are unconstrained and carry identity links; `nu`
#' is positive and carries `link_nu`. So, unlike [mvgaussian_distrib()], this
#' family's link scale is not its parameter scale, and `scale = "link"`
#' multiplies each `nu` component by the chain-rule factor \eqn{\nu} under the
#' default log link.
#'
#' # What it is for
#'
#' A gaussian fitted to data with a few outlying rows inflates its covariance
#' to cover them. A \eqn{t} with \eqn{\nu} estimated does not: an observation
#' at Mahalanobis distance \eqn{q} enters every derivative through the weight
#' \eqn{c = (\nu+p)/(\nu+q)}, which falls away as \eqn{1/q}. The gaussian is
#' the limit \eqn{\nu \to \infty}, where \eqn{c \equiv 1} and nothing is
#' downweighted.
#'
#' # The expected information
#'
#' Closed form, from the family's own scale mixture; see
#' [distrib_expected_hessian.MvStudentTDistrib()] for the four blocks.
#' [fit_distrib()] therefore REJECTS an `approx` argument on
#' [fisher_scoring()] for this family, as it does for any family that computes
#' its expected information exactly.
#'
#' # Reading a fit
#'
#' [mv_summary()] reports the square roots of the diagonal of the SCALE matrix
#' and the correlations. The correlations are the response's as well, a
#' positive multiple of a matrix leaving them alone, but the diagonal
#' quantities are not standard deviations of the response and are named
#' `scale_sd_v1`, ..., `scale_sd_vp` to say so.
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
#'   whole number, finite and at least 1. Anything else throws an error.
#' @param sigma A \pkg{parameters7} parametrization of the scale matrix, of
#'   dimension `n_dim` and of full rank. Defaults to
#'   `parameters7::log_cholesky(n_dim)`. A rank-deficient parametrization is
#'   rejected: the density would not normalize.
#' @param link_nu The link carrying the degrees of freedom to the
#'   unconstrained scale, a `linkfunctions7::link` object. Defaults to
#'   `linkfunctions7::log_link()`, which keeps \eqn{\nu} strictly positive at
#'   every point of the free scale.
#'
#' @return An S7 object of class [MvStudentTDistrib], with `param` the
#'   parametrization supplied. Its `params` are `mu1`, ..., `mup`, the prefixed
#'   free names, then `nu`; `n_params` is `p + param@n_free + 1`; every
#'   `params_bounds` entry is \eqn{(-\infty, \infty)} except `nu`'s, which is
#'   \eqn{(0, \infty)}; and every link is the identity except `nu`'s.
#'
#' @seealso [mvgaussian_distrib()] for the limiting family,
#'   [student_t1_distrib()] for the one-dimensional case,
#'   [mv_marginal.MvStudentTDistrib()] for a marginal, which keeps the same
#'   \eqn{\nu}, and [fit_distrib()] to estimate one.
#'
#' @examples
#' d <- mvstudent_t_distrib(2)
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
#' # Below two degrees of freedom the covariance does not exist and the
#' # density does; below one the mean goes too.
#' theta$nu <- 1.5
#' c(density = distrib_pdf(d, c(0, 0), theta), variance = variance(d, theta)[1, 1])
#' theta$nu <- 0.8
#' c(density = distrib_pdf(d, c(0, 0), theta), mean = mean(d, theta)[1])
#'
#' # The heavy tail is what the family is for: at four standard-scale units
#' # out, a t at nu = 6 puts far more mass than a gaussian of the same scale.
#' theta$nu <- 6
#' g <- mvgaussian_distrib(2)
#' c(t = distrib_pdf(d, c(4.5, -4.3), theta),
#'   gaussian = distrib_pdf(g, c(4.5, -4.3), theta[1:5]))
#'
#' @export
mvstudent_t_distrib <- function(n_dim, sigma = NULL,
                                link_nu = linkfunctions7::log_link()) {
  if (!is.numeric(n_dim) || length(n_dim) != 1L || !is.finite(n_dim) ||
    n_dim < 1 || n_dim != round(n_dim)) {
    stop("'n_dim' must be a single positive integer.", call. = FALSE)
  }
  p <- as.integer(n_dim)
  s <- if (is.null(sigma)) parameters7::log_cholesky(p) else sigma

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
  # covariance; no precision form of this family exists yet.
  free_names <- mv_prefixed_names(s@free_names, inverted = FALSE)
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

  MvStudentTDistrib(
    distrib_name = sprintf("multivariate student t [%dd, sigma=%s]", p, s@param_name),
    dimension = "multivariate",
    n_dim = p,
    bounds = c(-Inf, Inf),
    params = params,
    params_interpretation = stats::setNames(
      c(rep("location", p), rep("scale", s@n_free), "degrees of freedom"),
      params
    ),
    n_params = n_par,
    params_bounds = stats::setNames(
      c(rep(list(c(-Inf, Inf)), n_par - 1L), list(c(0, Inf))), params
    ),
    link_params = stats::setNames(
      c(rep(list(ident), n_par - 1L), list(link_nu)), params
    ),
    param = s
  )
}


#' @title The Pieces a Multivariate t Evaluates From
#'
#' @description
#' Assembles, once per call, the location, the scale matrix, its inverse, the
#' log-determinant and the degrees of freedom from a flat parameter vector,
#' together with the matrix parametrization's derivative arrays when they are
#' asked for. Every method of [MvStudentTDistrib] calls this first and works
#' from the result, so a parameter vector is unpacked and a matrix factorized
#' once per call instead of once per component.
#'
#' @details
#' The scale matrix is always parametrized as itself, so there is no inversion
#' branch here and no sign to flip on the log-determinant. That is the whole
#' difference from [mvg_pieces()], whose parametrization may carry either side.
#'
#' @param distrib An [MvStudentTDistrib] object.
#' @param theta A named list of parameters, already aligned by the generic, or
#'   any list whose components are in `distrib@params` order.
#' @param derivs Logical of length 1. When `TRUE` the first derivative arrays
#'   of the scale matrix are computed and returned as `a`. Defaults to `FALSE`.
#' @param derivs2 Logical of length 1. When `TRUE` the second derivative arrays
#'   are computed as well and returned as `a2`; the first derivatives are not
#'   computed with them, so a caller needing both sets `derivs` too. Defaults
#'   to `FALSE`.
#'
#' @return A named list with `mu` (numeric of length \eqn{p}), `eta` (the
#'   matrix parametrization's free vector), `nu` (a single number), `p`, `s`
#'   (the parametrization itself), `sigma` and `sigma_inv` (\eqn{p \times p}
#'   matrices) and `logdet` (the log-determinant of \eqn{\Sigma}), plus `a` and
#'   `a2` when asked for: lists of \eqn{p \times p} matrices, `a` one per free
#'   value and `a2` one per unordered pair.
#'
#' @section Notation:
#' \eqn{\mu} is the location, \eqn{\Sigma} the scale matrix, \eqn{\nu} the
#' degrees of freedom and \eqn{\eta} the free vector of the matrix
#' parametrization.
#'
#' @seealso [mvstudent_t_distrib()] for the family, [mvt_weights()] for the
#'   quantities computed from these, and [mvg_pieces()] for the gaussian's
#'   twin.
#'
#' @examples
#' d <- mvstudent_t_distrib(2)
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

  out <- list(
    mu = mu, eta = eta, nu = nu, p = p, s = s,
    sigma = unname(parameters7::param_value(s, eta)),
    sigma_inv = unname(parameters7::param_solve(s, eta)),
    logdet = parameters7::param_logdet(s, eta)
  )
  if (derivs || derivs2) {
    out$a <- lapply(parameters7::param_d1(s, eta), unname)
  }
  if (derivs2) {
    out$a2 <- lapply(parameters7::param_d2(s, eta), unname)
  }
  out
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
#' with \eqn{q = (y-\mu)^\top \Sigma^{-1}(y-\mu)}, row by row, and
#' exponentiates it unless `log = TRUE`. The last logarithm is taken with
#' [base::log1p()], which is the difference between a number and a loss of
#' every significant digit when \eqn{q/\nu} is small: near the center of a
#' family with many degrees of freedom, `log(1 + q/nu)` cancels its argument
#' away.
#'
#' @param distrib An [MvStudentTDistrib] object, from
#'   [mvstudent_t_distrib()].
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
#' d <- mvstudent_t_distrib(2)
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
#' g <- mvgaussian_distrib(2)
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
  q <- rowSums((r %*% pc$sigma_inv) * r)
  nu <- pc$nu
  p <- pc$p
  ld <- lgamma((nu + p) / 2) - lgamma(nu / 2) -
    (p / 2) * base::log(nu * pi) - 0.5 * pc$logdet -
    ((nu + p) / 2) * base::log1p(q / nu)
  if (log) ld else exp(ld)
}


#' @title Multivariate Student t Generator
#' @name distrib_rng.MvStudentTDistrib
#'
#' @description
#' Draws from the family through its scale-mixture representation,
#' \eqn{y = \mu + L z\sqrt{\nu/g}} with \eqn{z} a vector of independent
#' standard normals, \eqn{g \sim \chi^2_\nu} independent of them, and
#' \eqn{LL^\top = \Sigma}. A \eqn{t} is a gaussian whose precision has been
#' multiplied by a gamma variate, and that same representation is what gives
#' the family its heavy tails and its closed-form expected information.
#'
#' The draws consume `n * p` normal variates and `n` chi-squared variates from
#' R's own generators, in that order, so the stream is reproducible under
#' [base::set.seed()].
#'
#' @param distrib An [MvStudentTDistrib] object, from
#'   [mvstudent_t_distrib()].
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
#' degrees of freedom and \eqn{L} a lower Cholesky factor of \eqn{\Sigma}.
#'
#' @seealso [distrib_pdf.MvStudentTDistrib()] for the density this draws from,
#'   [variance.MvStudentTDistrib()] for the covariance the sample approaches
#'   when it exists, and [distrib_rng()] for the generic.
#'
#' @examples
#' d <- mvstudent_t_distrib(2)
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
#' # At nu = 1.5 it does not, and the sample covariance is a number that means
#' # nothing: it grows with the sample instead of settling.
#' t2 <- theta; t2$nu <- 1.5
#' set.seed(3)
#' vapply(c(2000, 20000, 200000),
#'        function(m) var(distrib_rng(d, m, t2))[1, 1], numeric(1))
#'
#' @keywords internal
S7::method(distrib_rng, MvStudentTDistrib) <- function(distrib, n, theta, ...) {
  pc <- mvt_pieces(distrib, theta)
  p <- pc$p
  l <- t(chol(pc$sigma))
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
#' The weight is the whole of the family's resistance to outliers. At
#' \eqn{q = 0} it is \eqn{(\nu+p)/\nu}, and it decays like \eqn{1/q}, so an
#' observation far from the location contributes less to every derivative
#' instead of dragging the fit towards itself. Letting \eqn{\nu \to \infty}
#' sends it to one and recovers the gaussian, where nothing is downweighted.
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
#' d <- mvstudent_t_distrib(2)
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
  q <- rowSums(r * w)
  list(r = r, w = w, q = q, cw = (pc$nu + pc$p) / (pc$nu + q))
}


#' @title Multivariate Student t Score
#' @name distrib_gradient.MvStudentTDistrib
#'
#' @description
#' Computes the first derivatives of the log-density in closed form. With
#' \eqn{w = \Sigma^{-1}(y-\mu)}, \eqn{q = (y-\mu)^\top w} and
#' \eqn{c = (\nu+p)/(\nu+q)},
#' \deqn{\partial_\mu \ell = c\,w, \qquad
#'   \partial_{\eta_k}\ell = -\tfrac{1}{2}\partial_{\eta_k}\log\lvert\Sigma\rvert
#'   + \tfrac{c}{2}\, w^\top A_k w,}
#' \deqn{\partial_\nu \ell = \tfrac{1}{2}\left[
#'   \psi\!\left(\tfrac{\nu+p}{2}\right) - \psi\!\left(\tfrac{\nu}{2}\right)
#'   - \tfrac{p}{\nu} - \log\!\left(1+\tfrac{q}{\nu}\right)
#'   + \tfrac{(\nu+p)q}{\nu(\nu+q)}\right].}
#' The gaussian's score is the limit \eqn{c \to 1}. Every observation enters
#' the location and matrix components through that one weight, and its decay
#' with \eqn{q} is the whole of the family's resistance to a distant row.
#'
#' @details
#' The \eqn{\nu} component is NOT evaluated in the form printed above. Both of
#' its bracketed pairs cancel to leading order as \eqn{\nu} grows, so it is
#' assembled as \eqn{A_p(\nu) + D(q/\nu) + (p/\nu)u/(1+u)} instead, with
#' [mvt_A()] and [mvt_D()] carrying out each cancellation exactly. The direct
#' form loses every digit at large \eqn{\nu} and changes sign at
#' \eqn{\nu = 10^9}; see [mvt_A()] for the measurements.
#'
#' @param distrib An [MvStudentTDistrib] object, from
#'   [mvstudent_t_distrib()].
#' @param y An \eqn{n \times p} numeric matrix of observations. A vector of
#'   length \eqn{p} is read as a single observation.
#' @param theta A named list of parameters, each component a single number.
#' @param scale One of `"parameter"` (the default) or `"link"`, handled by the
#'   generic before dispatch. The two differ here: under the default log link
#'   the `nu` component is multiplied by \eqn{\nu}, and the rest are unchanged.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A named list with one numeric vector of length \eqn{n} per
#'   parameter, in `distrib@params` order: \eqn{p} location components, one per
#'   free value of the matrix parametrization, then `nu`.
#'
#' @section Notation:
#' \eqn{\mu} is the location, \eqn{\Sigma} the scale matrix, \eqn{\nu} the
#' degrees of freedom, \eqn{\eta} the free vector of the matrix
#' parametrization, \eqn{A_k = \partial\Sigma/\partial\eta_k}, \eqn{q} the
#' squared Mahalanobis distance, \eqn{c} the weight,
#' \eqn{\psi} the digamma function and \eqn{\ell} the log-density of one
#' observation.
#'
#' @seealso [distrib_hessian.MvStudentTDistrib()] for the observed curvature,
#'   [distrib_expected_hessian.MvStudentTDistrib()] for the information,
#'   [mvt_weights()] for the weight, [mvt_A()] for the cancellation, and
#'   [distrib_gradient()] for the generic.
#'
#' @examples
#' d <- mvstudent_t_distrib(2)
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
#' # The location component is the weighted whitened residual, so a distant
#' # row contributes less than a near one of the same direction.
#' z <- distributions7:::mvt_weights(y, distributions7:::mvt_pieces(d, theta))
#' all.equal(cbind(g$mu1, g$mu2), z$cw * z$w, check.attributes = FALSE)
#'
#' # The nu component decays as nu^-2, computed without cancellation.
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

  dld <- parameters7::param_dlogdet(pc$s, pc$eta)
  for (k in seq_along(pc$a)) {
    out[[p + k]] <- -0.5 * dld[[k]] +
      0.5 * z$cw * rowSums((z$w %*% pc$a[[k]]) * z$w)
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
#' in \eqn{\nu} that does not involve the data, in a form whose terms carry one
#' sign so that nothing cancels.
#'
#' @details
#' # Why the direct form fails
#'
#' As \eqn{\nu} grows the multivariate t tends to the multivariate gaussian and
#' every derivative in \eqn{\nu} vanishes, so \eqn{A_p} is a difference of
#' terms agreeing to leading order:
#' \eqn{\psi(\tfrac{\nu+p}{2}) - \psi(\tfrac{\nu}{2})} is \eqn{p/\nu}, and
#' \eqn{p/\nu} is what it is subtracted from. Measured at \eqn{p = 4}, the
#' direct form is wrong by a relative \eqn{4.9\times10^{-5}} at
#' \eqn{\nu = 10^6}, by 39 per cent at \eqn{10^8}, and it CHANGES SIGN at
#' \eqn{10^9}, reading \eqn{+3.3\times10^{-16}} where the exact value is
#' \eqn{-4.0\times10^{-18}}.
#'
#' # The repair needs no series for an even dimension
#'
#' \eqn{p} is an integer dimension, so the shift between the two digamma
#' arguments is a whole number of steps of the recurrence
#' \eqn{\psi(x+1) = \psi(x) + 1/x}. For even \eqn{p} that gives a sum whose
#' terms all carry one sign:
#' \deqn{A_p(\nu) = -\sum_{j=0}^{p/2-1} \frac{4j}{\nu(\nu+2j)}.}
#' It is exactly zero at \eqn{p = 2}, where the direct form returns noise at
#' \eqn{10^{-16}}.
#'
#' For odd \eqn{p} the shift is a half-integer, so the recurrence carries the
#' quantity onto the univariate [mvt_A1()], which keeps a series of its own
#' above a measured crossover.
#'
#' @param nu The degrees of freedom, a numeric vector of any length, strictly
#'   positive. Nothing is validated: a non-positive value reaches
#'   [base::digamma()] and returns `NaN`.
#' @param p The dimension, a single positive whole number. Even and odd take
#'   different branches.
#'
#' @return A numeric vector as long as `nu`, non-positive at every \eqn{p} and
#'   exactly zero at \eqn{p = 2}.
#'
#' @section Notation:
#' \eqn{\nu} is the degrees of freedom, \eqn{p} the dimension and \eqn{\psi}
#' the digamma function.
#'
#' @seealso [mvt_A1()] for the univariate case the odd branch reduces to,
#'   [mvt_T()] for the same treatment of the trigamma pair, [mvt_D()] for the
#'   third cancellation, and [distrib_gradient.MvStudentTDistrib()] for the
#'   consumer.
#'
#' @examples
#' # At p = 2 the quantity is exactly zero and the direct form is noise.
#' direct <- function(nu, p) digamma((nu + p) / 2) - digamma(nu / 2) - p / nu
#' rbind(exact = distributions7:::mvt_A(c(3, 10, 1e6), 2),
#'       direct = direct(c(3, 10, 1e6), 2))
#'
#' # At p = 4 the two agree while the direct form has digits left, and part
#' # company at large nu, where the direct one changes sign.
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
  if (p %% 2L == 0L) {
    acc <- numeric(length(nu))
    for (j in seq_len(p %/% 2L) - 1L) {
      acc <- acc - 4 * j / (nu * (nu + 2 * j))
    }
    return(acc)
  }
  # odd p: the half-integer shift lands on the univariate A_1
  acc <- mvt_A1(nu)
  for (j in seq_len((p - 1L) %/% 2L) - 1L) {
    acc <- acc - 2 * (1 + 2 * j) / (nu * (nu + 1 + 2 * j))
  }
  acc
}


#' @title The Univariate Digamma Difference at One Degree of Freedom
#'
#' @description
#' Computes
#' \eqn{A_1(\nu) = \psi\!\left(\tfrac{\nu+1}{2}\right)
#' - \psi\!\left(\tfrac{\nu}{2}\right) - \tfrac{1}{\nu}}, the odd-dimension
#' base case of [mvt_A()]. The half-integer shift here leaves no recurrence to
#' telescope, so above \eqn{\nu = 200} the expansion
#' \eqn{1/(2\nu^2) - 1/(4\nu^4) + 1/(2\nu^6)} is used, and below it the direct
#' difference, which still has digits there.
#'
#' @details
#' The crossover of 200 is the one measured for the same expansion in this
#' package's `src/student_t.cpp`, and this is the one place in the package
#' where that series exists twice. The two copies are pinned against each other
#' in the tests. Across the switch the two branches agree to eleven figures:
#' at \eqn{\nu = 199, 200, 201} the values run
#' \eqn{1.2626\times10^{-5}}, \eqn{1.2500\times10^{-5}},
#' \eqn{1.2376\times10^{-5}} with no step.
#'
#' @param nu The degrees of freedom, a numeric vector of any length, strictly
#'   positive. The branch is taken elementwise.
#'
#' @return A numeric vector as long as `nu`, positive and decaying as
#'   \eqn{1/(2\nu^2)}.
#'
#' @section Notation:
#' \eqn{\nu} is the degrees of freedom and \eqn{\psi} the digamma function.
#'
#' @seealso [mvt_A()], which calls this on an odd dimension, and [mvt_T1()] for
#'   the trigamma twin.
#'
#' @examples
#' # No step across the crossover at nu = 200.
#' signif(distributions7:::mvt_A1(c(199, 200, 201)), 10)
#'
#' # Below it the branch is the direct difference and agrees exactly.
#' direct <- function(nu) digamma((nu + 1) / 2) - digamma(nu / 2) - 1 / nu
#' identical(distributions7:::mvt_A1(c(3, 50)), direct(c(3, 50)))
#'
#' # Above it the series is the accurate one: scaled by 2 nu^2 both approach 1,
#' # and the direct form starts to wander.
#' nu <- c(1e3, 1e4, 1e5)
#' rbind(series = distributions7:::mvt_A1(nu) * 2 * nu^2,
#'       direct = direct(nu) * 2 * nu^2)
#'
#' @keywords internal
mvt_A1 <- function(nu) {
  # psi((nu+1)/2) - psi(nu/2) - 1/nu = 1/(2 nu^2) - 1/(4 nu^4) + 1/(2 nu^6)
  # above the crossover measured for student_t.cpp's own copy of it
  out <- numeric(length(nu))
  big <- nu >= 200
  if (any(big)) {
    u2 <- 1 / nu[big]^2
    out[big] <- u2 * (0.5 + u2 * (-0.25 + u2 * 0.5))
  }
  if (any(!big)) {
    v <- nu[!big]
    out[!big] <- digamma((v + 1) / 2) - digamma(v / 2) - 1 / v
  }
  out
}


#' @title The Trigamma Difference of a Multivariate t, Without the Cancellation
#'
#' @description
#' Computes
#' \eqn{T_p(\nu) = \tfrac12\left[\psi'\!\left(\tfrac{\nu+p}{2}\right)
#' - \psi'\!\left(\tfrac{\nu}{2}\right)\right] + \tfrac{p}{\nu^2}}, the part of
#' the curvature in \eqn{\nu} that does not involve the data, written so that
#' its terms carry one sign. It is the second-order twin of [mvt_A()] and
#' cancels for the same reason: every derivative in \eqn{\nu} vanishes as the
#' family approaches the gaussian.
#'
#' @details
#' For even \eqn{p} the recurrence \eqn{\psi'(x+1) = \psi'(x) - 1/x^2}
#' telescopes into
#' \deqn{T_p(\nu) = \sum_{j=0}^{p/2-1} \frac{8j(\nu+j)}{\nu^2(\nu+2j)^2},}
#' exactly zero at \eqn{p = 2} where the direct form returns noise. For odd
#' \eqn{p} the half-integer shift carries the quantity onto the univariate
#' [mvt_T1()]. Measured at \eqn{p = 4}, the direct form is wrong by a relative
#' \eqn{4.4\times10^{-5}} at \eqn{\nu = 10^6} and by a factor of five at
#' \eqn{10^8}.
#'
#' @param nu The degrees of freedom, a numeric vector of any length, strictly
#'   positive. Nothing is validated.
#' @param p The dimension, a single positive whole number. Even and odd take
#'   different branches.
#'
#' @return A numeric vector as long as `nu`, non-negative at every \eqn{p} and
#'   exactly zero at \eqn{p = 2}.
#'
#' @section Notation:
#' \eqn{\nu} is the degrees of freedom, \eqn{p} the dimension and \eqn{\psi'}
#' the trigamma function.
#'
#' @seealso [mvt_A()] for the first-order twin, [mvt_T1()] for the univariate
#'   base case, and [distrib_hessian.MvStudentTDistrib()] for the consumer.
#'
#' @examples
#' # Exactly zero at p = 2, where the direct form is noise.
#' direct <- function(nu, p)
#'   0.5 * (trigamma((nu + p) / 2) - trigamma(nu / 2)) + p / nu^2
#' rbind(exact = distributions7:::mvt_T(c(3, 10, 1e6), 2),
#'       direct = direct(c(3, 10, 1e6), 2))
#'
#' # At p = 4 the two agree until the cancellation bites.
#' nu <- c(3, 1e4, 1e6, 1e8)
#' signif(rbind(exact = distributions7:::mvt_T(nu, 4),
#'              direct = direct(nu, 4)), 5)
#'
#' @keywords internal
mvt_T <- function(nu, p) {
  if (p %% 2L == 0L) {
    acc <- numeric(length(nu))
    for (j in seq_len(p %/% 2L) - 1L) {
      acc <- acc + 8 * j * (nu + j) / (nu^2 * (nu + 2 * j)^2)
    }
    return(acc)
  }
  # odd p: the half-integer shift is carried onto the univariate pair the
  # same way
  acc <- mvt_T1(nu)
  for (j in seq_len((p - 1L) %/% 2L) - 1L) {
    acc <- acc + 2 * (1 + 2 * j) * (2 * nu + 1 + 2 * j) /
      (nu^2 * (nu + 1 + 2 * j)^2)
  }
  acc
}


#' @title The Univariate Trigamma Difference at One Degree of Freedom
#'
#' @description
#' Computes
#' \eqn{T_1(\nu) = \tfrac12\left[\psi'\!\left(\tfrac{\nu+1}{2}\right)
#' - \psi'\!\left(\tfrac{\nu}{2}\right)\right] + \tfrac{1}{\nu^2}}, the
#' odd-dimension base case of [mvt_T()]. Above \eqn{\nu = 100} it uses the
#' expansion \eqn{-1/\nu^3 + 1/\nu^5 - 3/\nu^7}, and below it the direct
#' difference.
#'
#' @details
#' The series comes from the same duplication formula this package's
#' `src/student_t.cpp` uses, halved and with the \eqn{1/\nu^2} folded in.
#' Across the switch the two branches agree: at \eqn{\nu = 99, 100, 101} the
#' values run \eqn{-1.0305\times10^{-6}}, \eqn{-9.9990\times10^{-7}},
#' \eqn{-9.7049\times10^{-7}} with no step.
#'
#' Note the sign. \eqn{T_1} is NEGATIVE, where \eqn{T_p} for even \eqn{p} is
#' non-negative; an odd dimension assembles from this base case plus a
#' non-negative sum.
#'
#' @param nu The degrees of freedom, a numeric vector of any length, strictly
#'   positive. The branch is taken elementwise.
#'
#' @return A numeric vector as long as `nu`, negative and decaying as
#'   \eqn{-1/\nu^3}.
#'
#' @section Notation:
#' \eqn{\nu} is the degrees of freedom and \eqn{\psi'} the trigamma function.
#'
#' @seealso [mvt_T()], which calls this on an odd dimension, and [mvt_A1()] for
#'   the digamma twin.
#'
#' @examples
#' # No step across the crossover at nu = 100.
#' signif(distributions7:::mvt_T1(c(99, 100, 101)), 8)
#'
#' # Scaled by nu^3 the series approaches -1.
#' nu <- c(1e2, 1e3, 1e4)
#' distributions7:::mvt_T1(nu) * nu^3
#'
#' @keywords internal
mvt_T1 <- function(nu) {
  # [psi'((nu+1)/2) - psi'(nu/2)]/2 + 1/nu^2 = -1/nu^3 + 1/nu^5 - 3/nu^7 + ...
  # from the same duplication student_t.cpp's t_S() uses, halved and with the
  # 1/nu^2 folded in
  out <- numeric(length(nu))
  big <- nu >= 100
  if (any(big)) {
    u <- 1 / nu[big]; u2 <- u * u
    out[big] <- u2 * u * (-1 + u2 * (1 - u2 * 3))
  }
  if (any(!big)) {
    v <- nu[!big]
    out[!big] <- 0.5 * (trigamma((v + 1) / 2) - trigamma(v / 2)) + 1 / v^2
  }
  out
}


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
#' degrees of freedom, so it is small exactly where the family is close to the
#' gaussian or the observation is close to the location, so it is the ordinary
#' case and not an edge one. Taken directly at \eqn{u = 10^{-8}} the difference
#' loses about nine significant figures.
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
#' @seealso [mvt_A()] and [mvt_T()] for the other two cancellations, and
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

#' @title Multivariate Student t Observed Hessian
#' @name distrib_hessian.MvStudentTDistrib
#'
#' @description
#' Computes the second derivatives of the log-density in closed form, by
#' differentiating the score of
#' [distrib_gradient.MvStudentTDistrib()] once more. Every block picks up a
#' term in \eqn{\partial c/\partial\cdot}, the weight
#' \eqn{c = (\nu+p)/(\nu+q)} depending on the observation through \eqn{q}. That
#' dependence is what separates the family from the gaussian, where \eqn{c} is
#' one and those terms are absent. Writing \eqn{w = \Sigma^{-1}(y-\mu)} and
#' \eqn{d = \nu + q}, the location block is
#' \deqn{\ell^{(\mu_a\mu_b)} = \frac{2c\,w_a w_b}{d} - c\,(\Sigma^{-1})_{ab},}
#' and \eqn{\partial c/\partial\nu = (q-p)/d^2} is what carries \eqn{\nu} into
#' every mixed block.
#'
#' @details
#' The pure \eqn{\nu} component has three leading terms that cancel as
#' \eqn{\nu} grows. They collapse exactly into [mvt_T()], with no series and no
#' crossover for an even dimension; see [mvt_A()] for what the direct form
#' costs.
#'
#' The quantity \eqn{w^\top A_k w} appears in three of the blocks and is formed
#' once per free value, as is \eqn{\Sigma^{-1}A_k}.
#'
#' @param distrib An [MvStudentTDistrib] object, from
#'   [mvstudent_t_distrib()].
#' @param y An \eqn{n \times p} numeric matrix of observations. A vector of
#'   length \eqn{p} is read as a single observation.
#' @param theta A named list of parameters, each component a single number.
#' @param scale One of `"parameter"` (the default) or `"link"`, handled by the
#'   generic before dispatch. The two differ here, `nu` carrying a log link by
#'   default.
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
#' parametrization, \eqn{A_k} and \eqn{A_{kl}} its first and second derivative
#' arrays, \eqn{q} the squared Mahalanobis distance, \eqn{c} the weight and
#' \eqn{\ell^{(ij)}} the second derivative of the log-density in parameters
#' \eqn{i} and \eqn{j}.
#'
#' @seealso [distrib_expected_hessian.MvStudentTDistrib()] for the expectation
#'   of this matrix, which is closed form,
#'   [distrib_gradient.MvStudentTDistrib()] for the score, [mvt_T()] for the
#'   cancellation in the \eqn{\nu} component, and [distrib_hessian()] for the
#'   generic.
#'
#' @examples
#' d <- mvstudent_t_distrib(2)
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

  d2ld <- parameters7::param_d2logdet(pc$s, pc$eta)
  spair <- param_pair_lookup(pc$s)
  sa <- lapply(pc$a, function(ak) si %*% ak)
  # w' A_k w at every observation, once: it is the derivative of q in the
  # matrix directions and appears in three of the blocks.
  wak <- lapply(pc$a, function(ak) rowSums((z$w %*% ak) * z$w))
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
      cw * wak[[k]] * z$w[, a] / den - cw * (z$w %*% t(sa[[k]]))[, a]
    } else if (a <= p && b == n_par) {
      z$w[, a] * dc_dnu
    } else if (a < n_par && b < n_par) {
      k <- a - p
      l <- b - p
      idx <- spair[[paste(min(k, l), max(k, l), sep = ":")]]
      -0.5 * d2ld[[idx]] +
        cw * wak[[l]] * wak[[k]] / (2 * den) +
        0.5 * cw * (rowSums((z$w %*% pc$a2[[idx]]) * z$w) -
          2 * rowSums((z$w %*% pc$a[[l]] %*% si %*% pc$a[[k]]) * z$w))
    } else if (a < n_par && b == n_par) {
      0.5 * wak[[a - p]] * dc_dnu
    } else {
      # d2l/dnu2
      # the first three terms collapse into mvt_T(), exactly and with no
      # series: see mvt_A()
      0.5 * (
        mvt_T(nu, p) +
          z$q / (nu * den) -
          z$q * (nu^2 + 2 * p * nu + p * z$q) / (nu * den)^2
      )
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
#' polygamma, so no sampling and no quadrature runs, and two calls at the same
#' parameters return the same numbers. No component depends on the data, so
#' every returned vector is constant across rows and `y` is read for its row
#' count alone.
#'
#' @details
#' # Why the expectations separate
#'
#' Writing \eqn{q = z^\top\Sigma^{-1}z} and \eqn{c = (\nu+p)/(\nu+q)}, the
#' variable \eqn{v = q/(q+\nu)} is EXACTLY \eqn{\mathrm{Beta}(p/2, \nu/2)} and
#' is independent of the direction \eqn{z/\lVert z\rVert}, which is uniform on
#' the sphere. Each expectation therefore factors into a radial part and an
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
#' Each block reduces to the gaussian's as \eqn{\nu \to \infty}:
#' \eqn{I_{\mu\mu} \to \Sigma^{-1}}, \eqn{I_{kl} \to T_{kl}/2}, and
#' \eqn{I_{\nu\nu} \to 0}, the degrees of freedom ceasing to be identified in
#' the limit. Measured at \eqn{p = 2} on an unstructured matrix, \eqn{I_{\nu\nu}}
#' runs \eqn{2.8\times10^{-3}}, \eqn{8.2\times10^{-6}},
#' \eqn{4.9\times10^{-9}}, \eqn{1.3\times10^{-14}} at
#' \eqn{\nu = 6, 30, 200, 5000}, which is why a standard error for \eqn{\nu} is
#' worth reading only where the tail is genuinely heavy.
#'
#' @param distrib An [MvStudentTDistrib] object, from
#'   [mvstudent_t_distrib()].
#' @param y An \eqn{n \times p} numeric matrix of observations. Only its row
#'   count is used: the expectation is taken over the law, so no observation
#'   enters any component.
#' @param theta A named list of parameters, each component a single number.
#' @param scale One of `"parameter"` (the default) or `"link"`, handled by the
#'   generic before dispatch. The two differ here, `nu` carrying a log link by
#'   default.
#' @param approx Ignored: the expectation is exact and no approximation
#'   strategy is consulted. [fit_distrib()] rejects a `fisher_scoring(approx =)`
#'   for this family for that reason.
#' @param nsim Ignored, for the same reason.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A named list of numeric vectors of length \eqn{n}, keyed as
#'   [`hess_names(distrib@params)`][hess_names], each vector constant. The
#'   returned quantities are the expected SECOND DERIVATIVES, so the
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
#' d <- mvstudent_t_distrib(2)
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
    approx = c("bartlett", "integrate", "mc", "opg"), nsim = 10000, ...) {
  y <- as_mv_matrix(distrib, y)
  pc <- mvt_pieces(distrib, theta, derivs = TRUE)
  n <- nrow(y)
  p <- pc$p
  nu <- pc$nu
  si <- pc$sigma_inv
  sa <- lapply(pc$a, function(ak) si %*% ak)
  tk <- vapply(sa, function(m) sum(diag(m)), numeric(1))
  npr <- nu + p

  i_nunu <- (trigamma(nu / 2) - trigamma(npr / 2)) / 4 -
    p / (nu * npr) + p / (2 * nu * (npr + 2))

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
#' \deqn{\frac{\partial \ell}{\partial y} = -c\,\Sigma^{-1}(y-\mu),
#'   \qquad c = \frac{\nu+p}{\nu+q},}
#' one row per observation: the gaussian's expression with the family's weight
#' in front of it. The weight is what bounds the influence of a distant row.
#' A gaussian's response gradient grows without limit as the observation moves
#' away; this one rises, turns and decays like \eqn{1/\lVert y\rVert}, which is
#' the redescending score of a heavy-tailed family.
#'
#' It is also minus the score in the location, this being a location family, so
#' the same numbers appear in [distrib_gradient.MvStudentTDistrib()].
#'
#' @param distrib An [MvStudentTDistrib] object, from
#'   [mvstudent_t_distrib()].
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
#' d <- mvstudent_t_distrib(2)
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
#' g0 <- mvgaussian_distrib(2)
#' rbind(t = distrib_grad_y(d, far, theta)[, 1],
#'       gaussian = distrib_grad_y(g0, far, theta[1:5])[, 1])
#'
#' @keywords internal
S7::method(distrib_grad_y, MvStudentTDistrib) <- function(distrib, y, theta, ...) {
  y <- as_mv_matrix(distrib, y)
  pc <- mvt_pieces(distrib, theta)
  z <- mvt_weights(y, pc)
  -z$cw * z$w
}


#' @title The Pieces Every Mixed Response Derivative of the t Is Written In
#'
#' @description
#' Assembles, once per call, the first and second derivatives in the parameters
#' of the four quantities every mixed response derivative of this family is
#' built from: \eqn{s = \nu + q}, the whitened residual \eqn{w}, the inverse
#' scale matrix \eqn{\Sigma^{-1}}, and the two scalars \eqn{c = (\nu+p)/s} and
#' \eqn{d = (\nu+p)/s^2} that follow from \eqn{s}. The response derivatives are
#' \eqn{\ell^{(y)} = -c\,w} and
#' \eqn{\ell^{(yy)} = -c\,\Sigma^{-1} + 2d\,ww^\top}, so differentiating either
#' in the parameters is a product rule over these pieces and nothing else.
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
#' One middle matrix \eqn{B_{kl}} serves all three, which keeps the assembly
#' short. \eqn{c} and \eqn{d} then follow by the quotient rule, \eqn{\nu+p}
#' being linear in \eqn{\nu} and constant in everything else.
#'
#' @param distrib An [MvStudentTDistrib] object, from
#'   [mvstudent_t_distrib()].
#' @param y An \eqn{n \times p} numeric matrix of observations, already coerced
#'   by [as_mv_matrix()].
#' @param theta A named list of parameters, each component a single number.
#'
#' @return A named list of the quantities above, INDEXED BY PARAMETER POSITION:
#'   the \eqn{p} locations, then the matrix parametrization's free values, then
#'   \eqn{\nu}. Its components are `w`, `si`, `cw`, `dw` and `q` (the values),
#'   `sa`, `wa`, `Sa`, `ca`, `da` (lists of first derivatives, one per
#'   parameter), `npar`, and `pair(a, b)`, a function returning the second-order
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
#' d <- mvstudent_t_distrib(2)
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
  v <- lapply(pc$a, function(ak) (w %*% ak) %*% si)
  gq <- lapply(pc$a, function(ak) rowSums((w %*% ak) * w))
  P <- lapply(pc$a, function(ak) si %*% ak %*% si)

  sa <- wa <- Sa <- vector("list", npar)
  Na <- numeric(npar)
  for (j in seq_len(p)) {
    sa[[j]] <- -2 * w[, j]
    wa[[j]] <- matrix(-si[j, ], n, p, byrow = TRUE)
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

  ca <- lapply(seq_len(npar), function(a) Na[a] / s - N * sa[[a]] / s^2)
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
      w_ab <- matrix(P[[k]][a, ], n, p, byrow = TRUE)
    } else if (a > p && a <= p + nf && b <= p + nf) {
      ka <- a - p
      kb <- b - p
      mid <- pc$a[[ka]] %*% si %*% pc$a[[kb]] +
        pc$a[[kb]] %*% si %*% pc$a[[ka]] - mvt_a2(pc, ka, kb)
      s_ab <- rowSums((w %*% mid) * w)
      w_ab <- (w %*% mid) %*% si
      S_ab <- si %*% mid %*% si
    }
    # nu is linear in the recursion: every pair naming it has vanishing
    # second derivatives of s, w and the inverse, and only c and d move
    c_ab <- -Na[a] * sa[[b]] / s^2 - Na[b] * sa[[a]] / s^2 -
      N * s_ab / s^2 + 2 * N * sa[[a]] * sa[[b]] / s^3
    d_ab <- -2 * Na[a] * sa[[b]] / s^3 - 2 * Na[b] * sa[[a]] / s^3 -
      2 * N * s_ab / s^3 + 6 * N * sa[[a]] * sa[[b]] / s^4
    list(s = s_ab, w = w_ab, S = S_ab, c = c_ab, d = d_ab)
  }

  list(p = p, n = n, npar = npar, si = si, w = w, s = s,
       cw = N / s, dw = N / s^2, sa = sa, wa = wa, Sa = Sa,
       ca = ca, da = da, pair = pair)
}

#' The Second Derivative of the Scale Matrix, by Position
#' @description As [mvg_a2()], for the Student t's pieces.
#' @param pc The result of [mvt_pieces()] with `derivs2`.
#' @param k,l Positions among the structure's free values.
#' @return A \eqn{p \times p} numeric matrix.
#' @keywords internal
mvt_a2 <- function(pc, k, l) {
  nm <- pc$s@free_names
  ij <- sort(c(k, l))
  pc$a2[[paste(nm[ij], collapse = ":")]]
}

#' @title An Outer Product Per Observation
#'
#' @description
#' Builds the \eqn{p \times p \times n} array whose \eqn{i}th slice is
#' \eqn{a_ib_i^\top}, the outer product of row \eqn{i} of `a` with row \eqn{i}
#' of `b`. The whole array is assembled by two column replications and one
#' permutation, with no loop over the observations, so the response derivatives
#' of a multivariate family stay affordable at a sample size worth fitting.
#'
#' @param a,b Numeric matrices with the same dimensions, \eqn{n} rows and
#'   \eqn{p} columns. Nothing is validated: mismatched dimensions recycle
#'   silently, as they would in the arithmetic underneath.
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
#' Unlike the gaussian's, no component here is constant across rows: the
#' location block carries the observation through \eqn{c} even though its
#' gaussian part does not.
#'
#' @details
#' Nothing here is obstructed. The log-density carries no distribution
#' function, only `lgamma`, a logarithm and a quadratic form, each elementary
#' in \eqn{\nu}. As \eqn{\nu \to \infty} the weight goes to one and its
#' derivatives to zero, so every component becomes the gaussian's; the tests
#' compare them against that limit, and at \eqn{\nu = 10^7} the two agree to
#' \eqn{2.6\times10^{-6}}.
#'
#' @param distrib An [MvStudentTDistrib] object, from
#'   [mvstudent_t_distrib()].
#' @param y An \eqn{n \times p} numeric matrix of observations. A vector of
#'   length \eqn{p} is read as a single observation.
#' @param theta A named list of parameters, each component a single number.
#' @param scale One of `"parameter"` (the default) or `"link"`, handled by the
#'   generic before dispatch. The two differ here, `nu` carrying a log link by
#'   default.
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
#' d <- mvstudent_t_distrib(2)
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
#' g <- mvgaussian_distrib(2)
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
#' The shape differs from the gaussian's, and the difference is the family's
#' defining property. There the response Hessian is \eqn{-\Sigma^{-1}}, one
#' matrix for the whole sample, so its derivative is one matrix per parameter.
#' Here it depends on the observation through the weight, so this method and
#' [distrib_hess_y_hess.MvStudentTDistrib()] return ONE MATRIX PER ROW. Nor
#' does any component vanish: the gaussian's location components are exactly
#' zero and none of these is.
#'
#' This is one of the three derivatives a marginal criterion reads when the
#' family stands as a prior over a coefficient block.
#'
#' @param distrib An [MvStudentTDistrib] object, from
#'   [mvstudent_t_distrib()].
#' @param y An \eqn{n \times p} numeric matrix of observations. A vector of
#'   length \eqn{p} is read as a single observation.
#' @param theta A named list of parameters, each component a single number.
#' @param scale One of `"parameter"` (the default) or `"link"`, handled by the
#'   generic before dispatch. The two differ here.
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
#' d <- mvstudent_t_distrib(2)
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
#' distrib_cross2_y(mvgaussian_distrib(2), y, theta[1:5])$mu1
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
#' with every piece from [mvt_dpieces()]. No pair vanishes here, where a
#' gaussian's pair of location parameters is exactly zero: its response
#' gradient is linear in the location and this one is not.
#'
#' @param distrib An [MvStudentTDistrib] object, from
#'   [mvstudent_t_distrib()].
#' @param y An \eqn{n \times p} numeric matrix of observations. A vector of
#'   length \eqn{p} is read as a single observation.
#' @param theta A named list of parameters, each component a single number.
#' @param scale One of `"parameter"` (the default) or `"link"`, handled by the
#'   generic before dispatch. The two differ here.
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
#' d <- mvstudent_t_distrib(2)
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
#' distrib_grad_y_hess(mvgaussian_distrib(2), y, theta[1:5])$mu1_mu2
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
#' [mvt_dpieces()]. No pair vanishes, where for a gaussian every pair naming a
#' location parameter is exactly zero.
#'
#' @param distrib An [MvStudentTDistrib] object, from
#'   [mvstudent_t_distrib()].
#' @param y An \eqn{n \times p} numeric matrix of observations. A vector of
#'   length \eqn{p} is read as a single observation.
#' @param theta A named list of parameters, each component a single number.
#' @param scale One of `"parameter"` (the default) or `"link"`, handled by the
#'   generic before dispatch. The two differ here.
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
#' d <- mvstudent_t_distrib(2)
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
#' one \eqn{p \times p} matrix PER OBSERVATION. The gaussian's response
#' Hessian is one matrix for the whole sample; this one moves with the
#' observation, and far enough out its rank-one term makes the curvature
#' positive along the direction of the residual, which is the redescending
#' score seen at second order.
#'
#' @param distrib An [MvStudentTDistrib] object, from
#'   [mvstudent_t_distrib()].
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
#' d <- mvstudent_t_distrib(2)
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
#' distrib_hess_y(mvgaussian_distrib(2), far, theta[1:5])[1, 1]
#'
#' @keywords internal
S7::method(distrib_hess_y, MvStudentTDistrib) <- function(distrib, y, theta, ...) {
  y <- as_mv_matrix(distrib, y)
  pc <- mvt_pieces(distrib, theta)
  z <- mvt_weights(y, pc)
  n <- nrow(y)
  p <- pc$p
  out <- array(0, dim = c(p, p, n))
  for (i in seq_len(n)) {
    wi <- z$w[i, ]
    out[, , i] <- -z$cw[i] * pc$sigma_inv +
      (2 * z$cw[i] / (pc$nu + z$q[i])) * tcrossprod(wi)
  }
  out
}


#' @title Mean of a Multivariate Student t
#' @name mean.MvStudentTDistrib
#'
#' @description
#' Returns \eqn{\mathbb{E}[Y] = \mu} for \eqn{\nu > 1}, and a vector of `NaN`
#' otherwise. Below one degree of freedom the defining integral does not
#' converge: the density decays like \eqn{\lVert y\rVert^{-(\nu+p)}} and the
#' first absolute moment integrates \eqn{\lVert y\rVert^{1-\nu-p}} over a shell
#' of surface \eqn{\lVert y\rVert^{p-1}}, which is finite only for
#' \eqn{\nu > 1}. The location is a parameter and the mean is a moment; at
#' \eqn{\nu \le 1} the first exists and the second does not, so `NaN` is the
#' answer and \eqn{\mu} would be the wrong one.
#'
#' @param x An [MvStudentTDistrib] object, from [mvstudent_t_distrib()].
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
#' d <- mvstudent_t_distrib(2)
#' theta <- list(mu1 = 0.5, mu2 = -0.3, sigma_log_L1 = 0.1,
#'               sigma_log_L2 = -0.2, sigma_L2.1 = 0.4, nu = 6)
#'
#' mean(d, theta)
#'
#' # Below one degree of freedom the mean does not exist, while the location
#' # and the density both do.
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
#' @param x An [MvStudentTDistrib] object, from [mvstudent_t_distrib()].
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
#' d <- mvstudent_t_distrib(2)
#' theta <- list(mu1 = 0.5, mu2 = -0.3, sigma_log_L1 = 0.1,
#'               sigma_log_L2 = -0.2, sigma_L2.1 = 0.4, nu = 6)
#'
#' variance(d, theta)
#' mv_sigma(d, theta)
#' all.equal(variance(d, theta), (6 / 4) * mv_sigma(d, theta))
#'
#' # The factor blows up as nu falls to 2 and the covariance goes with it.
#' vapply(c(10, 3, 2.5, 2.05, 2), function(nu) {
#'   t2 <- theta; t2$nu <- nu
#'   variance(d, t2)[1, 1]
#' }, numeric(1))
#'
#' # Correlations are unaffected: a positive multiple leaves them alone.
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
#' Returns the marginal law of a subset of coordinates, which for this family
#' is again a Student t: the subvector of the location, the corresponding block
#' of the scale matrix, and THE SAME degrees of freedom. \eqn{\nu} does not
#' change with the dimension, and that is why the family is closed under
#' marginalization: the mixing variable of the scale-mixture representation is
#' shared by every coordinate, so it survives integrating any of them out.
#'
#' The returned object is a fresh [mvstudent_t_distrib()] of the reduced
#' dimension on an unstructured scale matrix, so the free values are recomputed
#' from the block by `parameters7::param_free()` and are not a subset of the
#' original's.
#'
#' @param distrib An [MvStudentTDistrib] object, from
#'   [mvstudent_t_distrib()].
#' @param theta A named list of parameters, each component a single number.
#' @param which An integer vector of coordinates to keep, between 1 and
#'   \eqn{p}. Duplicates and out-of-range values are not checked and reach the
#'   matrix subsetting.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A named list with `distrib`, an `MvStudentTDistrib` of dimension
#'   `length(which)`, and `theta`, its parameters as a named list.
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
#' d <- mvstudent_t_distrib(3)
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
  md <- mvstudent_t_distrib(length(which))
  eta <- parameters7::param_free(md@param, unname(sg))
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
#'   [mvstudent_t_distrib()].
#' @param theta A named list of parameters, each component a single number.
#'   The location components and `nu` are ignored.
#'
#' @return A \eqn{p \times p} symmetric positive definite numeric matrix, with
#'   both dimnames `v1`, ..., `vp`.
#'
#' @seealso [variance.MvStudentTDistrib()] for the covariance,
#'   [mv_sigma.MvGaussianDistrib()], where the two coincide,
#'   [mv_summary()], which reports the square roots of this matrix's diagonal
#'   as `scale_sd_v1` and so on, and [mv_sigma()] for the generic.
#'
#' @examples
#' d <- mvstudent_t_distrib(2)
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
#' Draws a parameter vector for testing: each location uniform on
#' \eqn{(-1, 1)}, each free value of the matrix parametrization uniform on
#' \eqn{(-0.4, 0.4)}, and the degrees of freedom uniform on \eqn{(3, 12)}. The
#' matrix band puts the scale matrix near the identity, as it does for the
#' gaussian; the \eqn{\nu} band is chosen so that the family is genuinely
#' heavy-tailed and every moment [check_distrib()] compares still exists, the
#' variance needing \eqn{\nu > 2} and the kurtosis \eqn{\nu > 4}.
#'
#' @param distrib An [MvStudentTDistrib] object, from
#'   [mvstudent_t_distrib()].
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
#' d <- mvstudent_t_distrib(2)
#'
#' set.seed(11)
#' unlist(generate_random_theta(d))
#'
#' # The three bands, over 400 draws.
#' set.seed(12)
#' round(apply(replicate(400, unlist(generate_random_theta(d))), 1, range), 2)
#'
#' # Every draw keeps nu above 4, so the covariance and the kurtosis of every
#' # coordinate exist.
#' set.seed(13)
#' all(replicate(50, is.finite(variance(d, generate_random_theta(d))[1, 1])))
#'
#' @keywords internal
S7::method(generate_random_theta, MvStudentTDistrib) <- function(distrib, ...) {
  p <- distrib@n_dim
  s <- distrib@param
  as.list(stats::setNames(
    c(stats::runif(p, -1, 1), stats::runif(s@n_free, -0.4, 0.4),
      stats::runif(1, 3, 12)),
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
#' admissible \eqn{\nu} and the mean only for \eqn{\nu > 1}, which is why the
#' generic is named for a location: [mean.MvStudentTDistrib()] returns `NaN`
#' below one degree of freedom and this returns the location regardless.
#'
#' @param distrib An [MvStudentTDistrib] object, from
#'   [mvstudent_t_distrib()].
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
#' d <- mvstudent_t_distrib(2)
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
