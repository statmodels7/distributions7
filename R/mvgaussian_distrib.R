#' @include multivariate.R cross_derivatives.R cross2_derivatives.R cross_theta2_derivatives.R
NULL

#' @title Multivariate Gaussian Distribution Class
#' @name MvGaussianDistrib
#'
#' @description
#' The S7 class of the gaussian family on \eqn{\mathbb{R}^p}, with density
#' \deqn{f(y) = (2\pi)^{-p/2}\lvert\Sigma\rvert^{-1/2}
#' \exp\!\left\{-\tfrac{1}{2}(y-\mu)^\top\Sigma^{-1}(y-\mu)\right\}.} The mean
#' \eqn{\mu} contributes \eqn{p} scalar parameters and the matrix is carried by
#' a \pkg{parameters7} parametrization, whose free values become scalar
#' parameters in their turn. It inherits from `multivariate_distrib`, so the
#' response is an \eqn{n \times p} matrix and the distribution function and the
#' quantile function signal an error.
#'
#' Build one with [mvgaussian1_distrib()], which fills the properties in and
#' checks that the matrix parametrization has full rank. This page documents
#' the raw S7 constructor, which validates neither the rank nor the agreement
#' between `n_dim` and the parametrization's dimension.
#'
#' @param param A \pkg{parameters7} parametrization of the matrix, inheriting
#'   from `parameters7::parameter`. Its `n_free` free values are flattened
#'   into scalar parameters of the distribution.
#' @param inverted Logical of length 1. `TRUE` when `param` carries the
#'   precision \eqn{\Omega = \Sigma^{-1}} and `FALSE` when it carries the
#'   covariance \eqn{\Sigma}. It decides which derivative arrays are the
#'   parametrization's own and which are obtained from them (see
#'   [mv_matrix_pieces()]), the sign of the log-determinant, the prefix of the
#'   free names, and the quantities that [mv_summary()] reports. The same matrix
#'   \eqn{\Sigma} gives the same law either way.
#' @param n_dim The dimension \eqn{p} of one observation. A single positive
#'   integer.
#' @inheritParams distrib
#'
#' @return An S7 object of class `MvGaussianDistrib`, inheriting from
#'   `multivariate_distrib` and from `distrib`. Beyond the parent's
#'   `distrib_name`, `dimension`, `bounds`, `params`,
#'   `params_interpretation`, `n_params`, `params_bounds`, `link_params`,
#'   `params_smooth` and `n_dim`, it carries `param` and `inverted` as
#'   supplied. For an object built by [mvgaussian1_distrib()] at \eqn{p = 2} on
#'   an unstructured covariance, `params` is
#'   `c("mu1", "mu2", "sigma_log_L1", "sigma_log_L2", "sigma_L2.1")`, every
#'   `params_bounds` entry is \eqn{(-\infty, \infty)} and every link is the
#'   identity.
#'
#' @seealso [mvgaussian1_distrib()] to build one, [mvstudent_t1_distrib()] for
#'   the heavy-tailed sibling, [mv_sigma()] for the covariance a parameter
#'   vector describes, and [distrib_gradient.MvGaussianDistrib()] for the
#'   closed-form score.
#'
#' @section Methods:
#' Registered on this class, each linked to the page that documents it:
#' - the distribution: [`distrib_pdf()`][distrib_pdf.MvGaussianDistrib] and
#'   [`distrib_rng()`][distrib_rng.MvGaussianDistrib];
#' - the derivatives in the parameters:
#'   [`distrib_gradient()`][distrib_gradient.MvGaussianDistrib],
#'   [`distrib_hessian()`][distrib_hessian.MvGaussianDistrib],
#'   [`distrib_expected_hessian()`][distrib_expected_hessian.MvGaussianDistrib],
#'   [`distrib_deriv3()`][distrib_deriv3.MvGaussianDistrib] and
#'   [`distrib_deriv4()`][distrib_deriv4.MvGaussianDistrib];
#' - the derivatives in the response:
#'   [`distrib_grad_y()`][distrib_grad_y.MvGaussianDistrib] and
#'   [`distrib_hess_y()`][distrib_hess_y.MvGaussianDistrib];
#' - the mixed derivatives in the response and the parameters:
#'   [`distrib_cross_y()`][distrib_cross_y.MvGaussianDistrib],
#'   [`distrib_cross2_y()`][distrib_cross2_y.MvGaussianDistrib],
#'   [`distrib_grad_y_hess()`][distrib_grad_y_hess.MvGaussianDistrib] and
#'   [`distrib_hess_y_hess()`][distrib_hess_y_hess.MvGaussianDistrib];
#' - the moments: [`mean()`][mean.MvGaussianDistrib] and
#'   [`variance()`][variance.MvGaussianDistrib];
#' - the multivariate summaries:
#'   [`mv_location()`][mv_location.MvGaussianDistrib],
#'   [`mv_sigma()`][mv_sigma.MvGaussianDistrib],
#'   [`mv_marginal()`][mv_marginal.MvGaussianDistrib] and
#'   [`mv_derived()`][mv_derived.MvGaussianDistrib];
#' - the validation and fitting helpers:
#'   [`generate_random_theta()`][generate_random_theta.MvGaussianDistrib] and
#'   [`distrib_start()`][distrib_start.MvGaussianDistrib].
#'
#' Every other generic is inherited from [multivariate_distrib()], whose methods
#' signal an error for the distribution function, the quantile function, the
#' skewness and the kurtosis.
#'
#' @examples
#' d <- mvgaussian1_distrib(2)
#' S7::S7_inherits(d, multivariate_distrib)
#'
#' # Five scalar parameters: two means and the three free values of the
#' # log-Cholesky covariance, prefixed by the matrix that they describe.
#' d@params
#' d@n_dim
#' d@param@free_names
#'
#' # Every parameter is already unconstrained, so every link is the identity
#' # and the positive definiteness lives inside the matrix parametrization.
#' vapply(d@link_params, function(l) l@link_name, character(1))
#'
#' # The precision form differs in one property and in the prefix.
#' o <- mvgaussian2_distrib(2, parameters7::log_cholesky(2))
#' c(covariance = d@inverted, precision = o@inverted)
#' o@params
#'
#' @export
MvGaussianDistrib <- S7::new_class("MvGaussianDistrib",
  parent = multivariate_distrib,
  properties = list(
    param = parameters7::parameter,
    inverted = S7::class_logical
  )
)


#' @title S7 Classes for the Two Multivariate Gaussian Parametrizations
#'
#' @description
#' `MvGaussian1Distrib` is the family whose matrix parametrization is the
#' covariance and `MvGaussian2Distrib` the family whose matrix
#' parametrization is the precision. Both inherit from [MvGaussianDistrib],
#' which carries every method; the two differ in what their free values
#' describe, and the computations are shared.
#'
#' @details
#' The two are separate families for the reason that [gaussian1_distrib()],
#' [gaussian2_distrib()] and [gaussian3_distrib()] are separate families: a
#' parametrization is a modeling decision, and a predictor imposed on the
#' coordinates of \eqn{\Sigma} is a different model from the same predictor
#' imposed on the coordinates of \eqn{\Omega}. Where the matrix parametrization
#' is closed under inversion the two describe the same laws in different
#' coordinates; where it is not, they describe different laws.
#'
#' @inheritParams MvGaussianDistrib
#'
#' @return An object of class `MvGaussian1Distrib` or `MvGaussian2Distrib`, a
#'   subclass of [MvGaussianDistrib] adding no properties of its own.
#'
#' @seealso [mvgaussian1_distrib()] and [mvgaussian2_distrib()], the
#'   constructors, and [MvGaussianDistrib] for the properties and the methods.
#'
#' @examples
#' c(covariance = S7::S7_inherits(mvgaussian1_distrib(2), MvGaussian1Distrib),
#'   precision = S7::S7_inherits(mvgaussian2_distrib(2), MvGaussian2Distrib))
#'
#' # Both are multivariate gaussians, the class that carries the methods.
#' S7::S7_inherits(mvgaussian2_distrib(2), MvGaussianDistrib)
#'
#' @export
MvGaussian1Distrib <- S7::new_class("MvGaussian1Distrib",
  parent = MvGaussianDistrib
)

#' @rdname MvGaussian1Distrib
#' @export
MvGaussian2Distrib <- S7::new_class("MvGaussian2Distrib",
  parent = MvGaussianDistrib
)


#' @title Construct a Multivariate Gaussian Distribution
#'
#' @description
#' Two families on \eqn{\mathbb{R}^p}, differing in which matrix of the law
#' their parametrization describes: `mvgaussian1_distrib()` parametrizes the
#' covariance \eqn{\Sigma} and `mvgaussian2_distrib()` the precision
#' \eqn{\Omega = \Sigma^{-1}}. The mean is a vector of \eqn{p} free parameters
#' and the matrix is carried by a \pkg{parameters7} parametrization whose free
#' values become ordinary scalar parameters, so every generic of the `distrib`
#' class applies to the object, with `theta` a flat named list of numbers. Both
#' default to an unstructured matrix in the log-Cholesky parametrization, which
#' has `p * (p + 1) / 2` free values.
#'
#' @details
#' # The two families
#'
#' The numbering is the package's convention, the one that separates
#' [gaussian1_distrib()], [gaussian2_distrib()] and [gaussian3_distrib()]: one
#' name per parametrization, because in a modeling framework a parametrization
#' decides what a linear predictor acts on. Here the choice can change the model
#' as well as the coordinates. A matrix parametrization is closed under
#' inversion when the inverse of every matrix that it produces is one of its
#' own; [parameters7::log_cholesky()], [parameters7::matrix_log()],
#' [parameters7::diagonal_matrix()], [parameters7::scalar_matrix()],
#' [parameters7::compound_symmetry()] and [parameters7::dr_prod()] are, and
#' [parameters7::correlation_matrix()], [parameters7::ar1()] and
#' [parameters7::autoregressive()] are not.
#'
#' Where it is closed the two families describe the same laws in different
#' coordinates, and a fit reaches the same maximum. Where it is not closed they
#' are different models: the inverse of an AR(1) covariance is tridiagonal and
#' is not an AR(1) matrix at any parameters, so [parameters7::ar1()] on the
#' covariance and on the precision give different maximized log-likelihoods.
#'
#' Both constructors accept any parametrization of the right dimension and of
#' full rank. `mvgaussian2_distrib(p, parameters7::ar1(p))` is the model whose
#' precision has the entries \eqn{\sigma^2\rho^{|i-j|}}, with \eqn{\sigma^2}
#' the scale of `ar1()`, and it is a valid law. The AR(1) process written
#' on the precision side is `mvgaussian2_distrib(p, parameters7::ar1_inv(p))`,
#' whose parametrization is that inverse.
#'
#' # The arithmetic of the two sides
#'
#' The derivatives of both families are computed from the derivative arrays of
#' \eqn{\Sigma} and of \eqn{\Sigma^{-1}}. The parametrization supplies one of
#' the two, and [mv_matrix_pieces()] obtains the other: through
#' [parameters7::param_inv_d1()] and its siblings under a covariance
#' parametrization, and by the chain rule for an inverse under a precision one.
#'
#' # Parameters
#'
#' The mean contributes `mu1`, ..., `mup`. The matrix contributes the
#' parametrization's free values, prefixed by the matrix that they describe:
#' `sigma_` for a covariance and `omega_` for a precision. A two-dimensional
#' gaussian on an unstructured covariance therefore has the five parameters
#' `mu1`, `mu2`, `sigma_log_L1`, `sigma_log_L2` and `sigma_L2.1`, while the
#' same parametrization on the precision gives `omega_log_L1` and the rest. The
#' prefix distinguishes the two models in a printed table, since the name of a
#' free value gives the construction of the matrix and does not say which
#' matrix it is.
#'
#' All of these parameters are unconstrained, so all of the links are the
#' identity. The constraint that makes the matrix positive definite is carried
#' by the matrix parametrization, which needs no link to express it. The
#' parameter scale and the link scale therefore coincide, and
#' `scale = "link"` returns the same numbers as `scale = "parameter"`.
#'
#' # Reading a fit
#'
#' The free values are coordinates of the parametrization and have no direct
#' interpretation. [mv_summary()] carries the variance matrix of a fit onto the
#' standard deviations and the correlations by the delta method, and `print()`
#' shows those. A precision parametrization reports instead the conditional
#' variances and the partial correlations, which it describes directly:
#' \eqn{1/\Omega_{jj} = \operatorname{Var}(Y_j \mid Y_{-j})}, whose ratio to
#' the marginal variance is \eqn{1 - R_j^2} for the regression of that
#' coordinate on all the others.
#'
#' # Rank
#'
#' A rank-deficient parametrization is rejected with an error. A singular
#' covariance gives a law supported on a subspace, with no density against
#' Lebesgue measure; a singular precision gives a quadratic form that is flat
#' along its null space and does not normalize. Such a parametrization can
#' serve as a penalty but does not define a density.
#'
#' # The response
#'
#' An \eqn{n \times p} matrix, one row per observation. A plain vector of
#' length \eqn{p} is read as a single observation. A parameter may not vary from
#' observation to observation here, the parametrization describing one matrix
#' for the whole sample.
#'
#' @section The distribution:
#' \deqn{f(y) = (2\pi)^{-p/2}\lvert \Sigma \rvert^{-1/2}\exp\!\left\{-\tfrac{1}{2}(y-\mu)^\top\Sigma^{-1}(y-\mu)\right\}}
#' on \eqn{y \in \mathbb{R}^{p}}, with
#' \deqn{\mathbb{E}[Y] = \mu, \qquad \operatorname{Var}(Y) = \Sigma.}
#'
#' @section Notation:
#' \eqn{\mu} is the mean vector, \eqn{\Sigma} the covariance, \eqn{\Omega} the
#' precision, \eqn{p} the dimension of one observation and \eqn{n} the number
#' of observations. \eqn{\eta} is the free vector of the matrix
#' parametrization, the unconstrained coordinates that an optimizer moves.
#'
#' @param n_dim The dimension \eqn{p} of one observation. A single positive
#'   whole number, finite and at least 1. Any other value signals an error.
#' @param sigma A \pkg{parameters7} parametrization of the covariance, of
#'   dimension `n_dim` and of full rank. Defaults to
#'   `parameters7::log_cholesky(n_dim)`.
#' @param omega A \pkg{parameters7} parametrization of the precision, of
#'   dimension `n_dim` and of full rank. Defaults to
#'   `parameters7::log_cholesky(n_dim)`.
#'
#' @return An S7 object of class [MvGaussian1Distrib] or [MvGaussian2Distrib],
#'   with `param` the parametrization supplied and `inverted` recording which
#'   side it carries. Its `params` are `mu1`, ..., `mup` followed by the
#'   prefixed free names, `n_params` is `p + param@n_free`, every entry of
#'   `params_bounds` is \eqn{(-\infty, \infty)} and every link is the identity.
#'
#' @seealso [gaussian1_distrib()] for the one-dimensional case,
#'   [mvstudent_t1_distrib()] for the heavy-tailed sibling,
#'   [mv_summary()] for the quantities that a fit reports,
#'   [fit_distrib()] to estimate one,
#'   [parameters7::log_cholesky()] for the default parametrization, and
#'   [parameters7::inverse_of()] for writing a structure on the other side.
#'
#' @examples
#' d <- mvgaussian1_distrib(2)
#' d
#'
#' theta <- list(mu1 = 0.5, mu2 = -0.3, sigma_log_L1 = 0.1,
#'               sigma_log_L2 = -0.2, sigma_L2.1 = 0.4)
#'
#' # The covariance that the free values describe, positive definite for any
#' # five numbers.
#' S <- mv_sigma(d, theta)
#' round(S, 4)
#' eigen(S, only.values = TRUE)$values
#'
#' # The density, against the formula written out.
#' y <- rbind(c(0, 0), c(1, -1))
#' mu <- c(0.5, -0.3)
#' ref <- apply(y, 1, function(r)
#'   -0.5 * (2 * log(2 * pi) + log(det(S)) + t(r - mu) %*% solve(S) %*% (r - mu)))
#' all.equal(distrib_pdf(d, y, theta, log = TRUE), as.numeric(ref))
#'
#' # A structured covariance spends fewer parameters: two variances, no
#' # correlation.
#' mvgaussian1_distrib(2, parameters7::diagonal_matrix(2))@params
#'
#' # The precision form is the other family, and its free values are prefixed
#' # by the matrix that they describe.
#' o <- mvgaussian2_distrib(2)
#' o@params
#' th_o <- list(mu1 = 0, mu2 = 0, omega_log_L1 = 0.1,
#'              omega_log_L2 = -0.2, omega_L2.1 = 0.4)
#' Om <- parameters7::param_value(o@param, unlist(th_o)[3:5])
#' all.equal(mv_sigma(o, th_o), solve(Om), check.attributes = FALSE)
#'
#' # The conditional variance of the first coordinate given the second, and
#' # the marginal variance beside it.
#' c(conditional = 1 / Om[1, 1], marginal = mv_sigma(o, th_o)[1, 1])
#'
#' # An AR(1) covariance and an AR(1) precision are different models, and the
#' # name of the family records which one is meant.
#' c(covariance = mvgaussian1_distrib(4, parameters7::ar1(4))@distrib_name,
#'   precision = mvgaussian2_distrib(4, parameters7::ar1(4))@distrib_name)
#'
#' @export
mvgaussian1_distrib <- function(n_dim, sigma = NULL) {
  mvgaussian_build(n_dim, sigma, inverted = FALSE, cls = MvGaussian1Distrib)
}

#' @rdname mvgaussian1_distrib
#' @export
mvgaussian2_distrib <- function(n_dim, omega = NULL) {
  mvgaussian_build(n_dim, omega, inverted = TRUE, cls = MvGaussian2Distrib)
}


#' @title Build One of the Two Multivariate Gaussian Families
#'
#' @description
#' The body that [mvgaussian1_distrib()] and [mvgaussian2_distrib()] share. The
#' two differ in the class that they instantiate, in the `inverted` property,
#' in the prefix of their free names, in the tag of `distrib_name` (`sigma=` or
#' `omega=`) and in what `params_interpretation` calls the matrix. Every method
#' is registered on [MvGaussianDistrib], the parent of both.
#'
#' @param n_dim The dimension of one observation.
#' @param s The matrix parametrization, or `NULL` for
#'   `parameters7::log_cholesky(n_dim)`.
#' @param inverted `TRUE` where `s` carries the precision.
#' @param cls The S7 class to instantiate.
#'
#' @return An object of class `cls`.
#'
#' @seealso [mvgaussian1_distrib()], the constructors that it serves.
#'
#' @keywords internal
mvgaussian_build <- function(n_dim, s, inverted, cls) {
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
      "  density. A singular covariance is supported on a subspace and a\n",
      "  singular precision does not normalize; either way the law has no\n",
      "  density against Lebesgue measure. Such a structure is a penalty, not\n",
      "  a distribution."
    ), s@rank, s@dimension), call. = FALSE)
  }

  mu_names <- paste0("mu", seq_len(p))
  free_names <- mv_prefixed_names(s@free_names, inverted)
  clash <- intersect(mu_names, free_names)
  if (length(clash)) {
    stop(sprintf(paste0(
      "The matrix parameter's free value '%s' has the same name as a mean component.\n",
      "  Parameter names must be unique, since every derivative component is\n",
      "  keyed by them."
    ), clash[1L]), call. = FALSE)
  }

  params <- c(mu_names, free_names)
  n_par <- length(params)

  cls(
    # No spaces inside the brackets: print() capitalises the first letter after
    # every space in a distribution's name, and "Covariance Log_cholesky" is
    # not what the matrix parameter is called. Same convention as truncated().
    distrib_name = sprintf(
      "multivariate gaussian [%dd, %s=%s]", p,
      if (inverted) "omega" else "sigma", s@param_name
    ),
    dimension = "multivariate",
    n_dim = p,
    bounds = c(-Inf, Inf),
    params = params,
    params_interpretation = stats::setNames(
      c(rep("mean", p), rep(
        if (inverted) "precision" else "covariance", s@n_free
      )),
      params
    ),
    n_params = n_par,
    # Every parameter is already unconstrained: the mean is free, and the
    # parameter's free values are free by construction. The constraint that
    # makes the matrix positive definite is inside the matrix parameter, which is
    # exactly why it does not need a link to carry it.
    params_bounds = stats::setNames(
      rep(list(c(-Inf, Inf)), n_par), params
    ),
    link_params = stats::setNames(
      rep(list(linkfunctions7::identity_link()), n_par), params
    ),
    param = s,
    inverted = inverted
  )
}


#' @title The Pieces From Which a Multivariate Gaussian Evaluates
#'
#' @description
#' Assembles, once per call, the mean, the covariance, the inverse covariance
#' and the log-determinant from a flat parameter vector, together with the
#' derivative arrays of \eqn{\Sigma} and of \eqn{\Sigma^{-1}} when they are
#' requested. The density, generator, derivative and moment methods of
#' [MvGaussianDistrib] call this first and work from the result, so a parameter
#' vector is unpacked and a matrix factorized once per call and not once per
#' component.
#'
#' @details
#' The matrix quantities come from [mv_matrix_pieces()], which the Student t
#' reads as well. Whichever side the parametrization carries, the result holds
#' \eqn{\Sigma}, \eqn{\Sigma^{-1}}, the log-determinant of \eqn{\Sigma} and the
#' derivatives of both matrices.
#'
#' @param distrib An [MvGaussianDistrib] object.
#' @param theta A named list of parameters, already aligned by the generic, or
#'   any list whose components are in `distrib@params` order.
#' @param derivs Logical of length 1. When `TRUE` the first derivative arrays
#'   are computed. Defaults to `FALSE`.
#' @param derivs2 Logical of length 1. When `TRUE` the second derivative arrays
#'   are computed, and the first ones with them, so `derivs` need not also be
#'   set. Defaults to `FALSE`.
#'
#' @return A named list with `mu` (numeric of length \eqn{p}), `eta` (the
#'   matrix parametrization's free vector), `p`, `s` (the parametrization
#'   itself) and the elements of [mv_matrix_pieces()]: `sigma` and `sigma_inv`
#'   (\eqn{p \times p} matrices), `logdet` (the log-determinant of
#'   \eqn{\Sigma}) and `inverted`; with derivatives requested also `a` and `ai`
#'   (the first derivative arrays of \eqn{\Sigma} and of \eqn{\Sigma^{-1}}, one
#'   per free value) and `dlogdet`; with `derivs2` also `a2`, `ai2` (one array
#'   per unordered pair in `parameters7::param_tuple_indices()` order) and
#'   `d2logdet`.
#'
#' @section Notation:
#' \eqn{\Sigma} is the covariance, \eqn{\Omega = \Sigma^{-1}} the precision and
#' \eqn{\eta} the free vector of the matrix parametrization.
#'
#' @seealso [mvgaussian1_distrib()] for the family, [mv_matrix_pieces()] for the
#'   matrix quantities, and [distrib_gradient.MvGaussianDistrib()] for a
#'   consumer.
#'
#' @examples
#' d <- mvgaussian1_distrib(2)
#' theta <- list(mu1 = 0.5, mu2 = -0.3, sigma_log_L1 = 0.1,
#'               sigma_log_L2 = -0.2, sigma_L2.1 = 0.4)
#' pc <- distributions7:::mvg_pieces(d, theta, derivs = TRUE)
#' names(pc)
#'
#' # sigma_inv is the inverse, and logdet its log-determinant.
#' all.equal(pc$sigma %*% pc$sigma_inv, diag(2), check.attributes = FALSE)
#' all.equal(pc$logdet, log(det(pc$sigma)))
#'
#' # One derivative array per free value of the parametrization.
#' length(pc$a)
#' round(pc$a[[3]], 4)
#'
#' # The same pieces from the precision side describe the same covariance.
#' o <- mvgaussian2_distrib(2, parameters7::log_cholesky(2))
#' th_o <- list(mu1 = 0.5, mu2 = -0.3, omega_log_L1 = 0.1,
#'              omega_log_L2 = -0.2, omega_L2.1 = 0.4)
#' po <- distributions7:::mvg_pieces(o, th_o)
#' all.equal(po$logdet, -parameters7::param_logdet(o@param, unlist(th_o)[3:5]))
#'
#' @keywords internal
mvg_pieces <- function(distrib, theta, derivs = FALSE, derivs2 = FALSE) {
  p <- distrib@n_dim
  s <- distrib@param
  v <- mv_flat_theta(distrib, theta)
  mu <- v[seq_len(p)]
  eta <- v[p + seq_len(s@n_free)]

  # the matrix, its inverse and their arrays, transported once where the
  # parametrization carries the precision; written in mv_matrix_pieces() so
  # that the Student t reads the same arithmetic and not a second copy of it
  c(list(mu = unname(mu), eta = eta, p = p, s = s),
    mv_matrix_pieces(s, eta, distrib@inverted, derivs, derivs2))
}


#' @title Mean Vector of a Multivariate Gaussian
#' @name mv_location.MvGaussianDistrib
#'
#' @description
#' Returns the mean vector \eqn{\mu}, which for this family is the first \eqn{p}
#' parameters read off `theta` in order. The method is [mv_leading_location()],
#' shared with the Student t, whose location is also its leading parameters. For
#' the gaussian \eqn{\mu} is both the center of the density and the expectation
#' of the law, which does not hold for every family.
#'
#' @param distrib An [MvGaussianDistrib] object, from [mvgaussian1_distrib()].
#' @param theta A named list of parameters, one number each. Only the \eqn{p}
#'   mean components are read.
#'
#' @return A numeric vector of length \eqn{p}, named `v1`, ..., `vp` after the
#'   coordinates of the response. The parameter names `mu1`, ..., `mup` do not
#'   appear on it.
#'
#' @seealso [mean.MvGaussianDistrib()], which returns the same vector as the
#'   expectation of the law, [mv_sigma.MvGaussianDistrib()] for the matrix,
#'   and [mv_location()] for the generic.
#'
#' @examples
#' d <- mvgaussian1_distrib(3)
#' theta <- as.list(stats::setNames(c(1, -2, 0.5, rep(0, 6)), d@params))
#'
#' mv_location(d, theta)
#'
#' # For a gaussian the location is also the expectation, so the two agree.
#' all.equal(unname(mv_location(d, theta)), unname(mean(d, theta)))
#'
#' @keywords internal
S7::method(mv_location, MvGaussianDistrib) <- mv_leading_location

#' @title Covariance of a Multivariate Gaussian
#' @name mv_sigma.MvGaussianDistrib
#'
#' @description
#' Returns the covariance \eqn{\Sigma}, assembled from the matrix
#' parametrization's free values. Where the parametrization carries the
#' precision the matrix is inverted first, so the result is the covariance
#' either way. For this family the matrix is both the scale matrix of the
#' density and the variance of the law, so [variance.MvGaussianDistrib()]
#' returns the same matrix; the two part company in the heavy-tailed sibling.
#'
#' @param distrib An [MvGaussianDistrib] object, from [mvgaussian1_distrib()].
#' @param theta A named list of parameters, one number each. The \eqn{p} mean
#'   components are ignored.
#'
#' @return A \eqn{p \times p} symmetric positive definite numeric matrix, with
#'   both dimnames `v1`, ..., `vp`.
#'
#' @seealso [variance.MvGaussianDistrib()] for the same matrix read as a
#'   moment, [mv_location.MvGaussianDistrib()] for the mean,
#'   [mv_summary()] for the standard deviations and correlations of a fit, and
#'   [mv_sigma()] for the generic.
#'
#' @examples
#' d <- mvgaussian1_distrib(2)
#' theta <- list(mu1 = 0.5, mu2 = -0.3, sigma_log_L1 = 0.1,
#'               sigma_log_L2 = -0.2, sigma_L2.1 = 0.4)
#' round(mv_sigma(d, theta), 4)
#'
#' # For a gaussian the scale matrix is the variance.
#' all.equal(mv_sigma(d, theta), variance(d, theta))
#'
#' # A precision parametrization reports the covariance too, inverted here.
#' o <- mvgaussian2_distrib(2, parameters7::log_cholesky(2))
#' th_o <- list(mu1 = 0, mu2 = 0, omega_log_L1 = 0.1,
#'              omega_log_L2 = -0.2, omega_L2.1 = 0.4)
#' Om <- parameters7::param_value(o@param, unlist(th_o)[3:5])
#' all.equal(mv_sigma(o, th_o), solve(Om), check.attributes = FALSE)
#'
#' @keywords internal
S7::method(mv_sigma, MvGaussianDistrib) <- function(distrib, theta) {
  pc <- mvg_pieces(distrib, theta)
  nm <- paste0("v", seq_len(distrib@n_dim))
  dimnames(pc$sigma) <- list(nm, nm)
  pc$sigma
}


#' @title Residuals of a Multivariate Gaussian and Their Image Under the Precision
#'
#' @description
#' Computes the centered response \eqn{r_i = y_i - \mu} and its image under the
#' inverse covariance, \eqn{w_i = \Sigma^{-1} r_i}. The score in the mean is
#' \eqn{w}, the response gradient is \eqn{-w}, and the derivatives in the
#' matrix are quadratic forms in \eqn{r}. Both are formed once per call.
#'
#' @param y An \eqn{n \times p} numeric matrix of observations, already
#'   coerced by [as_mv_matrix()].
#' @param pc The result of [mvg_pieces()], from which `mu` and `sigma_inv` are
#'   read.
#'
#' @return A named list with `r` and `w`, each an \eqn{n \times p} numeric
#'   matrix. Row \eqn{i} of `w` is \eqn{\Sigma^{-1}(y_i - \mu)}, the
#'   right-multiplication of `r` by `sigma_inv` giving the same rows because
#'   \eqn{\Sigma^{-1}} is symmetric.
#'
#' @section Notation:
#' \eqn{\mu} is the mean, \eqn{\Sigma} the covariance and \eqn{y_i} the
#' \eqn{i}-th row of the response.
#'
#' @seealso [mvg_pieces()] for the argument, and
#'   [distrib_gradient.MvGaussianDistrib()] and
#'   [distrib_hessian.MvGaussianDistrib()] for consumers.
#'
#' @examples
#' d <- mvgaussian1_distrib(2)
#' theta <- list(mu1 = 0.5, mu2 = -0.3, sigma_log_L1 = 0.1,
#'               sigma_log_L2 = -0.2, sigma_L2.1 = 0.4)
#' pc <- distributions7:::mvg_pieces(d, theta)
#' y <- rbind(c(0, 0), c(1, -1), c(2, 0.5))
#' res <- distributions7:::mvg_residuals(y, pc)
#'
#' # r is the centered response and w its image under the inverse covariance.
#' res$r
#' all.equal(res$w, res$r %*% solve(pc$sigma), check.attributes = FALSE)
#'
#' # The quadratic form of the log-density is the row sum of r * w.
#' all.equal(rowSums(res$r * res$w),
#'           mahalanobis(y, pc$mu, pc$sigma))
#'
#' # The score in the mean is w itself, one row per observation.
#' g <- distrib_gradient(d, y, theta)
#' all.equal(cbind(g$mu1, g$mu2), res$w, check.attributes = FALSE)
#'
#' @keywords internal
mvg_residuals <- function(y, pc) {
  r <- sweep(y, 2L, pc$mu, "-")
  list(r = r, w = r %*% pc$sigma_inv)
}


#' @title Multivariate Gaussian Density
#' @name distrib_pdf.MvGaussianDistrib
#'
#' @description
#' Computes the multivariate gaussian log-density
#' \deqn{\ell = -\frac{p}{2}\log 2\pi - \frac{1}{2}\log\lvert\Sigma\rvert
#'   - \frac{1}{2}(y-\mu)^\top \Sigma^{-1} (y-\mu),}
#' row by row, and exponentiates it unless `log = TRUE`. The log-determinant
#' comes from the matrix parametrization's own `param_logdet()`, so no
#' determinant is formed. The quadratic form is the squared norm of the residual
#' whitened by the parametrization's factor (see [mv_quad()]), which keeps its
#' digits where the covariance is nearly singular. The matrix quantities are
#' computed once per call, since the matrix does not vary with the response.
#'
#' @param distrib An [MvGaussianDistrib] object, from [mvgaussian1_distrib()].
#' @param y An \eqn{n \times p} numeric matrix of observations, one row each.
#'   A plain numeric vector of length \eqn{p} is read as a single observation.
#'   Every point of \eqn{\mathbb{R}^p} is in the support, so no row is
#'   rejected. A matrix with zero rows returns `numeric(0)`.
#' @param theta A named list of parameters, each component a single number:
#'   `mu1`, ..., `mup` and the matrix parametrization's prefixed free values.
#'   A parameter may not vary by observation here, and a component longer than
#'   one is an error.
#' @param log Logical of length 1. When `TRUE` the log-density is returned,
#'   which stays finite for a point the density itself underflows at. Defaults
#'   to `FALSE`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of length \eqn{n}, one density per row of `y`.
#'
#' @seealso [distrib_gradient.MvGaussianDistrib()] for the derivatives of this
#'   log-density, [distrib_rng.MvGaussianDistrib()] to draw from it,
#'   [mv_sigma.MvGaussianDistrib()] for the covariance it uses, and
#'   [distrib_pdf()] for the generic.
#'
#' @examples
#' d <- mvgaussian1_distrib(2)
#' theta <- list(mu1 = 0.5, mu2 = -0.3, sigma_log_L1 = 0.1,
#'               sigma_log_L2 = -0.2, sigma_L2.1 = 0.4)
#' y <- rbind(c(0, 0), c(1, -1), c(0.5, -0.3))
#'
#' distrib_pdf(d, y, theta)
#'
#' # Against the formula written out with an explicit inverse.
#' S <- mv_sigma(d, theta)
#' mu <- c(0.5, -0.3)
#' ref <- -0.5 * (2 * log(2 * pi) + log(det(S)) + mahalanobis(y, mu, S))
#' all.equal(distrib_pdf(d, y, theta, log = TRUE), as.numeric(ref))
#'
#' # A vector of length p is one observation, and the mode is the mean.
#' distrib_pdf(d, c(0.5, -0.3), theta, log = TRUE)
#'
#' # Far out in the tail the density underflows and its logarithm does not.
#' distrib_pdf(d, c(60, -60), theta)
#' distrib_pdf(d, c(60, -60), theta, log = TRUE)
#'
#' @keywords internal
S7::method(distrib_pdf, MvGaussianDistrib) <- function(distrib, y, theta,
                                                       log = FALSE, ...) {
  y <- as_mv_matrix(distrib, y)
  pc <- mvg_pieces(distrib, theta)
  if (!nrow(y)) return(numeric(0))
  res <- mvg_residuals(y, pc)
  q <- mv_quad(res$r, pc$s, pc$eta, distrib@inverted)
  ld <- -0.5 * (pc$p * base::log(2 * pi) + pc$logdet + q)
  if (log) ld else exp(ld)
}


#' @title Multivariate Gaussian Generator
#' @name distrib_rng.MvGaussianDistrib
#'
#' @description
#' Draws from the multivariate gaussian by the standard factorization, \eqn{y =
#' \mu + L z} with \eqn{z} a vector of independent standard normals and \eqn{L
#' L^\top = \Sigma}. The factor is read from the parametrization (see
#' [mv_sampling_factor()]): its own factor where it carries the covariance, the
#' inverse transpose of that factor where it carries the precision. A covariance
#' that is nearly singular, which `chol()` of the assembled matrix can reject,
#' is sampled without a decomposition. The draws consume `n * p` values from R's
#' own generator through [stats::rnorm()], so the stream is reproducible under
#' [base::set.seed()].
#'
#' @param distrib An [MvGaussianDistrib] object, from [mvgaussian1_distrib()].
#' @param n The number of observations to draw. A single non-negative whole
#'   number; `n = 0` returns a matrix with zero rows.
#' @param theta A named list of parameters, each component a single number.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return An \eqn{n \times p} numeric matrix, one draw per row, with column
#'   names `v1`, ..., `vp` and no row names.
#'
#' @seealso [distrib_pdf.MvGaussianDistrib()] for the density this draws from,
#'   [mv_sigma.MvGaussianDistrib()] for the covariance, and [distrib_rng()]
#'   for the generic.
#'
#' @examples
#' d <- mvgaussian1_distrib(2)
#' theta <- list(mu1 = 0.5, mu2 = -0.3, sigma_log_L1 = 0.1,
#'               sigma_log_L2 = -0.2, sigma_L2.1 = 0.4)
#'
#' set.seed(1)
#' y <- distrib_rng(d, 4, theta)
#' y
#'
#' # The sample moments approach the parameters, at the usual 1/sqrt(n) rate.
#' set.seed(2)
#' big <- distrib_rng(d, 20000, theta)
#' round(colMeans(big), 3)
#' round(var(big), 3)
#' round(mv_sigma(d, theta), 3)
#'
#' @keywords internal
S7::method(distrib_rng, MvGaussianDistrib) <- function(distrib, n, theta, ...) {
  pc <- mvg_pieces(distrib, theta)
  p <- pc$p
  l <- mv_sampling_factor(pc$s, pc$eta, distrib@inverted)
  z <- matrix(stats::rnorm(n * p), nrow = n, ncol = p)
  out <- z %*% t(l)
  out <- sweep(out, 2L, pc$mu, "+")
  colnames(out) <- paste0("v", seq_len(p))
  out
}


#' @title Multivariate Gaussian Score
#' @name distrib_gradient.MvGaussianDistrib
#'
#' @description
#' Computes the first derivatives of the log-density in closed form. With \eqn{r
#' = y - \mu}, \eqn{w = \Sigma^{-1}r} and \eqn{P_k} the derivative of
#' \eqn{\Sigma^{-1}} in the \eqn{k}-th free value of the matrix parametrization,
#' \deqn{\frac{\partial \ell}{\partial \mu} = w, \qquad \frac{\partial
#' \ell}{\partial \eta_k} = -\frac{1}{2}\frac{\partial
#' \log\lvert\Sigma\rvert}{\partial \eta_k} - \frac{1}{2}\, r^\top P_k r,} which
#' is \eqn{\tfrac{1}{2} w^\top A_k w} in the second term, with \eqn{A_k} the
#' derivative of \eqn{\Sigma}. The first term of the matrix component is the
#' parametrization's own `param_dlogdet()`, so no trace is formed here; its sign
#' is flipped where the parametrization carries the precision. \eqn{P_k} is read
#' from [mv_matrix_pieces()] and not formed as \eqn{-\Sigma^{-1}A_k\Sigma^{-1}},
#' a product that cancels where the covariance is nearly singular.
#'
#' @param distrib An [MvGaussianDistrib] object, from [mvgaussian1_distrib()].
#' @param y An \eqn{n \times p} numeric matrix of observations. A vector of
#'   length \eqn{p} is read as a single observation.
#' @param theta A named list of parameters, each component a single number.
#' @param scale One of `"parameter"` (the default) or `"link"`. The method
#'   returns the parameter scale and the generic applies the link scale to its
#'   result. Every link of this family is the identity, so the two scales
#'   coincide and the results are the same numbers to the bit.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A named list with one numeric vector of length \eqn{n} per
#'   parameter, in `distrib@params` order: \eqn{p} mean components followed by
#'   one per free value of the matrix parametrization.
#'
#' @section Notation:
#' \eqn{\mu} is the mean, \eqn{\Sigma} the covariance, \eqn{\eta} the free
#' vector of the matrix parametrization, \eqn{A_k =
#' \partial\Sigma/\partial\eta_k}, \eqn{P_k =
#' \partial\Sigma^{-1}/\partial\eta_k} and \eqn{\ell} the log-density of one
#' observation.
#'
#' @seealso [distrib_hessian.MvGaussianDistrib()] for the observed curvature,
#'   [distrib_expected_hessian.MvGaussianDistrib()] for the information,
#'   [distrib_grad_y.MvGaussianDistrib()] for the derivative in the response,
#'   and [distrib_gradient()] for the generic.
#'
#' @examples
#' d <- mvgaussian1_distrib(2)
#' theta <- list(mu1 = 0.5, mu2 = -0.3, sigma_log_L1 = 0.1,
#'               sigma_log_L2 = -0.2, sigma_L2.1 = 0.4)
#' set.seed(1)
#' y <- distrib_rng(d, 25, theta)
#'
#' g <- distrib_gradient(d, y, theta)
#' names(g)
#' vapply(g, sum, numeric(1))
#'
#' # Against a numerical derivative of the log-likelihood, which shares no
#' # arithmetic with the closed form.
#' ll <- function(v) {
#'   t2 <- as.list(v); names(t2) <- d@params
#'   sum(distrib_pdf(d, y, t2, log = TRUE))
#' }
#' max(abs(vapply(g, sum, numeric(1)) - numDeriv::grad(ll, unlist(theta))))
#'
#' # The mean component is w, one row per observation.
#' S <- mv_sigma(d, theta)
#' w <- sweep(y, 2, c(0.5, -0.3)) %*% solve(S)
#' all.equal(cbind(g$mu1, g$mu2), w, check.attributes = FALSE)
#'
#' # Every link is the identity, so the link scale changes nothing.
#' identical(g, distrib_gradient(d, y, theta, scale = "link"))
#'
#' @keywords internal
S7::method(distrib_gradient, MvGaussianDistrib) <- function(distrib, y, theta,
                                                            scale = c("parameter", "link"),
                                                            ...) {
  y <- as_mv_matrix(distrib, y)
  pc <- mvg_pieces(distrib, theta, derivs = TRUE)
  n <- nrow(y)
  res <- mvg_residuals(y, pc)
  w <- res$w
  r <- res$r

  out <- vector("list", distrib@n_params)
  names(out) <- distrib@params
  for (j in seq_len(pc$p)) out[[j]] <- w[, j]

  # w' A_k w = -r' d_k(Sigma^-1) r, read through the chart's derivative of the
  # inverse: the sandwich Sigma^-1 A_k Sigma^-1 cancels near a singular matrix
  dld <- pc$dlogdet
  for (k in seq_along(pc$ai)) {
    out[[pc$p + k]] <- -0.5 * dld[[k]] -
      0.5 * rowSums((r %*% pc$ai[[k]]) * r)
  }
  out
}


#' @title Multivariate Gaussian Observed Hessian
#' @name distrib_hessian.MvGaussianDistrib
#'
#' @description
#' Computes the second derivatives of the log-density in closed form. With
#' \eqn{r = y - \mu} and \eqn{P_k}, \eqn{P_{kl}} the first and second
#' derivatives of \eqn{\Sigma^{-1}} in the free values, \deqn{\ell^{(\mu_a
#' \mu_b)} = -(\Sigma^{-1})_{ab}, \qquad \ell^{(\mu_a \eta_k)} = (P_k r)_a,}
#' \deqn{\ell^{(\eta_k \eta_l)} = -\tfrac{1}{2}\frac{\partial^2
#' \log\lvert\Sigma\rvert}{\partial\eta_k \partial\eta_l} - \tfrac{1}{2}\,
#' r^\top P_{kl} r.} The mean block does not depend on the response, so it is
#' constant across rows; the mixed and matrix blocks do. \eqn{P_k} and
#' \eqn{P_{kl}} are read from [mv_matrix_pieces()], exact for
#' [parameters7::log_cholesky()] where the covariance is nearly singular;
#' written as products of \eqn{\Sigma^{-1}} and the derivatives of \eqn{\Sigma},
#' the matrix block cancels there for an observation away from the degenerate
#' direction.
#'
#' @param distrib An [MvGaussianDistrib] object, from [mvgaussian1_distrib()].
#' @param y An \eqn{n \times p} numeric matrix of observations. A vector of
#'   length \eqn{p} is read as a single observation.
#' @param theta A named list of parameters, each component a single number.
#' @param scale One of `"parameter"` (the default) or `"link"`. The method
#'   returns the parameter scale and the generic applies the link scale to its
#'   result. Every link of this family is the identity, so the two scales
#'   coincide.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A named list of numeric vectors of length \eqn{n}, one per unordered
#'   pair of parameters, keyed as [`hess_names(distrib@params)`][hess_names]:
#'   the diagonal first, then the off-diagonal pairs.
#'
#' @section Notation:
#' \eqn{\mu} is the mean, \eqn{\Sigma} the covariance, \eqn{\eta} the free
#' vector of the matrix parametrization, \eqn{P_k} and \eqn{P_{kl}} the first
#' and second derivatives of \eqn{\Sigma^{-1}}, and \eqn{\ell^{(ij)}} the second
#' derivative of the log-density in parameters \eqn{i} and \eqn{j}.
#'
#' @seealso [distrib_expected_hessian.MvGaussianDistrib()] for the expectation
#'   of this matrix, which is simpler, [distrib_gradient.MvGaussianDistrib()]
#'   for the score, and [distrib_hessian()] for the generic.
#'
#' @examples
#' d <- mvgaussian1_distrib(2)
#' theta <- list(mu1 = 0.5, mu2 = -0.3, sigma_log_L1 = 0.1,
#'               sigma_log_L2 = -0.2, sigma_L2.1 = 0.4)
#' set.seed(1)
#' y <- distrib_rng(d, 25, theta)
#'
#' H <- distrib_hessian(d, y, theta)
#' names(H)
#'
#' # The mean block is minus the inverse covariance and does not move with y.
#' c(H$mu1_mu1[1], H$mu1_mu2[1])
#' -solve(mv_sigma(d, theta))[1, ]
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
#' @keywords internal
S7::method(distrib_hessian, MvGaussianDistrib) <- function(distrib, y, theta,
                                                           scale = c("parameter", "link"),
                                                           ...) {
  y <- as_mv_matrix(distrib, y)
  pc <- mvg_pieces(distrib, theta, derivs2 = TRUE)
  n <- nrow(y)
  p <- pc$p
  r <- mvg_residuals(y, pc)$r
  si <- pc$sigma_inv

  d2ld <- pc$d2logdet
  # P_k = d_k(Sigma^-1) and P_kl, read from the chart: the mixed block is
  # (P_k r)_a and the matrix block -1/2 r' P_kl r, neither formed as a
  # sandwich of Sigma^-1, which cancels near a singular matrix
  rp <- lapply(pc$ai, function(ak) r %*% ak)
  spair <- param_pair_lookup(pc$s)

  pairs <- mv_hess_indices(distrib)
  out <- vector("list", length(pairs))
  names(out) <- hess_names(distrib@params)

  for (i in seq_along(pairs)) {
    a <- pairs[[i]][1L]
    b <- pairs[[i]][2L]
    out[[i]] <- if (a <= p && b <= p) {
      rep(-si[a, b], n)
    } else if (a <= p) {
      rp[[b - p]][, a]
    } else if (b <= p) {
      rp[[a - p]][, b]
    } else {
      k <- a - p
      l <- b - p
      idx <- spair[[paste(min(k, l), max(k, l), sep = ":")]]
      -0.5 * d2ld[[idx]] - 0.5 * rowSums((r %*% pc$ai2[[idx]]) * r)
    }
  }
  out
}


#' @title Index Pairs Behind the Hessian Keys of a Multivariate Distribution
#'
#' @description
#' Returns, for each key of [hess_names()], the two positions in
#' `distrib@params` that the key names, in the order in which the keys come. A
#' multivariate method branches on whether a component belongs to the mean
#' block, the matrix block or the mixed block, which depends on position.
#'
#' @details
#' [hess_pairs()] answers the same question with parameter names, which suits a
#' univariate method looking a closed-form component up in a table. Here the
#' names are matched back to positions, so the two orderings come from one
#' enumeration and cannot drift apart.
#'
#' @param distrib A [multivariate_distrib()] object, or any `distrib` whose
#'   `params` are the parameter names.
#'
#' @return A list of integer vectors of length 2, one per unordered pair, as
#'   long as `hess_names(distrib@params)` and in that order. Each pair is
#'   `(a, b)` with `a` and `b` positions in `distrib@params`.
#'
#' @seealso [hess_names()] for the keys, [hess_pairs()] for the same
#'   enumeration as names, and [distrib_hessian.MvGaussianDistrib()] for the
#'   consumer.
#'
#' @examples
#' d <- mvgaussian1_distrib(2)
#' idx <- distributions7:::mv_hess_indices(d)
#' length(idx) == length(hess_names(d@params))
#'
#' # Key by key, the positions name the parameters the key is built from.
#' head(data.frame(key = hess_names(d@params),
#'                 a = vapply(idx, `[`, integer(1), 1L),
#'                 b = vapply(idx, `[`, integer(1), 2L)), 6)
#'
#' # Which is the same enumeration hess_pairs() gives as names.
#' all.equal(lapply(idx, function(p) d@params[p]),
#'           distributions7:::hess_pairs(d@params),
#'           check.attributes = FALSE)
#'
#' @keywords internal
mv_hess_indices <- function(distrib) {
  params <- distrib@params
  lapply(hess_pairs(params), function(pr) match(pr, params))
}


#' @title Where Each Pair of Free Values Sits in a Parametrization's Second Derivatives
#'
#' @description
#' Builds a lookup from a pair of free-value positions to the position of the
#' matching component of the parametrization's second derivative arrays
#' (`param_d2()`, `param_inv_d2()` and `param_d2logdet()` share the
#' enumeration). A Hessian method walks the unordered pairs of the
#' distribution's own parameters and finds, for each matrix pair, the array that
#' belongs to it; the two enumerations are not the same, because the
#' distribution's parameters carry the mean components in front.
#'
#' @details
#' The keys are built from \pkg{parameters7}'s own enumeration rather than by
#' taking a component name apart, for the reason that package documents: a free
#' value whose label contains the separator splits into the wrong number of
#' pieces, so a name is not a safe route back to a pair of indices.
#'
#' @param s A \pkg{parameters7} parametrization, from which
#'   `param_tuple_indices()` supplies the enumeration.
#'
#' @return A named list of single integers, one per unordered pair of free
#'   values, keyed `"k:l"` with \eqn{k \le l}. Its length is
#'   `s@n_free * (s@n_free + 1) / 2`.
#'
#' @seealso [mv_hess_indices()] for the other half of the same bookkeeping and
#'   [distrib_hessian.MvGaussianDistrib()] for the consumer.
#'
#' @examples
#' s <- parameters7::log_cholesky(2)
#' lk <- distributions7:::param_pair_lookup(s)
#' unlist(lk)
#'
#' # The lookup is a permutation of the positions of param_d2()'s components.
#' d2 <- parameters7::param_d2(s, c(0.1, -0.2, 0.4))
#' setequal(unlist(lk), seq_along(d2))
#'
#' # The array it finds for the pair (1, 3) is the mixed second derivative in
#' # those two free values, against a difference of the first derivatives.
#' eta <- c(0.1, -0.2, 0.4)
#' h <- 1e-5
#' num <- (parameters7::param_d1(s, eta + c(0, 0, h))[[1]] -
#'         parameters7::param_d1(s, eta - c(0, 0, h))[[1]]) / (2 * h)
#' max(abs(d2[[lk[["1:3"]]]] - num))
#'
#' @keywords internal
param_pair_lookup <- function(s) {
  idx <- parameters7::param_tuple_indices(s)
  keys <- vapply(idx, function(kl) {
    paste(min(kl), max(kl), sep = ":")
  }, character(1))
  stats::setNames(as.list(seq_along(idx)), keys)
}


#' @title Multivariate Gaussian Expected Information
#' @name distrib_expected_hessian.MvGaussianDistrib
#'
#' @description
#' Computes the expectation of the observed Hessian in closed form, which for
#' this family is simpler than the observed matrix itself:
#' \deqn{\mathbb{E}[\ell^{(\mu_a \mu_b)}] = -(\Sigma^{-1})_{ab}, \qquad
#'   \mathbb{E}[\ell^{(\mu_a \eta_k)}] = 0, \qquad
#'   \mathbb{E}[\ell^{(\eta_k \eta_l)}] =
#'   -\tfrac{1}{2}\operatorname{tr}(\Sigma^{-1} A_k \Sigma^{-1} A_l).}
#' The components do not depend on the data, so every returned vector is
#' constant across rows and `y` is used only for its row count.
#'
#' @details
#' # Why the mixed block vanishes
#'
#' Each mixed component is a linear function of \eqn{w = \Sigma^{-1}(y-\mu)},
#' and \eqn{\mathbb{E}[w] = 0}. The mean parameters and the matrix parameters
#' are therefore orthogonal, and the information matrix is block diagonal
#' whatever the matrix parametrization is. A Fisher scoring step in the mean
#' therefore does not depend on the matrix parameters, and a step in the matrix
#' parameters does not depend on the mean.
#'
#' # The derivative arrays used
#'
#' The observed matrix block carries the second derivatives through both the
#' log-determinant term and the quadratic form. Under the expectation the two
#' terms cancel and only the first derivatives remain, so this method computes
#' no second derivative array.
#'
#' @param distrib An [MvGaussianDistrib] object, from [mvgaussian1_distrib()].
#' @param y An \eqn{n \times p} numeric matrix of observations. Only its row
#'   count is used: the expectation is taken over the law, so no observation
#'   enters any component.
#' @param theta A named list of parameters, each component a single number.
#' @param scale One of `"parameter"` (the default) or `"link"`. The method
#'   returns the parameter scale and the generic applies the link scale to its
#'   result. Every link of this family is the identity, so the two scales
#'   coincide.
#' @param approx Ignored: the expectation is exact and no approximation strategy
#'   is consulted. It is present so that the signature matches the generic's.
#' @param nsim Ignored, for the same reason.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A named list of numeric vectors of length \eqn{n}, keyed as
#'   [`hess_names(distrib@params)`][hess_names], each vector constant.
#'
#' @section Notation:
#' \eqn{\mu} is the mean, \eqn{\Sigma} the covariance, \eqn{\eta} the free
#' vector of the matrix parametrization, \eqn{A_k = \partial\Sigma/\partial\eta_k},
#' and \eqn{\ell^{(ij)}} the second derivative of the log-density in parameters
#' \eqn{i} and \eqn{j}.
#'
#' @seealso [distrib_hessian.MvGaussianDistrib()] for the observed matrix,
#'   [distrib_gradient.MvGaussianDistrib()] for the score, [fit_distrib()],
#'   whose Fisher scoring inverts this matrix, and
#'   [distrib_expected_hessian()] for the generic.
#'
#' @examples
#' d <- mvgaussian1_distrib(2)
#' theta <- list(mu1 = 0.5, mu2 = -0.3, sigma_log_L1 = 0.1,
#'               sigma_log_L2 = -0.2, sigma_L2.1 = 0.4)
#' y <- matrix(0, 3, 2)
#'
#' EH <- distrib_expected_hessian(d, y, theta)
#'
#' # The mean block is minus the inverse covariance.
#' c(EH$mu1_mu1[1], EH$mu1_mu2[1])
#' -solve(mv_sigma(d, theta))[1, ]
#'
#' # Every mixed mean-matrix component is exactly zero.
#' mixed <- grep("^mu[0-9]+_sigma", names(EH), value = TRUE)
#' vapply(EH[mixed], function(z) z[1], numeric(1))
#'
#' # And the closed form is what averaging the observed Hessian converges to.
#' set.seed(4)
#' big <- distrib_rng(d, 50000, theta)
#' emp <- vapply(distrib_hessian(d, big, theta), mean, numeric(1))
#' round(rbind(sampled = emp,
#'             closed = vapply(EH, function(z) z[1], numeric(1))), 3)
#'
#' @keywords internal
S7::method(distrib_expected_hessian, MvGaussianDistrib) <- function(
    distrib, y, theta, scale = c("parameter", "link"),
    approx = c("opg", "bartlett", "integrate", "mc"), nsim = 10000, ...) {
  y <- as_mv_matrix(distrib, y)
  pc <- mvg_pieces(distrib, theta, derivs = TRUE)
  n <- nrow(y)
  p <- pc$p
  si <- pc$sigma_inv
  sa <- lapply(pc$a, function(ak) si %*% ak)

  pairs <- mv_hess_indices(distrib)
  out <- vector("list", length(pairs))
  names(out) <- hess_names(distrib@params)

  for (i in seq_along(pairs)) {
    a <- pairs[[i]][1L]
    b <- pairs[[i]][2L]
    v <- if (a <= p && b <= p) {
      -si[a, b]
    } else if (a <= p || b <= p) {
      0
    } else {
      -0.5 * sum(sa[[a - p]] * t(sa[[b - p]]))
    }
    out[[i]] <- rep(v, n)
  }
  out
}


#' @title Multivariate Gaussian Response Gradient
#' @name distrib_grad_y.MvGaussianDistrib
#'
#' @description
#' Computes the derivative of the log-density in the response,
#' \deqn{\frac{\partial \ell}{\partial y} = -\Sigma^{-1}(y - \mu),}
#' one row per observation. This is \eqn{-w}, with \eqn{w = \Sigma^{-1}(y -
#' \mu)}, the negative of the score in the mean: for a location family the two
#' derivatives differ only in sign. For a univariate family the derivative in
#' the response is a numeric vector of length \eqn{n}; here it is an \eqn{n
#' \times p} matrix.
#'
#' @param distrib An [MvGaussianDistrib] object, from [mvgaussian1_distrib()].
#' @param y An \eqn{n \times p} numeric matrix of observations. A vector of
#'   length \eqn{p} is read as a single observation and gives a one-row result.
#' @param theta A named list of parameters, each component a single number.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return An \eqn{n \times p} numeric matrix, row \eqn{i} holding
#'   \eqn{\partial\ell_i/\partial y_i}.
#'
#' @seealso [distrib_hess_y.MvGaussianDistrib()] for the second derivative,
#'   [distrib_cross_y.MvGaussianDistrib()] for the mixed block, and
#'   [distrib_grad_y()] for the generic.
#'
#' @examples
#' d <- mvgaussian1_distrib(2)
#' theta <- list(mu1 = 0.5, mu2 = -0.3, sigma_log_L1 = 0.1,
#'               sigma_log_L2 = -0.2, sigma_L2.1 = 0.4)
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
#' # And it is minus the score in the mean, this being a location family.
#' g <- distrib_gradient(d, y, theta)
#' all.equal(distrib_grad_y(d, y, theta), -cbind(g$mu1, g$mu2),
#'           check.attributes = FALSE)
#'
#' @keywords internal
S7::method(distrib_grad_y, MvGaussianDistrib) <- function(distrib, y, theta, ...) {
  y <- as_mv_matrix(distrib, y)
  pc <- mvg_pieces(distrib, theta)
  -mvg_residuals(y, pc)$w
}


#' @title Multivariate Gaussian Response Hessian
#' @name distrib_hess_y.MvGaussianDistrib
#'
#' @description
#' Computes the second derivative of the log-density in the response,
#' \deqn{\frac{\partial^2 \ell}{\partial y\, \partial y^\top} = -\Sigma^{-1}.}
#' The log-density is quadratic in \eqn{y}, so the matrix is the same at every
#' observation and the response is not read. One \eqn{p \times p} matrix is
#' returned for the whole sample, and not one matrix per row.
#'
#' @param distrib An [MvGaussianDistrib] object, from [mvgaussian1_distrib()].
#' @param y An \eqn{n \times p} numeric matrix of observations. It is not read.
#' @param theta A named list of parameters, each component a single number.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A \eqn{p \times p} numeric matrix, symmetric and negative definite.
#'
#' @seealso [distrib_grad_y.MvGaussianDistrib()] for the first derivative,
#'   [distrib_cross2_y.MvGaussianDistrib()] for its derivative in the
#'   parameters, and [distrib_hess_y()] for the generic.
#'
#' @examples
#' d <- mvgaussian1_distrib(2)
#' theta <- list(mu1 = 0.5, mu2 = -0.3, sigma_log_L1 = 0.1,
#'               sigma_log_L2 = -0.2, sigma_L2.1 = 0.4)
#' set.seed(1)
#' y <- distrib_rng(d, 4, theta)
#'
#' distrib_hess_y(d, y, theta)
#' -solve(mv_sigma(d, theta))
#'
#' # One matrix whatever the sample size, the same at every observation.
#' identical(distrib_hess_y(d, y, theta), distrib_hess_y(d, y[1, ], theta))
#'
#' # Against a numerical Hessian at one observation.
#' max(abs(distrib_hess_y(d, y, theta) -
#'         numDeriv::hessian(function(z) distrib_pdf(d, z, theta, log = TRUE),
#'                           y[1, ])))
#'
#' @keywords internal
S7::method(distrib_hess_y, MvGaussianDistrib) <- function(distrib, y, theta, ...) {
  pc <- mvg_pieces(distrib, theta)
  -pc$sigma_inv
}


#' @title Multivariate Gaussian Mixed Response-Parameter Derivatives
#' @name distrib_cross_y.MvGaussianDistrib
#'
#' @description
#' Computes \eqn{\partial^2\ell/\partial y\, \partial\theta_k}, one \eqn{n
#' \times p} matrix per parameter. With \eqn{w = \Sigma^{-1}(y - \mu)} the
#' response gradient is \eqn{-w}, and differentiating it in the mean and in the
#' free values of the matrix parametrization gives \deqn{\frac{\partial^2
#' \ell}{\partial y \,\partial \mu_j} = \Sigma^{-1}e_j, \qquad \frac{\partial^2
#' \ell}{\partial y \,\partial \eta_k} = -P_k (y - \mu) = \Sigma^{-1}A_k w,}
#' with \eqn{A_k} and \eqn{P_k} the derivatives of \eqn{\Sigma} and of
#' \eqn{\Sigma^{-1}}. The mean block is column \eqn{j} of \eqn{\Sigma^{-1}} at
#' every observation; the matrix block carries the observation, and is computed
#' from \eqn{P_k} as [mv_matrix_pieces()] supplies it.
#'
#' @details
#' The result has one \eqn{n \times p} matrix per parameter, the shape of the
#' response, so the derivative in \eqn{\theta_k} of each row of the response
#' gradient is read row by row.
#'
#' Every link of this family is the identity, so `scale = "link"` and
#' `scale = "parameter"` give the same numbers.
#'
#' @param distrib An [MvGaussianDistrib] object, from [mvgaussian1_distrib()].
#' @param y An \eqn{n \times p} numeric matrix of observations. A vector of
#'   length \eqn{p} is read as a single observation.
#' @param theta A named list of parameters, each component a single number.
#' @param scale One of `"parameter"` (the default) or `"link"`. The method
#'   returns the parameter scale and the generic applies the link scale to its
#'   result. Every link of this family is the identity, so the two scales
#'   coincide.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A named list of \eqn{n \times p} numeric matrices, one per
#'   parameter, in `distrib@params` order.
#'
#' @section Notation:
#' \eqn{\mu} is the mean, \eqn{\Sigma} the covariance, \eqn{\eta} the free
#' vector of the matrix parametrization, \eqn{A_k =
#' \partial\Sigma/\partial\eta_k}, \eqn{P_k =
#' \partial\Sigma^{-1}/\partial\eta_k}, \eqn{e_j} the \eqn{j}-th standard basis
#' vector and \eqn{\ell} the log-density of one observation.
#'
#' @seealso [distrib_grad_y.MvGaussianDistrib()], whose derivative in the
#'   parameters this is, [distrib_cross2_y.MvGaussianDistrib()] for the next
#'   order, and [distrib_cross_y()] for the generic.
#'
#' @examples
#' d <- mvgaussian1_distrib(2)
#' theta <- list(mu1 = 0.5, mu2 = -0.3, sigma_log_L1 = 0.1,
#'               sigma_log_L2 = -0.2, sigma_L2.1 = 0.4)
#' set.seed(1)
#' y <- distrib_rng(d, 4, theta)
#'
#' cy <- distrib_cross_y(d, y, theta)
#' names(cy)
#' dim(cy$sigma_L2.1)
#'
#' # The mean block is a row of the inverse covariance, the same at every
#' # observation.
#' cy$mu1
#' solve(mv_sigma(d, theta))[1, ]
#'
#' # The matrix block against a difference of the response gradient.
#' h <- 1e-5
#' tp <- theta; tp$sigma_L2.1 <- tp$sigma_L2.1 + h
#' tm <- theta; tm$sigma_L2.1 <- tm$sigma_L2.1 - h
#' max(abs(cy$sigma_L2.1 -
#'         (distrib_grad_y(d, y, tp) - distrib_grad_y(d, y, tm)) / (2 * h)))
#'
#' @keywords internal
S7::method(distrib_cross_y, MvGaussianDistrib) <-
  function(distrib, y, theta, scale = c("parameter", "link"), ...) {
    y <- as_mv_matrix(distrib, y)
    p <- distrib@n_dim
    pc <- mvg_pieces(distrib, theta, derivs = TRUE)
    r <- mvg_residuals(y, pc)$r
    si <- pc$sigma_inv
    n <- nrow(y)

    mu_part <- lapply(seq_len(p), function(j) {
      matrix(rep(si[j, ], each = n), nrow = n, ncol = p)
    })
    # rows: w_i' A_k Sigma^-1 = -r_i' d_k(Sigma^-1), read through the chart's
    # derivative of the inverse rather than as a sandwich of Sigma^-1
    eta_part <- lapply(pc$ai, function(ak) -(r %*% ak))

    stats::setNames(c(mu_part, eta_part), distrib@params)
  }


#' @title Multivariate Gaussian Third Derivative in Two Responses and One Parameter
#' @name distrib_cross2_y.MvGaussianDistrib
#'
#' @description
#' Computes \eqn{\partial^3\ell/\partial y\,\partial y^\top \partial\theta_k},
#' one \eqn{p \times p} matrix per parameter. The response Hessian is
#' \eqn{-\Sigma^{-1}}, which does not involve the mean, so every mean component
#' is exactly the zero matrix; the matrix components are
#' \deqn{\frac{\partial^3\ell}{\partial y\,\partial y^\top\partial\eta_k}
#'   = B_k = \Sigma^{-1}A_k\Sigma^{-1},}
#' with \eqn{A_k = \partial\Sigma/\partial\eta_k}, which is \eqn{-P_k} and is
#' read as such from [mv_matrix_pieces()]. The components do not depend on the
#' observation, so one matrix is returned per parameter and `y` is not read.
#'
#' @details
#' This is one of the three derivatives that a marginal criterion reads when the
#' family is a prior over a coefficient block, the other two being
#' [distrib_hess_y_hess.MvGaussianDistrib()] and
#' [distrib_grad_y_hess.MvGaussianDistrib()].
#'
#' @param distrib An [MvGaussianDistrib] object, from [mvgaussian1_distrib()].
#' @param y An \eqn{n \times p} numeric matrix of observations. Not read: no
#'   component of this derivative depends on the response.
#' @param theta A named list of parameters, each component a single number.
#' @param scale One of `"parameter"` (the default) or `"link"`. The method
#'   returns the parameter scale and the generic applies the link scale to its
#'   result. Every link of this family is the identity, so the two scales
#'   coincide.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A named list of \eqn{p \times p} numeric matrices, one per
#'   parameter, in `distrib@params` order. The \eqn{p} mean components are the
#'   zero matrix.
#'
#' @section Notation:
#' \eqn{\Sigma} is the covariance, \eqn{\eta} the free vector of the matrix
#' parametrization, \eqn{A_k = \partial\Sigma/\partial\eta_k} and \eqn{\ell}
#' the log-density of one observation.
#'
#' @seealso [distrib_hess_y_hess.MvGaussianDistrib()] and
#'   [distrib_grad_y_hess.MvGaussianDistrib()] for the two fourth- and
#'   third-order siblings, [distrib_hess_y.MvGaussianDistrib()], whose
#'   derivative in the parameters this is, and [distrib_cross2_y()] for the
#'   generic.
#'
#' @examples
#' d <- mvgaussian1_distrib(2)
#' theta <- list(mu1 = 0.5, mu2 = -0.3, sigma_log_L1 = 0.1,
#'               sigma_log_L2 = -0.2, sigma_L2.1 = 0.4)
#' y <- matrix(0, 4, 2)
#'
#' c2 <- distrib_cross2_y(d, y, theta)
#' names(c2)
#'
#' # Every mean component is the zero matrix, the response Hessian carrying
#' # no mean at all.
#' c2$mu1
#'
#' # A matrix component against a difference of the response Hessian.
#' h <- 1e-5
#' tp <- theta; tp$sigma_L2.1 <- tp$sigma_L2.1 + h
#' tm <- theta; tm$sigma_L2.1 <- tm$sigma_L2.1 - h
#' max(abs(c2$sigma_L2.1 -
#'         (distrib_hess_y(d, y, tp) - distrib_hess_y(d, y, tm)) / (2 * h)))
#'
#' @keywords internal
S7::method(distrib_cross2_y, MvGaussianDistrib) <-
  function(distrib, y, theta, scale = c("parameter", "link"), ...) {
    p <- distrib@n_dim
    pc <- mvg_pieces(distrib, theta, derivs = TRUE)
    si <- pc$sigma_inv
    zero <- matrix(0, p, p)
    stats::setNames(
      c(rep(list(zero), p), lapply(pc$ai, function(ak) -ak)),
      distrib@params)
  }

#' @title Multivariate Gaussian Fourth Derivative in Two Responses and Two Parameters
#' @name distrib_hess_y_hess.MvGaussianDistrib
#'
#' @description
#' Computes \eqn{\partial^4\ell/\partial y\,\partial y^\top
#' \partial\theta_a\partial\theta_b}, one \eqn{p \times p} matrix per unordered
#' pair of parameters. The response Hessian \eqn{-\Sigma^{-1}} does not involve
#' the mean, so a pair naming any mean parameter gives the zero matrix; for a
#' pair of free values,
#' \deqn{\frac{\partial^4\ell}{\partial y\,\partial y^\top
#'     \partial\eta_k\partial\eta_l}
#'     = \Sigma^{-1}A_{kl}\Sigma^{-1}
#'       - \Sigma^{-1}\!\left(A_l\Sigma^{-1}A_k
#'         + A_k\Sigma^{-1}A_l\right)\!\Sigma^{-1},}
#' with \eqn{A_k} and \eqn{A_{kl}} the first and second derivatives of
#' \eqn{\Sigma}. The matrix is \eqn{-P_{kl}}, the second derivative of
#' \eqn{\Sigma^{-1}}, and is read as such from [mv_matrix_pieces()]. The
#' components do not depend on the observation, so `y` is not read.
#'
#' @param distrib An [MvGaussianDistrib] object, from [mvgaussian1_distrib()].
#' @param y An \eqn{n \times p} numeric matrix of observations. It is not read.
#' @param theta A named list of parameters, each component a single number.
#' @param scale One of `"parameter"` (the default) or `"link"`. The method
#'   returns the parameter scale and the generic applies the link scale to its
#'   result. Every link of this family is the identity, so the two scales
#'   coincide.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A named list of \eqn{p \times p} numeric matrices, keyed as
#'   [`hess_names(distrib@params)`][hess_names]. Every key naming a mean
#'   parameter holds the zero matrix.
#'
#' @section Notation:
#' \eqn{\Sigma} is the covariance, \eqn{\eta} the free vector of the matrix
#' parametrization, \eqn{A_k} and \eqn{A_{kl}} its first and second derivative
#' arrays, and \eqn{\ell} the log-density of one observation.
#'
#' @seealso [distrib_cross2_y.MvGaussianDistrib()] for the same derivative one
#'   parameter down, [distrib_grad_y_hess.MvGaussianDistrib()] for its sibling
#'   with one response index, and [distrib_hess_y_hess()] for the generic.
#'
#' @examples
#' d <- mvgaussian1_distrib(2)
#' theta <- list(mu1 = 0.5, mu2 = -0.3, sigma_log_L1 = 0.1,
#'               sigma_log_L2 = -0.2, sigma_L2.1 = 0.4)
#' y <- matrix(0, 4, 2)
#'
#' hh <- distrib_hess_y_hess(d, y, theta)
#' length(hh) == length(hess_names(d@params))
#'
#' # Any pair naming a mean gives the zero matrix.
#' hh$mu1_sigma_L2.1
#'
#' # A matrix pair against a second difference of the response Hessian.
#' h <- 1e-4
#' f <- function(a, b) {
#'   t2 <- theta
#'   t2$sigma_log_L1 <- t2$sigma_log_L1 + a
#'   t2$sigma_L2.1 <- t2$sigma_L2.1 + b
#'   distrib_hess_y(d, y, t2)
#' }
#' num <- (f(h, h) - f(h, -h) - f(-h, h) + f(-h, -h)) / (4 * h * h)
#' round(hh$sigma_log_L1_sigma_L2.1, 6)
#' round(num, 6)
#'
#' @keywords internal
S7::method(distrib_hess_y_hess, MvGaussianDistrib) <-
  function(distrib, y, theta, scale = c("parameter", "link"), ...) {
    p <- distrib@n_dim
    pc <- mvg_pieces(distrib, theta, derivs2 = TRUE)
    si <- pc$sigma_inv
    nm <- distrib@params
    zero <- matrix(0, p, p)
    stats::setNames(lapply(hess_names(nm), function(k) {
      ij <- hess_pairs(nm)[[k]]
      a <- match(ij[1L], nm)
      b <- match(ij[2L], nm)
      if (a <= p || b <= p) return(zero)
      ka <- a - p
      kb <- b - p
      -mv_ai2(pc, ka, kb)
    }), hess_names(nm))
  }

#' @title Multivariate Gaussian Third Derivative in One Response and Two Parameters
#' @name distrib_grad_y_hess.MvGaussianDistrib
#'
#' @description
#' Computes \eqn{\partial^3\ell/\partial y\,\partial\theta_a\partial\theta_b},
#' one \eqn{n \times p} matrix per unordered pair of parameters. The response
#' gradient \eqn{-\Sigma^{-1}(y-\mu)} is linear in the mean, so a pair naming
#' two mean parameters gives the zero matrix. The mixed pairs and the pure
#' matrix pairs are
#' \deqn{\frac{\partial^3\ell}{\partial y\,\partial\mu_j\partial\eta_k}
#'     = -B_k e_j, \qquad
#'   \frac{\partial^3\ell}{\partial y\,\partial\eta_k\partial\eta_l}
#'     = \Sigma^{-1}A_{kl}w
#'       - \Sigma^{-1}\!\left(A_l\Sigma^{-1}A_k
#'         + A_k\Sigma^{-1}A_l\right)\!w,}
#' with \eqn{B_k = \Sigma^{-1}A_k\Sigma^{-1}} and
#' \eqn{w = \Sigma^{-1}(y-\mu)}; the two are \eqn{P_k e_j} and \eqn{-P_{kl}(y -
#' \mu)}, computed from the derivatives of \eqn{\Sigma^{-1}} that
#' [mv_matrix_pieces()] supplies. Only the pure matrix pairs carry the
#' observation; the mixed pairs repeat one row.
#'
#' @param distrib An [MvGaussianDistrib] object, from [mvgaussian1_distrib()].
#' @param y An \eqn{n \times p} numeric matrix of observations. A vector of
#'   length \eqn{p} is read as a single observation.
#' @param theta A named list of parameters, each component a single number.
#' @param scale One of `"parameter"` (the default) or `"link"`. The method
#'   returns the parameter scale and the generic applies the link scale to its
#'   result. Every link of this family is the identity, so the two scales
#'   coincide.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A named list of \eqn{n \times p} numeric matrices, keyed as
#'   [`hess_names(distrib@params)`][hess_names]. Every key naming two mean
#'   parameters holds a matrix of zeros.
#'
#' @section Notation:
#' \eqn{\mu} is the mean, \eqn{\Sigma} the covariance, \eqn{\eta} the free
#' vector of the matrix parametrization, \eqn{A_k} and \eqn{A_{kl}} its first
#' and second derivative arrays, \eqn{e_j} the \eqn{j}-th standard basis vector
#' and \eqn{\ell} the log-density of one observation.
#'
#' @seealso [distrib_cross2_y.MvGaussianDistrib()] and
#'   [distrib_hess_y_hess.MvGaussianDistrib()] for the other two derivatives that
#'   a marginal criterion reads, and [distrib_grad_y_hess()] for the generic.
#'
#' @examples
#' d <- mvgaussian1_distrib(2)
#' theta <- list(mu1 = 0.5, mu2 = -0.3, sigma_log_L1 = 0.1,
#'               sigma_log_L2 = -0.2, sigma_L2.1 = 0.4)
#' set.seed(1)
#' y <- distrib_rng(d, 4, theta)
#'
#' gh <- distrib_grad_y_hess(d, y, theta)
#' dim(gh$sigma_log_L1_sigma_L2.1)
#'
#' # Two mean parameters give exactly zero, the response gradient being linear
#' # in the mean.
#' gh$mu1_mu2
#'
#' # A mixed pair repeats one row; a matrix pair does not.
#' gh$mu1_sigma_L2.1
#' round(gh$sigma_log_L1_sigma_L2.1, 4)
#'
#' # Against a second difference of the response gradient.
#' h <- 1e-4
#' f <- function(a, b) {
#'   t2 <- theta
#'   t2$sigma_log_L1 <- t2$sigma_log_L1 + a
#'   t2$sigma_L2.1 <- t2$sigma_L2.1 + b
#'   distrib_grad_y(d, y, t2)
#' }
#' max(abs(gh$sigma_log_L1_sigma_L2.1 -
#'         (f(h, h) - f(h, -h) - f(-h, h) + f(-h, -h)) / (4 * h * h)))
#'
#' @keywords internal
S7::method(distrib_grad_y_hess, MvGaussianDistrib) <-
  function(distrib, y, theta, scale = c("parameter", "link"), ...) {
    y <- as_mv_matrix(distrib, y)
    p <- distrib@n_dim
    n <- nrow(y)
    pc <- mvg_pieces(distrib, theta, derivs2 = TRUE)
    r <- mvg_residuals(y, pc)$r
    nm <- distrib@params
    zero <- matrix(0, n, p)
    stats::setNames(lapply(hess_names(nm), function(k) {
      ij <- hess_pairs(nm)[[k]]
      a <- match(ij[1L], nm)
      b <- match(ij[2L], nm)
      # the response gradient is -Sigma^-1 r, linear in the mean, so a
      # component naming two means vanishes
      if (a <= p && b <= p) return(zero)
      if (a <= p || b <= p) {
        j <- if (a <= p) a else b
        kk <- if (a <= p) b - p else a - p
        # d_k(Sigma^-1) e_j, the same row at every observation
        m <- pc$ai[[kk]]
        return(matrix(rep(m[, j], each = n), nrow = n, ncol = p))
      }
      ka <- a - p
      kb <- b - p
      -(r %*% mv_ai2(pc, ka, kb))
    }), hess_names(nm))
  }

# Commented out on 2026-10-12 (T4): the observed Hessian reads the chart's second derivatives of the inverse,
# so the second derivative of the covariance has no caller here.
# Kept rather than deleted.
# @title The Second Derivative of the Covariance, by Position
#
# @description
# Reads the component of `parameters7::param_d2()` belonging to the pair of
# free values \eqn{(k, l)}. The key is constructed from the sorted pair of
# free names, never parsed out of a stored key, so a free value whose own
# label contains the separator cannot split into the wrong number of pieces.
# The result is symmetric in its two positions.
#
# @param pc The result of [mvg_pieces()] called with `derivs2 = TRUE`, whose
#   `a2` component holds the second derivative arrays and whose `s` component
#   supplies the free names.
# @param k,l Positions among the matrix parametrization's free values, each a
#   single whole number between 1 and `s@n_free`. Their order does not matter.
#
# @return A \eqn{p \times p} symmetric numeric matrix,
#   \eqn{\partial^2\Sigma/\partial\eta_k\partial\eta_l}.
#
# @seealso [mvg_pieces()] for the argument and [param_pair_lookup()] for the
#   other route to the same enumeration.
#
# @examples
# d <- mvgaussian1_distrib(2)
# theta <- list(mu1 = 0, mu2 = 0, sigma_log_L1 = 0.1,
#               sigma_log_L2 = -0.2, sigma_L2.1 = 0.4)
# pc <- distributions7:::mvg_pieces(d, theta, derivs2 = TRUE)
#
# distributions7:::mvg_a2(pc, 1L, 3L)
#
# # Symmetric in the two positions.
# identical(distributions7:::mvg_a2(pc, 1L, 3L),
#           distributions7:::mvg_a2(pc, 3L, 1L))
#
# # And it is the mixed second derivative, against a difference of the first.
# h <- 1e-5
# eta <- c(0.1, -0.2, 0.4)
# s <- d@param
# num <- (parameters7::param_d1(s, eta + c(0, 0, h))[[1]] -
#         parameters7::param_d1(s, eta - c(0, 0, h))[[1]]) / (2 * h)
# max(abs(distributions7:::mvg_a2(pc, 1L, 3L) - num))
#
# @keywords internal
# mvg_a2 <- function(pc, k, l) {
#   nm <- pc$s@free_names
#   ij <- sort(c(k, l))
#   pc$a2[[paste(nm[ij], collapse = ":")]]
# }


#' @title The Second Derivative of the Inverse, by Position
#'
#' @description
#' Reads the component of the second derivative array of \eqn{\Sigma^{-1}}
#' that [mv_matrix_pieces()] supplies, for the pair of free values
#' \eqn{(k, l)}. The key is built from the sorted pair of free names, so the
#' result does not depend on the order of the two positions. Shared by the
#' gaussian and the Student t.
#'
#' @param pc Pieces carrying `ai2`, from [mv_matrix_pieces()] with `derivs2`.
#' @param k,l Positions among the free values of the matrix parametrization.
#'
#' @return A \eqn{p \times p} numeric matrix.
#'
#' @keywords internal
mv_ai2 <- function(pc, k, l) {
  nm <- pc$s@free_names
  ij <- sort(c(k, l))
  unname(pc$ai2[[paste(nm[ij], collapse = ":")]])
}

#' @title Mean of a Multivariate Gaussian
#' @name mean.MvGaussianDistrib
#'
#' @description
#' Returns \eqn{\mathbb{E}[Y] = \mu}, the mean vector. For this family the
#' expectation is a parameter and needs no integration: the density is
#' symmetric about \eqn{\mu} and every coordinate has finite moments of every
#' order, so the first moment exists at every parameter value the constructor
#' admits.
#'
#' @param x An [MvGaussianDistrib] object, from [mvgaussian1_distrib()].
#' @param theta A named list of parameters, each component a single number.
#'   Only the \eqn{p} mean components are read.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of length \eqn{p}, named `v1`, ..., `vp`.
#'
#' @seealso [variance.MvGaussianDistrib()] for the second moment,
#'   [mv_location.MvGaussianDistrib()], which returns the same vector as the
#'   density's center, and [base::mean()] for the generic.
#'
#' @examples
#' d <- mvgaussian1_distrib(2)
#' theta <- list(mu1 = 0.5, mu2 = -0.3, sigma_log_L1 = 0.1,
#'               sigma_log_L2 = -0.2, sigma_L2.1 = 0.4)
#'
#' mean(d, theta)
#'
#' # Which is what a large sample average approaches.
#' set.seed(2)
#' round(colMeans(distrib_rng(d, 20000, theta)), 3)
#'
#' @keywords internal
S7::method(mean, MvGaussianDistrib) <- function(x, theta, ...) {
  mv_location(x, theta)
}


#' @title Variance of a Multivariate Gaussian
#' @name variance.MvGaussianDistrib
#'
#' @description
#' Returns \eqn{\operatorname{Var}(Y) = \Sigma}, the covariance matrix the
#' matrix parametrization carries, inverted first where that parametrization
#' carries the precision. For this family the scale matrix of the density is
#' the covariance of the law, so the value agrees with
#' [mv_sigma.MvGaussianDistrib()]; the two part company in the heavy-tailed
#' sibling, where the scale matrix exists at every \eqn{\nu} and the covariance
#' does not.
#'
#' The result is a matrix, the second central moment of a vector being a matrix,
#' where a univariate family returns a numeric vector.
#'
#' @param x An [MvGaussianDistrib] object, from [mvgaussian1_distrib()].
#' @param theta A named list of parameters, each component a single number.
#'   The \eqn{p} mean components are ignored.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A \eqn{p \times p} symmetric positive definite numeric matrix, with
#'   both dimnames `v1`, ..., `vp`.
#'
#' @seealso [mean.MvGaussianDistrib()] for the first moment,
#'   [mv_sigma.MvGaussianDistrib()] for the same matrix read as the density's
#'   scale, [variance.MvStudentTDistrib()] for the family where the two
#'   differ, and [variance()] for the generic.
#'
#' @examples
#' d <- mvgaussian1_distrib(2)
#' theta <- list(mu1 = 0.5, mu2 = -0.3, sigma_log_L1 = 0.1,
#'               sigma_log_L2 = -0.2, sigma_L2.1 = 0.4)
#'
#' variance(d, theta)
#'
#' # The scale matrix and the covariance are one matrix for this family.
#' identical(variance(d, theta), mv_sigma(d, theta))
#'
#' # Which a large sample covariance approaches.
#' set.seed(2)
#' round(var(distrib_rng(d, 20000, theta)), 3)
#'
#' @keywords internal
S7::method(variance, MvGaussianDistrib) <- function(x, theta, ...) {
  mv_sigma(x, theta)
}


#' @title Random Parameters for a Multivariate Gaussian
#' @name generate_random_theta.MvGaussianDistrib
#'
#' @description
#' Draws a parameter vector for testing: each mean uniform on \eqn{(-1, 1)} and
#' each free value of the matrix parametrization uniform on \eqn{(-0.4, 0.4)},
#' which for the default log-Cholesky parametrization puts the matrix near the
#' identity. The method of the base class draws every unbounded parameter
#' uniformly on \eqn{(-3, 3)}, which on a log-Cholesky diagonal gives variances
#' between \eqn{e^{-6}} and \eqn{e^{6}} and a poorly conditioned covariance; the
#' narrow band keeps the matrices that [check_distrib()] draws for this family
#' well conditioned.
#'
#' @param distrib An [MvGaussianDistrib] object, from [mvgaussian1_distrib()].
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A named list of `distrib@n_params` single numbers, named and ordered
#'   as `distrib@params`.
#'
#' @seealso [check_distrib()], which draws parameters this way, and
#'   [generate_random_theta()] for the generic.
#'
#' @examples
#' d <- mvgaussian1_distrib(2)
#'
#' set.seed(11)
#' unlist(generate_random_theta(d))
#'
#' # The matrix components stay inside (-0.4, 0.4) and the means inside
#' # (-1, 1), so the covariance drawn is never far from the identity.
#' set.seed(12)
#' many <- replicate(500, unlist(generate_random_theta(d)))
#' round(apply(many, 1, range), 2)
#'
#' # Every draw gives a positive definite covariance, the parametrization
#' # having no boundary to reach.
#' set.seed(13)
#' all(replicate(50, min(eigen(mv_sigma(d, generate_random_theta(d)),
#'                             only.values = TRUE)$values) > 0))
#'
#' @keywords internal
S7::method(generate_random_theta, MvGaussianDistrib) <- function(distrib, ...) {
  p <- distrib@n_dim
  s <- distrib@param
  as.list(stats::setNames(
    c(stats::runif(p, -1, 1), stats::runif(s@n_free, -0.4, 0.4)),
    distrib@params
  ))
}


#' @title Precision Derivative Arrays of a Multivariate Gaussian
#'
#' @description
#' Supplies the derivative arrays of the precision \eqn{P = \Sigma^{-1}} in the
#' free values of the matrix parametrization, to orders 1 through 4, as an
#' accessor keyed by an index multiset. Under a precision parametrization they
#' are the parametrization's own `param_d1()` to `param_d4()`; under a
#' covariance parametrization they are its derivatives of the inverse,
#' [parameters7::param_inv_d1()] to [parameters7::param_inv_d4()]. For
#' [parameters7::log_cholesky()] these are exact where the covariance is nearly
#' singular; for another parametrization they are the ordered-partition sum of
#' [parameters7::inverse_of()].
#'
#' @details
#' The accessor returns \eqn{P} itself for the empty multiset. The gaussian
#' does not use it; the multivariate Student t does, for a partition block of
#' mean indices that carries no matrix index. The function takes the pieces,
#' so that the Student t passes it the arrays of its scale matrix and the
#' package keeps one copy of the construction.
#'
#' @param pc The pieces, as returned by [mvg_pieces()] or `mvt_pieces()`: any
#'   list carrying `s` (the parametrization), `eta` (its free vector) and
#'   `sigma_inv`.
#' @param order The highest order wanted, a single whole number from 1 to 4.
#'   Every array up to that order is computed.
#' @param inverted Logical of length 1. `TRUE` when the free values
#'   parametrize the precision. Defaults to `FALSE`.
#'
#' @return A named list with `get`, a function of one integer vector returning
#'   the \eqn{p \times p} array for that multiset; `sign_ld`, the coefficient
#'   that the log-determinant term carries in a derivative of the log-density
#'   (\eqn{-1/2} for a covariance, \eqn{+1/2} for a precision); and `pc`, the
#'   pieces as supplied.
#'
#' @section Notation:
#' \eqn{\Sigma} is the covariance, \eqn{P = \Sigma^{-1}} the precision,
#' \eqn{\eta} the free vector of the matrix parametrization and \eqn{P_t} the
#' derivative of the precision over the multiset \eqn{t}.
#'
#' @seealso [mvg_higher()] and [mvt_higher()] for the consumers, and
#'   [mvg_pieces()] for the argument.
#'
#' @examples
#' d <- mvgaussian1_distrib(2)
#' theta <- list(mu1 = 0, mu2 = 0, sigma_log_L1 = 0.1,
#'               sigma_log_L2 = -0.2, sigma_L2.1 = 0.4)
#' pc <- distributions7:::mvg_pieces(d, theta, derivs2 = TRUE)
#' pt <- distributions7:::mvg_ptensors(pc, 2L)
#'
#' P <- pc$sigma_inv
#' A1 <- pc$a[[1]]
#'
#' # The empty multiset is the precision itself.
#' all.equal(pt$get(integer(0)), P)
#'
#' # At order one it is the derivative of an inverse, -P A P.
#' all.equal(pt$get(1L), -P %*% A1 %*% P)
#'
#' # A precision parametrization reads its own arrays and flips the sign that
#' # the log-determinant term carries.
#' c(covariance = pt$sign_ld,
#'   precision = distributions7:::mvg_ptensors(pc, 2L, inverted = TRUE)$sign_ld)
#'
#' @keywords internal
mvg_ptensors <- function(pc, order, inverted = FALSE) {
  s <- pc$s
  eta <- pc$eta
  # the precision's derivatives: the chart's own arrays when it carries the
  # precision, the chart's derivatives of its inverse when it carries the
  # covariance (exact for log_cholesky(), the ordered-partition sum of
  # inverse_of() otherwise)
  dfun <- if (inverted) {
    list(parameters7::param_d1, parameters7::param_d2,
         parameters7::param_d3, parameters7::param_d4)
  } else {
    list(parameters7::param_inv_d1, parameters7::param_inv_d2,
         parameters7::param_inv_d3, parameters7::param_inv_d4)
  }

  key <- function(t) paste(sort(t), collapse = ",")
  amap <- list()
  for (o in seq_len(order)) {
    arr <- dfun[[o]](s, eta)
    idx <- if (o == 1L) {
      lapply(seq_len(s@n_free), identity)
    } else {
      parameters7::param_tuple_indices(s, o)
    }
    for (i in seq_along(idx)) {
      amap[[paste0(o, ":", key(idx[[i]]))]] <- unname(unclass(arr[[i]]))
    }
  }
  P <- unname(pc$sigma_inv)
  # the empty multiset is the precision itself; the gaussian never asks for
  # it, but the Student t does, a partition block of pure mean indices
  # carrying no matrix index at all
  get <- function(t) {
    if (!length(t)) return(P)
    amap[[paste0(length(t), ":", key(t))]]
  }
  list(get = get, sign_ld = if (inverted) 0.5 else -0.5, pc = pc)
}


# Commented out on 2026-10-12 (T4): mvg_ptensors() reads the precision's derivative arrays from the
# chart (parameters7::param_inv_d1() to param_inv_d4()), so this has no caller.
# The same ordered-partition sum is parameters7's inv_partition_sum().
# Kept rather than deleted.
# @title Ordered Partitions of a Set of Positions
#
# @description
# Enumerates all ordered partitions of a set of positions into nonempty
# blocks: every set partition, in every ordering of its blocks. This is how
# the derivative of an inverse distributes its differentiations,
# \deqn{P_t = \sum (-1)^q\, P A_{B_1} P \cdots A_{B_q} P,}
# the blocks appearing in a product where order matters even though the
# partition itself is unordered. The count is the ordered Bell number: 1, 3,
# 13, 75 at sizes 1 to 4.
#
# @details
# The set partitions are built by inserting the last element into each block
# of each partition of the rest, and into a block of its own; the orderings by
# inserting the last index at each position of each permutation of the rest.
# Both recursions are written out here, the sizes involved being at most 4.
# `numericals7::set_partitions()` supplies the unordered enumeration the rest
# of the toolkit uses, and would still leave the block orderings to do.
#
# @param pos An integer vector of positions, of length 1 to 4 in practice.
#   Its values are carried through unchanged, so it may hold any labels the
#   caller indexes with.
#
# @return A list of ordered partitions. Each is a list of integer vectors, the
#   blocks in the order the product takes them, and the blocks of one
#   partition together hold every element of `pos` exactly once.
#
# @section Notation:
# \eqn{P} is a precision matrix, \eqn{A_B} the derivative of its inverse in
# the free values \eqn{B}, and \eqn{q} the number of blocks of a partition.
#
# @seealso [mvg_ptensors()], the only consumer, and
#   [numericals7::set_partitions()] for the unordered enumeration the toolkit
#   uses elsewhere.
#
# @examples
# # Two positions give three ordered partitions: one block, then the two
# # orderings of two singleton blocks.
# distributions7:::mv_ordered_partitions(1:2)
#
# # The counts are the ordered Bell numbers.
# vapply(1:4, function(n)
#   length(distributions7:::mv_ordered_partitions(seq_len(n))), integer(1))
#
# # Every partition covers the positions exactly once.
# all(vapply(distributions7:::mv_ordered_partitions(1:3),
#            function(p) identical(sort(unlist(p)), 1:3), TRUE))
#
# @keywords internal
# mv_ordered_partitions <- function(pos) {
#   if (length(pos) == 1L) return(list(list(pos)))
#   # set partitions by recursive insertion of the last element
#   set_parts <- function(v) {
#     if (length(v) == 1L) return(list(list(v)))
#     prev <- set_parts(v[-length(v)])
#     out <- list()
#     x <- v[length(v)]
#     for (pp in prev) {
#       for (j in seq_along(pp)) {
#         q <- pp
#         q[[j]] <- c(q[[j]], x)
#         out[[length(out) + 1L]] <- q
#       }
#       out[[length(out) + 1L]] <- c(pp, list(x))
#     }
#     out
#   }
#   perms <- function(n) {
#     if (n == 1L) return(list(1L))
#     out <- list()
#     for (pr in perms(n - 1L)) {
#       for (j in seq_len(n)) {
#         out[[length(out) + 1L]] <- append(pr, n, after = j - 1L)
#       }
#     }
#     out
#   }
#   out <- list()
#   for (pp in set_parts(pos)) {
#     q <- length(pp)
#     for (pr in perms(q)) {
#       out[[length(out) + 1L]] <- pp[pr]
#     }
#   }
#   out
# }


#' @title The Closed-Form Higher Derivatives of a Multivariate Gaussian
#'
#' @description
#' Computes every third- or fourth-order derivative of the log-density in the
#' parameters. It enumerates the parameter tuples the way [deriv_names()] does,
#' splits each into mean indices and matrix indices, and reads the surviving
#' cases off the gaussian's algebra. Writing \eqn{r = y - \mu} and \eqn{P_t}
#' for the precision's derivative array over the matrix indices \eqn{t}, a
#' tuple with
#'
#' - three or more mean indices is exactly zero, the quadratic form being
#'   quadratic in \eqn{\mu};
#' - two mean indices \eqn{(i, j)} gives \eqn{-P_t[i, j]};
#' - one mean index \eqn{i} gives \eqn{(r P_t)_i}, one value per observation;
#' - no mean index gives
#'   \eqn{\pm\tfrac{1}{2}\,\partial^{\lvert t\rvert}\log\lvert M\rvert
#'   - \tfrac{1}{2}\, r^\top P_t\, r}, the sign following which side the
#'   parametrization carries.
#'
#' @details
#' The log-determinant derivatives come from `parameters7::param_d3logdet()` and
#' `param_d4logdet()`, and the precision arrays from [mvg_ptensors()]. With
#' `expected = TRUE` the response enters through \eqn{r} and \eqn{r^\top P_t r},
#' whose expectations are 0 and \eqn{\operatorname{tr}(P_t\Sigma)}, so the
#' expectation is exact.
#'
#' @param distrib An [MvGaussianDistrib] object, from [mvgaussian1_distrib()].
#' @param y An \eqn{n \times p} numeric matrix of observations. A vector of
#'   length \eqn{p} is read as a single observation.
#' @param theta A named list of parameters, each component a single number.
#' @param order The derivative order, `3L` or `4L`.
#' @param expected Logical of length 1. When `TRUE` the expectation of each
#'   component is returned. Defaults to `FALSE`.
#'
#' @return A named list of numeric vectors of length \eqn{n}, one per
#'   derivative component, keyed and ordered as
#'   `deriv_names(distrib@params, order)`. There are 35 components at order 3
#'   and 70 at order 4 for a two-dimensional gaussian on an unstructured
#'   covariance.
#'
#' @section Notation:
#' \eqn{\mu} is the mean, \eqn{M} the matrix the parametrization carries,
#' \eqn{P = \Sigma^{-1}} the precision, \eqn{\eta} the free vector,
#' \eqn{r = y - \mu} the centered response and \eqn{P_t} the precision's
#' derivative array over the multiset of free values \eqn{t}.
#'
#' @seealso [distrib_deriv3.MvGaussianDistrib()] and
#'   [distrib_deriv4.MvGaussianDistrib()], the two methods that it serves, and
#'   [mvg_ptensors()] for the arrays that it reads.
#'
#' @examples
#' d <- mvgaussian1_distrib(2)
#' theta <- list(mu1 = 0.5, mu2 = -0.3, sigma_log_L1 = 0.1,
#'               sigma_log_L2 = -0.2, sigma_L2.1 = 0.4)
#' set.seed(1)
#' y <- distrib_rng(d, 4, theta)
#'
#' d3 <- distributions7:::mvg_higher(d, y, theta, 3L)
#' length(d3)
#' length(distributions7:::mvg_higher(d, y, theta, 4L))
#'
#' # Three mean indices vanish exactly.
#' vapply(d3[grep("^mu[0-9]+_mu[0-9]+_mu[0-9]+$", names(d3))],
#'        function(z) max(abs(z)), numeric(1))
#'
#' # Against one stencil on the analytic Hessian, which shares no algebra
#' # with the expansion.
#' h <- 1e-4
#' tp <- theta; tp$sigma_L2.1 <- tp$sigma_L2.1 + h
#' tm <- theta; tm$sigma_L2.1 <- tm$sigma_L2.1 - h
#' c(exact = sum(d3[["sigma_log_L1_sigma_log_L1_sigma_L2.1"]]),
#'   stencil = (sum(distrib_hessian(d, y, tp)[["sigma_log_L1_sigma_log_L1"]]) -
#'              sum(distrib_hessian(d, y, tm)[["sigma_log_L1_sigma_log_L1"]])) /
#'             (2 * h))
#'
#' @keywords internal
mvg_higher <- function(distrib, y, theta, order, expected = FALSE) {
  y <- as_mv_matrix(distrib, y)
  n <- nrow(y)
  p <- distrib@n_dim
  pt <- mvg_ptensors(mvg_pieces(distrib, theta), order,
                     inverted = distrib@inverted)
  pc <- pt$pc
  s <- pc$s
  r <- sweep(y, 2L, pc$mu, "-")

  ldfun <- switch(order - 2L,
    parameters7::param_d3logdet, parameters7::param_d4logdet
  )
  ld <- ldfun(s, pc$eta)
  ldkey <- vapply(
    parameters7::param_tuple_indices(s, order),
    function(t) paste(sort(t), collapse = ","), character(1)
  )

  idx <- deriv_indices(distrib@params, order)
  nms <- deriv_names(distrib@params, order)
  out <- vector("list", length(idx))
  names(out) <- nms

  for (i in seq_along(idx)) {
    t <- idx[[i]]
    is_mu <- t <= p
    n_mu <- sum(is_mu)
    if (n_mu >= 3L) {
      # the quadratic form is quadratic in mu, so a third mean derivative,
      # whatever else the tuple carries, is zero
      out[[i]] <- rep(0, n)
      next
    }
    et <- t[!is_mu] - p
    if (n_mu == 0L) {
      pt_mat <- pt$get(et)
      # E[r' P_t r] = tr(P_t Sigma), the response entering quadratically
      quad <- if (expected) {
        rep(sum(pt_mat * pc$sigma), n)
      } else {
        rowSums((r %*% pt_mat) * r)
      }
      ldc <- ld[[match(paste(sort(et), collapse = ","), ldkey)]]
      out[[i]] <- rep(pt$sign_ld * ldc, n) - 0.5 * quad
    } else if (n_mu == 1L && expected) {
      # linear in the residual, whose expectation is zero
      out[[i]] <- rep(0, n)
    } else if (n_mu == 1L) {
      iu <- t[is_mu]
      vecs <- r %*% pt$get(et)
      out[[i]] <- vecs[, iu]
    } else {
      ij <- t[is_mu]
      out[[i]] <- rep(-pt$get(et)[ij[1L], ij[2L]], n)
    }
  }
  out
}


#' @title Multivariate Gaussian Third Derivatives
#' @name distrib_deriv3.MvGaussianDistrib
#'
#' @description
#' Computes every third derivative of the log-density in the parameters, in
#' closed form. Writing \eqn{r = y - \mu} and \eqn{P_t} for the derivative of
#' the precision over the matrix indices of the tuple: three mean indices give
#' exactly zero, the log-density being quadratic in \eqn{\mu}; two give
#' \eqn{-P_t[i, j]}; one gives \eqn{(r P_t)_i}; none gives
#' \eqn{\pm\tfrac{1}{2}\,\partial^3\log\lvert M\rvert - \tfrac{1}{2} r^\top P_t r}.
#' The arrays \eqn{P_t} come from [mvg_ptensors()].
#'
#' @details
#' With `expected = TRUE` the expectation is exact: the response enters through
#' \eqn{r}, whose expectation is zero, and through \eqn{r^\top P_t r}, whose
#' expectation is \eqn{\operatorname{tr}(P_t\Sigma)}. `approx` and `nsim` are
#' not read.
#'
#' @param distrib An [MvGaussianDistrib] object, from [mvgaussian1_distrib()].
#' @param y An \eqn{n \times p} numeric matrix of observations. With
#'   `expected = TRUE` only its row count is used.
#' @param theta A named list of parameters, each component a single number.
#' @param expected Logical of length 1. When `TRUE` the expectation of each
#'   component is returned. Defaults to `FALSE`.
#' @param scale One of `"parameter"` (the default) or `"link"`. The method
#'   returns the parameter scale and the generic applies the link scale to its
#'   result. Every link of this family is the identity, so the two coincide.
#' @param approx Not read: the expectation is exact.
#' @param nsim Not read, for the same reason.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A named list of numeric vectors of length \eqn{n}, keyed and ordered
#'   as `deriv_names(distrib@params, 3)`. With `expected = TRUE` every vector
#'   is constant.
#'
#' @section Notation:
#' \eqn{\mu} is the mean, \eqn{M} the matrix that the parametrization carries,
#' \eqn{\Sigma} the covariance, \eqn{\eta} the free vector, \eqn{r = y - \mu}
#' the centered response and \eqn{P_t} the derivative of the precision over the
#' multiset \eqn{t}.
#'
#' @seealso [distrib_deriv4.MvGaussianDistrib()] for the next order,
#'   [distrib_hessian.MvGaussianDistrib()] for the second, [mvg_higher()] for
#'   the shared engine, and [distrib_deriv3()] for the generic.
#'
#' @examples
#' d <- mvgaussian1_distrib(2)
#' theta <- list(mu1 = 0.5, mu2 = -0.3, sigma_log_L1 = 0.1,
#'               sigma_log_L2 = -0.2, sigma_L2.1 = 0.4)
#' set.seed(1)
#' y <- distrib_rng(d, 4, theta)
#'
#' d3 <- distrib_deriv3(d, y, theta)
#' length(d3)
#'
#' # Any component naming three means is exactly zero.
#' d3[["mu1_mu1_mu2"]]
#'
#' # Against one stencil on the analytic Hessian.
#' h <- 1e-4
#' tp <- theta; tp$sigma_L2.1 <- tp$sigma_L2.1 + h
#' tm <- theta; tm$sigma_L2.1 <- tm$sigma_L2.1 - h
#' c(exact = sum(d3[["mu1_mu2_sigma_L2.1"]]),
#'   stencil = (sum(distrib_hessian(d, y, tp)[["mu1_mu2"]]) -
#'              sum(distrib_hessian(d, y, tm)[["mu1_mu2"]])) / (2 * h))
#'
#' # The expectation in closed form, against an average over a large sample.
#' e3 <- distrib_deriv3(d, y, theta, expected = TRUE)
#' set.seed(9)
#' big <- distrib_deriv3(d, distrib_rng(d, 50000, theta), theta)
#' c(closed = e3[["sigma_log_L1_sigma_log_L1_sigma_L2.1"]][1],
#'   sampled = mean(big[["sigma_log_L1_sigma_log_L1_sigma_L2.1"]]))
#'
#' @keywords internal
S7::method(distrib_deriv3, MvGaussianDistrib) <- function(distrib, y, theta, expected = FALSE, scale = c("parameter", "link"), approx = c("integrate", "bartlett", "mc", "opg"), nsim = 10000, ...) {
  mvg_higher(distrib, y, theta, 3L, expected = expected)
}

#' @title Multivariate Gaussian Fourth Derivatives
#' @name distrib_deriv4.MvGaussianDistrib
#'
#' @description
#' Computes every fourth derivative of the log-density in the parameters, in
#' closed form, by the algebra of [distrib_deriv3.MvGaussianDistrib()] one
#' order up: the tuple is split into mean indices and matrix indices, three or
#' more mean indices give exactly zero, and the rest are read off the
#' derivative of the precision \eqn{P_t} from [mvg_ptensors()] and the fourth
#' derivative of the log-determinant.
#'
#' @details
#' With `expected = TRUE` the expectation is exact, for the reason given at
#' order three. `approx` and `nsim` are not read.
#'
#' @param distrib An [MvGaussianDistrib] object, from [mvgaussian1_distrib()].
#' @param y An \eqn{n \times p} numeric matrix of observations. With
#'   `expected = TRUE` only its row count is used.
#' @param theta A named list of parameters, each component a single number.
#' @param expected Logical of length 1. When `TRUE` the expectation of each
#'   component is returned. Defaults to `FALSE`.
#' @param scale One of `"parameter"` (the default) or `"link"`. The method
#'   returns the parameter scale and the generic applies the link scale to its
#'   result. Every link of this family is the identity, so the two coincide.
#' @param approx Not read: the expectation is exact.
#' @param nsim Not read, for the same reason.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A named list of numeric vectors of length \eqn{n}, keyed and ordered
#'   as `deriv_names(distrib@params, 4)`. With `expected = TRUE` every vector
#'   is constant.
#'
#' @section Notation:
#' \eqn{\mu} is the mean, \eqn{M} the matrix that the parametrization carries,
#' \eqn{\eta} its free vector, \eqn{r = y - \mu} the centered response and
#' \eqn{P_t} the derivative of the precision over the multiset \eqn{t}.
#'
#' @seealso [distrib_deriv3.MvGaussianDistrib()] for the order below,
#'   [mvg_higher()] for the shared engine, and [distrib_deriv4()] for the
#'   generic.
#'
#' @examples
#' d <- mvgaussian1_distrib(2)
#' theta <- list(mu1 = 0.5, mu2 = -0.3, sigma_log_L1 = 0.1,
#'               sigma_log_L2 = -0.2, sigma_L2.1 = 0.4)
#' set.seed(1)
#' y <- distrib_rng(d, 4, theta)
#'
#' d4 <- distrib_deriv4(d, y, theta)
#' length(d4)
#'
#' # Four mean indices vanish, as three already did.
#' d4[["mu1_mu1_mu2_mu2"]]
#'
#' # Against one stencil on the analytic third order.
#' h <- 1e-4
#' tp <- theta; tp$sigma_L2.1 <- tp$sigma_L2.1 + h
#' tm <- theta; tm$sigma_L2.1 <- tm$sigma_L2.1 - h
#' c(exact = sum(d4[["mu1_mu2_sigma_L2.1_sigma_L2.1"]]),
#'   stencil = (sum(distrib_deriv3(d, y, tp)[["mu1_mu2_sigma_L2.1"]]) -
#'              sum(distrib_deriv3(d, y, tm)[["mu1_mu2_sigma_L2.1"]])) / (2 * h))
#'
#' @keywords internal
S7::method(distrib_deriv4, MvGaussianDistrib) <- function(distrib, y, theta, expected = FALSE, scale = c("parameter", "link"), approx = c("integrate", "bartlett", "mc", "opg"), nsim = 10000, ...) {
  mvg_higher(distrib, y, theta, 4L, expected = expected)
}
