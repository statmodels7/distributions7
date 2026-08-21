#' @include multivariate.R cross_derivatives.R cross2_derivatives.R cross_theta2_derivatives.R
NULL

#' Multivariate Gaussian Distribution
#'
#' @description
#' The S7 class of multivariate gaussian distributions, parametrized by a mean
#' vector and by a \pkg{parameters7} structure for the covariance or the
#' precision. Constructed by \code{\link{mvgaussian_distrib}}.
#'
#' @inheritParams multivariate_distrib
#' @param param The \pkg{parameters7} structure carrying the matrix.
#' @param inverted Whether the matrix parameter parametrizes the precision rather than
#'   the covariance.
#'
#' @return An object of class \code{MvGaussianDistrib}.
#'
#' @seealso \code{\link{mvgaussian_distrib}}
#'
#' @examples
#' S7::S7_inherits(mvgaussian_distrib(2), MvGaussianDistrib)
#'
#' @export
MvGaussianDistrib <- S7::new_class("MvGaussianDistrib",
  parent = multivariate_distrib,
  properties = list(
    param = parameters7::parameter,
    inverted = S7::class_logical
  )
)


#' Construct a Multivariate Gaussian Distribution
#'
#' @description
#' The gaussian distribution on \eqn{\mathbb{R}^p}, with the mean a vector of
#' \eqn{p} free parameters and the matrix carried by a structure from
#' \pkg{parameters7}.
#'
#' @details
#' Exactly one of \code{sigma} and \code{omega} may be given, and
#' the name of the argument decides which side of the model the matrix parameter
#' parametrizes: the covariance in the first case, the precision in the second.
#' One constructor returns one of two behaviors, in the manner of
#' \code{\link{truncated}}, which chooses between its continuous and discrete
#' classes from the arguments it is handed.
#'
#' The precision form is the cheaper one and is preferable where the
#' modeling allows it. Written in \eqn{\Omega}, the log-density, the score and
#' the Hessian are multiplications, and the first term of the score is the
#' parameter's own \code{param_dlogdet()}; written in \eqn{\Sigma} the same
#' quantities need a solve at every step.
#'
#' \strong{Parameters.} The mean contributes \code{mu1}, ..., \code{mup}, and
#' the matrix parameter contributes its free values prefixed by the matrix they
#' describe: \code{sigma_} for a covariance and \code{omega_} for a precision.
#' A two-dimensional gaussian on an unstructured covariance therefore has five
#' parameters, \code{mu1}, \code{mu2}, \code{sigma_log_L1},
#' \code{sigma_log_L2} and \code{sigma_L2.1}, while the same structure on the
#' precision gives \code{omega_log_L1} and the rest. The prefix is what
#' distinguishes the two models in a printed table, since the name of a free
#' value says how the matrix is built and not which matrix it is.
#'
#' All of the parameters are unconstrained, and their links are therefore the
#' identity: the constraint that makes the matrix positive definite lives
#' inside the matrix parameter, which is why it needs no link to express it. A
#' practical consequence is that the parameter scale and the link scale
#' coincide here, so \code{scale = "link"} changes nothing.
#'
#' \strong{Reading a fit.} The free values are coordinates, not quantities
#' anybody reads. \code{\link{mv_summary}} carries the fit's variance matrix
#' onto the standard deviations and correlations by the delta method, and
#' \code{print()} shows them; a precision parametrization also reports the
#' conditional variances and the partial correlations, which are what it
#' describes directly. The conditional variance is
#' \eqn{1/\Omega_{jj} = \mathrm{Var}(Y_j \mid Y_{-j})}, and its ratio to the
#' marginal variance is \eqn{1 - R_j^2} for the regression of that coordinate
#' on all the others.
#'
#' \strong{Rank.} A rank-deficient structure is rejected. A singular covariance
#' gives a law supported on a subspace, with no density against Lebesgue
#' measure, and a singular precision gives a quadratic form that is flat along
#' its null space and does not normalize. The two are different failures and
#' both are failures; a structure of that kind is a legitimate penalty and not
#' a legitimate density.
#'
#' \strong{The response} is an \eqn{n \times p} matrix, one row per
#' observation. A plain vector of length \eqn{p} is read as a single
#' observation.
#'
#' @section The distribution:
#' \deqn{f(y) = (2\pi)^{-p/2}\lvert \Sigma \rvert^{-1/2}\exp\!\left\{-\tfrac{1}{2}(y-\mu)'\Sigma^{-1}(y-\mu)\right\}}
#' on \eqn{y \in \mathbb{R}^{p}}, with
#'
#' \deqn{\mathbb{E}[Y] = \mu, \qquad \operatorname{Var}(Y) = \Sigma.}
#'
#' The matrix parameter is carried by \pkg{parameters7}, either as the
#' covariance \eqn{\Sigma} or as the precision
#' \eqn{\Omega = \Sigma^{-1}}, and its free values are flattened into
#' ordinary scalar parameters with identity links.
#'
#' @param n_dim The dimension \eqn{p}.
#' @param sigma A \pkg{parameters7} structure for the covariance.
#'   Defaults to \code{parameters7::log_cholesky(n_dim)} when neither structure
#'   is given.
#' @param omega A \pkg{parameters7} structure for the precision.
#'
#' @return An object of class \code{\link{MvGaussianDistrib}}.
#'
#' @seealso \code{\link{gaussian1_distrib}}, \code{\link{fit_distrib}},
#'   \code{\link[parameters7]{log_cholesky}}
#'
#' @examples
#' d <- mvgaussian_distrib(2)
#' d
#'
#' theta <- list(mu1 = 0, mu2 = 0, sigma_log_L1 = 0, sigma_log_L2 = 0, sigma_L2.1 = 0.5)
#' y <- rbind(c(0, 0), c(1, -1))
#' distrib_pdf(d, y, theta, log = TRUE)
#'
#' # the covariance the free values describe
#' mv_sigma(d, theta)
#'
#' # a diagonal covariance: two variances instead of three free values
#' mvgaussian_distrib(2, sigma = parameters7::diagonal_matrix(2))@params
#'
#' # or the precision, which is the cheaper parametrization
#' mvgaussian_distrib(2, omega = parameters7::log_cholesky(2))@params
#'
#' @export
mvgaussian_distrib <- function(n_dim, sigma = NULL, omega = NULL) {
  if (!is.numeric(n_dim) || length(n_dim) != 1L || !is.finite(n_dim) ||
    n_dim < 1 || n_dim != round(n_dim)) {
    stop("'n_dim' must be a single positive integer.", call. = FALSE)
  }
  p <- as.integer(n_dim)

  if (!is.null(sigma) && !is.null(omega)) {
    stop(paste0(
      "Give at most one of 'sigma' and 'omega'. They name the\n",
      "  two sides of the same model, and a distribution parametrized by both\n",
      "  would be over-determined."
    ), call. = FALSE)
  }
  inverted <- !is.null(omega)
  s <- if (inverted) omega else sigma
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

  MvGaussianDistrib(
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


#' The Pieces a Multivariate Gaussian Evaluates From
#'
#' @description
#' Assembles, once per call, the mean, the covariance, its inverse and its
#' log-determinant from a flat parameter vector, together with the matrix parameter's
#' derivative matrices when they are asked for.
#'
#' @details
#' Whichever side the matrix parameter parametrizes, the arithmetic below is written
#' in the covariance, so a precision structure is inverted once here rather
#' than at every use. The log-determinant follows the matrix parameter's own, with its
#' sign flipped for a precision, which is the one place the two forms differ in
#' anything but cost.
#'
#' @param distrib A \code{\link{MvGaussianDistrib}} object.
#' @param theta A named list of parameters, already aligned.
#' @param derivs Whether the matrix parameter's first derivative matrices are needed.
#' @param derivs2 Whether its second derivatives are needed as well.
#'
#' @return A list with \code{mu}, \code{sigma}, \code{sigma_inv},
#'   \code{logdet}, \code{eta}, and optionally \code{a} and \code{a2}, the
#'   derivatives of the covariance with respect to the free values.
#'
#' @keywords internal
mvg_pieces <- function(distrib, theta, derivs = FALSE, derivs2 = FALSE) {
  p <- distrib@n_dim
  s <- distrib@param
  v <- mv_flat_theta(distrib, theta)
  mu <- v[seq_len(p)]
  eta <- v[p + seq_len(s@n_free)]

  m <- parameters7::param_value(s, eta)
  ld <- parameters7::param_logdet(s, eta)

  if (distrib@inverted) {
    omega <- m
    sigma <- parameters7::param_solve(s, eta)
    sigma_inv <- omega
    logdet <- -ld
  } else {
    sigma <- m
    sigma_inv <- parameters7::param_solve(s, eta)
    logdet <- ld
  }

  out <- list(
    mu = unname(mu), sigma = unname(sigma), sigma_inv = unname(sigma_inv),
    logdet = logdet, eta = eta, p = p, s = s
  )

  if (derivs || derivs2) {
    d <- parameters7::param_d1(s, eta)
    if (distrib@inverted) {
      # dSigma/deta = -Sigma (dOmega/deta) Sigma, the derivative of an inverse.
      d <- lapply(d, function(ak) -(sigma %*% ak %*% sigma))
    }
    out$a <- lapply(d, unname)
  }
  if (derivs2) {
    d2 <- parameters7::param_d2(s, eta)
    if (distrib@inverted) {
      idx <- parameters7::param_tuple_indices(s)
      om1 <- parameters7::param_d1(s, eta)
      d2 <- lapply(seq_along(idx), function(i) {
        k <- idx[[i]][1L]
        l <- idx[[i]][2L]
        # Differentiating -S B_k S once more, with dS/deta_l = -S B_l S.
        sigma %*% (om1[[l]] %*% sigma %*% om1[[k]] +
          om1[[k]] %*% sigma %*% om1[[l]] - d2[[i]]) %*% sigma
      })
      names(d2) <- parameters7::param_tuple_names(s)
    }
    out$a2 <- lapply(d2, unname)
  }
  out
}


#' @title Mean of a Multivariate Gaussian
#' @name mv_location.MvGaussianDistrib
#' @description The first \eqn{p} parameters, which are the mean vector.
#' @param distrib A \code{\link{MvGaussianDistrib}} object.
#' @param theta A named list of parameters.
#' @return A named numeric vector of length \eqn{p}.
#' @keywords internal
S7::method(mv_location, MvGaussianDistrib) <- mv_leading_location

#' @title The Covariance a Multivariate Gaussian Carries
#' @name mv_sigma.MvGaussianDistrib
#' @description
#' The covariance, assembled from the matrix parameter and inverted first when the
#' structure parametrizes the precision.
#' @param distrib A \code{\link{MvGaussianDistrib}} object.
#' @param theta A named list of parameters.
#' @return A \eqn{p \times p} numeric matrix.
#' @keywords internal
S7::method(mv_sigma, MvGaussianDistrib) <- function(distrib, theta) {
  pc <- mvg_pieces(distrib, theta)
  nm <- paste0("v", seq_len(distrib@n_dim))
  dimnames(pc$sigma) <- list(nm, nm)
  pc$sigma
}


#' Residuals and Whitened Residuals
#'
#' @description
#' The centered response and its image under the inverse covariance, which are
#' what every derivative of a multivariate gaussian is written in.
#'
#' @param y An \eqn{n \times p} matrix.
#' @param pc The result of \code{\link{mvg_pieces}}.
#'
#' @return A list with \code{r}, the residuals, and \code{w}, the rows of
#'   \eqn{R \Sigma^{-1}}.
#'
#' @keywords internal
mvg_residuals <- function(y, pc) {
  r <- sweep(y, 2L, pc$mu, "-")
  list(r = r, w = r %*% pc$sigma_inv)
}


#' @title Multivariate Gaussian Density
#' @name distrib_pdf.MvGaussianDistrib
#' @description
#' \deqn{\ell = -\frac{p}{2}\log 2\pi - \frac{1}{2}\log|\Sigma|
#'   - \frac{1}{2}(y-\mu)^\top \Sigma^{-1} (y-\mu),}
#' evaluated row by row. The quadratic form goes through the matrix parameter's own
#' factor rather than through an explicit inverse.
#' @param distrib A \code{\link{MvGaussianDistrib}} object.
#' @param y An \eqn{n \times p} matrix of observations, or a vector of length
#'   \eqn{p} for one observation.
#' @param theta A named list of parameters.
#' @param log Logical; if \code{TRUE}, returns the log-density.
#' @return A numeric vector with one value per observation.
#' @keywords internal
S7::method(distrib_pdf, MvGaussianDistrib) <- function(distrib, y, theta,
                                                       log = FALSE, ...) {
  y <- as_mv_matrix(distrib, y)
  pc <- mvg_pieces(distrib, theta)
  if (!nrow(y)) return(numeric(0))
  res <- mvg_residuals(y, pc)
  q <- rowSums(res$r * res$w)
  ld <- -0.5 * (pc$p * base::log(2 * pi) + pc$logdet + q)
  if (log) ld else exp(ld)
}


#' @title Multivariate Gaussian Generator
#' @name distrib_rng.MvGaussianDistrib
#' @description
#' \eqn{\mu + L z} with \eqn{z} standard normal and \eqn{LL^\top = \Sigma}, the
#' factor taken from the matrix parameter where it parametrizes the covariance and
#' from a factorization of the inverse otherwise.
#' @param distrib A \code{\link{MvGaussianDistrib}} object.
#' @param n The number of observations to draw.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return An \eqn{n \times p} numeric matrix.
#' @keywords internal
S7::method(distrib_rng, MvGaussianDistrib) <- function(distrib, n, theta, ...) {
  pc <- mvg_pieces(distrib, theta)
  p <- pc$p
  l <- t(chol(pc$sigma))
  z <- matrix(stats::rnorm(n * p), nrow = n, ncol = p)
  out <- z %*% t(l)
  out <- sweep(out, 2L, pc$mu, "+")
  colnames(out) <- paste0("v", seq_len(p))
  out
}


#' @title Multivariate Gaussian Score
#' @name distrib_gradient.MvGaussianDistrib
#' @description
#' Closed form. With \eqn{w = \Sigma^{-1}(y - \mu)} and \eqn{A_k} the
#' derivative of \eqn{\Sigma} in the \eqn{k}-th free value of the matrix parameter,
#' \deqn{\frac{\partial \ell}{\partial \mu} = w, \qquad
#'   \frac{\partial \ell}{\partial \eta_k} =
#'   -\frac{1}{2}\frac{\partial \log|\Sigma|}{\partial \eta_k}
#'   + \frac{1}{2} w^\top A_k w.}
#' The first term of the second expression is the matrix parameter's own
#' \code{param_dlogdet()}, so no trace is formed here.
#' @param distrib A \code{\link{MvGaussianDistrib}} object.
#' @param y An \eqn{n \times p} matrix of observations.
#' @param theta A named list of parameters.
#' @param scale Handled by the generic; the two scales coincide here.
#' @param ... Unused.
#' @return A named list with one numeric vector per parameter.
#' @keywords internal
S7::method(distrib_gradient, MvGaussianDistrib) <- function(distrib, y, theta,
                                                            scale = c("parameter", "link"),
                                                            ...) {
  y <- as_mv_matrix(distrib, y)
  pc <- mvg_pieces(distrib, theta, derivs = TRUE)
  n <- nrow(y)
  res <- mvg_residuals(y, pc)
  w <- res$w

  out <- vector("list", distrib@n_params)
  names(out) <- distrib@params
  for (j in seq_len(pc$p)) out[[j]] <- w[, j]

  dld <- parameters7::param_dlogdet(pc$s, pc$eta)
  if (distrib@inverted) dld <- -dld
  for (k in seq_along(pc$a)) {
    out[[pc$p + k]] <- -0.5 * dld[[k]] +
      0.5 * rowSums((w %*% pc$a[[k]]) * w)
  }
  out
}


#' @title Multivariate Gaussian Observed Hessian
#' @name distrib_hessian.MvGaussianDistrib
#' @description
#' Closed form. With \eqn{w = \Sigma^{-1}(y-\mu)}, \eqn{A_k} and \eqn{A_{kl}}
#' the first and second derivatives of \eqn{\Sigma},
#' \deqn{\ell^{(\mu_a \mu_b)} = -(\Sigma^{-1})_{ab}, \qquad
#'   \ell^{(\mu_a \eta_k)} = -(\Sigma^{-1} A_k w)_a,}
#' \deqn{\ell^{(\eta_k \eta_l)} =
#'   -\tfrac{1}{2}\frac{\partial^2 \log|\Sigma|}{\partial\eta_k \partial\eta_l}
#'   + \tfrac{1}{2} w^\top A_{kl} w - w^\top A_l \Sigma^{-1} A_k w.}
#' @param distrib A \code{\link{MvGaussianDistrib}} object.
#' @param y An \eqn{n \times p} matrix of observations.
#' @param theta A named list of parameters.
#' @param scale Handled by the generic; the two scales coincide here.
#' @param ... Unused.
#' @return A named list keyed as \code{\link{hess_names}(distrib@params)}.
#' @keywords internal
S7::method(distrib_hessian, MvGaussianDistrib) <- function(distrib, y, theta,
                                                           scale = c("parameter", "link"),
                                                           ...) {
  y <- as_mv_matrix(distrib, y)
  pc <- mvg_pieces(distrib, theta, derivs2 = TRUE)
  n <- nrow(y)
  p <- pc$p
  w <- mvg_residuals(y, pc)$w
  si <- pc$sigma_inv

  d2ld <- parameters7::param_d2logdet(pc$s, pc$eta)
  if (distrib@inverted) d2ld <- -d2ld
  # Sigma^{-1} A_k, formed once: it appears in the mixed block and again in the
  # matrix block.
  sa <- lapply(pc$a, function(ak) si %*% ak)
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
      -(w %*% t(sa[[b - p]]))[, a]
    } else if (b <= p) {
      -(w %*% t(sa[[a - p]]))[, b]
    } else {
      k <- a - p
      l <- b - p
      idx <- spair[[paste(min(k, l), max(k, l), sep = ":")]]
      -0.5 * d2ld[[idx]] + 0.5 * rowSums((w %*% pc$a2[[idx]]) * w) -
        rowSums((w %*% pc$a[[l]] %*% si %*% pc$a[[k]]) * w)
    }
  }
  out
}


#' Index Pairs Behind the Hessian Keys of a Multivariate Distribution
#'
#' @description
#' The positions in \code{distrib@params} of each unordered pair, in the order
#' \code{\link{hess_names}} uses.
#'
#' @details
#' \code{\link{hess_pairs}} returns pairs of parameter NAMES, which is what a
#' univariate method wants when it looks a component up. A multivariate method
#' needs the positions instead, to tell a mean component from a matrix one, so
#' the names are matched back rather than the enumeration being written a
#' second time.
#'
#' @param distrib A \code{\link{multivariate_distrib}} object.
#'
#' @return A list of integer vectors of length 2.
#'
#' @keywords internal
mv_hess_indices <- function(distrib) {
  params <- distrib@params
  lapply(hess_pairs(params), function(pr) match(pr, params))
}


#' Where Each Pair of Free Values Sits in a Structure's Second Derivatives
#'
#' @description
#' A lookup from a pair of free-value positions to the position of the
#' corresponding component of \code{param_d2()}.
#'
#' @details
#' Built from \pkg{parameters7}'s own enumeration rather than by taking a
#' component key apart, for the reason that package documents: a free value
#' whose label contains the separator splits into the wrong number of pieces.
#'
#' @param s A \pkg{parameters7} structure.
#'
#' @return A named list of integers, keyed \code{"k:l"} with \eqn{k \le l}.
#'
#' @keywords internal
param_pair_lookup <- function(s) {
  idx <- parameters7::param_tuple_indices(s)
  keys <- vapply(idx, function(kl) {
    paste(min(kl), max(kl), sep = ":")
  }, character(1))
  stats::setNames(as.list(seq_along(idx)), keys)
}


#' @title Multivariate Gaussian Expected Hessian
#' @name distrib_expected_hessian.MvGaussianDistrib
#' @description
#' Closed form, and simpler than the observed one.
#' \deqn{\mathbb{E}[\ell^{(\mu_a \mu_b)}] = -(\Sigma^{-1})_{ab}, \qquad
#'   \mathbb{E}[\ell^{(\mu_a \eta_k)}] = 0, \qquad
#'   \mathbb{E}[\ell^{(\eta_k \eta_l)}] =
#'   -\tfrac{1}{2}\mathrm{tr}(\Sigma^{-1} A_k \Sigma^{-1} A_l).}
#' @details
#' The mixed block vanishes because \eqn{\mathbb{E}[w] = 0}, which is the
#' orthogonality of the mean and the covariance parameters that makes Fisher
#' scoring on this family so well behaved. The matrix block needs no second
#' derivative at all: the terms in \eqn{A_{kl}} cancel between the
#' log-determinant and the quadratic form.
#' @param distrib A \code{\link{MvGaussianDistrib}} object.
#' @param y An \eqn{n \times p} matrix of observations.
#' @param theta A named list of parameters.
#' @param scale Handled by the generic; the two scales coincide here.
#' @param approx Ignored: the expectation is exact.
#' @param nsim Ignored.
#' @param ... Unused.
#' @return A named list keyed as \code{\link{hess_names}(distrib@params)}.
#' @keywords internal
S7::method(distrib_expected_hessian, MvGaussianDistrib) <- function(
    distrib, y, theta, scale = c("parameter", "link"),
    approx = c("bartlett", "integrate", "mc", "opg"), nsim = 10000, ...) {
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
#' @description
#' \eqn{\partial \ell / \partial y = -\Sigma^{-1}(y - \mu)}, one row per
#' observation. The shape differs from the univariate case, where the
#' derivative in a scalar response is a vector: here it is an
#' \eqn{n \times p} matrix.
#' @param distrib A \code{\link{MvGaussianDistrib}} object.
#' @param y An \eqn{n \times p} matrix of observations.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return An \eqn{n \times p} numeric matrix.
#' @keywords internal
S7::method(distrib_grad_y, MvGaussianDistrib) <- function(distrib, y, theta, ...) {
  y <- as_mv_matrix(distrib, y)
  pc <- mvg_pieces(distrib, theta)
  -mvg_residuals(y, pc)$w
}


#' @title Multivariate Gaussian Response Hessian
#' @name distrib_hess_y.MvGaussianDistrib
#' @description
#' \eqn{\partial^2 \ell / \partial y \partial y^\top = -\Sigma^{-1}}, the same
#' matrix at every observation, so it is returned once rather than repeated.
#' @param distrib A \code{\link{MvGaussianDistrib}} object.
#' @param y An \eqn{n \times p} matrix of observations.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return A \eqn{p \times p} numeric matrix.
#' @keywords internal
S7::method(distrib_hess_y, MvGaussianDistrib) <- function(distrib, y, theta, ...) {
  pc <- mvg_pieces(distrib, theta)
  -pc$sigma_inv
}


#' @title Multivariate Gaussian Mixed Response-Parameter Derivatives
#' @name distrib_cross_y.MvGaussianDistrib
#' @description
#' \eqn{\partial^2 \ell / \partial y \partial \theta_k}, one \eqn{n \times p}
#' matrix per parameter. With \eqn{w = \Sigma^{-1}(y - \mu)} the response
#' gradient is \eqn{-w}, so differentiating it in the mean and in the free
#' values of the matrix parameter gives
#' \deqn{\frac{\partial^2 \ell}{\partial y \partial \mu_j} = \Sigma^{-1}e_j,
#'   \qquad
#'   \frac{\partial^2 \ell}{\partial y \partial \eta_k} = \Sigma^{-1}A_k w,}
#' with \eqn{A_k = \partial\Sigma/\partial\eta_k}. The mean block is the same
#' at every observation, the matrix block is not.
#' @details
#' The shape is the one a consumer needs: a penalty whose prior is this family
#' reads the block of \eqn{\partial^2\rho/\partial\beta\,\partial\theta_k} for
#' one hyperparameter at a time, and the coefficients of one group are the row
#' of \eqn{y} the density is read at.
#'
#' The link scale is the parameter scale here: the mean components carry the
#' identity link and so do the matrix parameter's free values, which are
#' unconstrained by construction.
#' @param distrib A \code{\link{MvGaussianDistrib}} object.
#' @param y An \eqn{n \times p} matrix of observations.
#' @param theta A named list of parameters.
#' @param scale Handled by the generic.
#' @param ... Unused.
#' @return A named list, one \eqn{n \times p} matrix per parameter.
#' @keywords internal
S7::method(distrib_cross_y, MvGaussianDistrib) <-
  function(distrib, y, theta, scale = c("parameter", "link"), ...) {
    y <- as_mv_matrix(distrib, y)
    p <- distrib@n_dim
    pc <- mvg_pieces(distrib, theta, derivs = TRUE)
    w <- mvg_residuals(y, pc)$w
    si <- pc$sigma_inv
    n <- nrow(y)

    mu_part <- lapply(seq_len(p), function(j) {
      matrix(si[j, ], nrow = n, ncol = p, byrow = TRUE)
    })
    # rows: (Sigma^-1 A_k w_i)' = w_i' A_k Sigma^-1, both matrices symmetric
    eta_part <- lapply(pc$a, function(ak) w %*% ak %*% si)

    stats::setNames(c(mu_part, eta_part), distrib@params)
  }


#' @title Multivariate Gaussian Higher Mixed Response Derivatives
#' @name distrib_cross2_y.MvGaussianDistrib
#' @description
#' The three derivatives a marginal criterion reads when this family is a
#' prior. Writing \eqn{B_k = \Sigma^{-1}A_k\Sigma^{-1}} with
#' \eqn{A_k = \partial\Sigma/\partial\eta_k},
#' \deqn{\frac{\partial^3\ell}{\partial y\,\partial y^\top\partial\eta_k}
#'     = B_k, \qquad
#'   \frac{\partial^4\ell}{\partial y\,\partial y^\top
#'     \partial\eta_k\partial\eta_l}
#'     = \Sigma^{-1}A_{kl}\Sigma^{-1}
#'       - \Sigma^{-1}\!\left(A_l\Sigma^{-1}A_k
#'         + A_k\Sigma^{-1}A_l\right)\!\Sigma^{-1},}
#' \deqn{\frac{\partial^3\ell}{\partial y\,\partial\mu_j\partial\eta_k}
#'     = -B_k e_j, \qquad
#'   \frac{\partial^3\ell}{\partial y\,\partial\eta_k\partial\eta_l}
#'     = \Sigma^{-1}A_{kl}w
#'       - \Sigma^{-1}\!\left(A_l\Sigma^{-1}A_k
#'         + A_k\Sigma^{-1}A_l\right)\!w.}
#' @details
#' The response Hessian is \eqn{-\Sigma^{-1}}, which does not depend on the
#' observation and does not depend on the mean at all, so every component of
#' the first two involving a mean is exactly zero and the rest are one matrix
#' rather than one per row. Only the third carries the observation, through
#' \eqn{w}.
#' @param distrib A \code{\link{MvGaussianDistrib}} object.
#' @param y An \eqn{n \times p} matrix of observations.
#' @param theta A named list of parameters.
#' @param scale Handled by the generic.
#' @param ... Unused.
#' @return A named list, keyed by parameter for \code{distrib_cross2_y} and by
#'   parameter pair for the other two.
#' @keywords internal
S7::method(distrib_cross2_y, MvGaussianDistrib) <-
  function(distrib, y, theta, scale = c("parameter", "link"), ...) {
    p <- distrib@n_dim
    pc <- mvg_pieces(distrib, theta, derivs = TRUE)
    si <- pc$sigma_inv
    zero <- matrix(0, p, p)
    stats::setNames(
      c(rep(list(zero), p), lapply(pc$a, function(ak) si %*% ak %*% si)),
      distrib@params)
  }

#' @rdname distrib_cross2_y.MvGaussianDistrib
#' @name distrib_hess_y_hess.MvGaussianDistrib
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
      aa <- pc$a[[ka]]
      ab <- pc$a[[kb]]
      si %*% (mvg_a2(pc, ka, kb) - aa %*% si %*% ab - ab %*% si %*% aa) %*% si
    }), hess_names(nm))
  }

#' @rdname distrib_cross2_y.MvGaussianDistrib
#' @name distrib_grad_y_hess.MvGaussianDistrib
#' @keywords internal
S7::method(distrib_grad_y_hess, MvGaussianDistrib) <-
  function(distrib, y, theta, scale = c("parameter", "link"), ...) {
    y <- as_mv_matrix(distrib, y)
    p <- distrib@n_dim
    n <- nrow(y)
    pc <- mvg_pieces(distrib, theta, derivs2 = TRUE)
    si <- pc$sigma_inv
    w <- mvg_residuals(y, pc)$w
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
        # -(Sigma^-1 A_k Sigma^-1) e_j, the same row at every observation
        m <- -(si %*% pc$a[[kk]] %*% si)
        return(matrix(m[, j], nrow = n, ncol = p, byrow = TRUE))
      }
      ka <- a - p
      kb <- b - p
      aa <- pc$a[[ka]]
      ab <- pc$a[[kb]]
      mid <- mvg_a2(pc, ka, kb) - aa %*% si %*% ab - ab %*% si %*% aa
      w %*% mid %*% si
    }), hess_names(nm))
  }

#' The Second Derivative of the Covariance, by Position
#'
#' @description
#' \code{param_d2} keyed by the pair of free values rather than by the string
#' the structure names it with, the key being CONSTRUCTED from the sorted
#' pair and never parsed out of a name.
#'
#' @param pc The result of \code{\link{mvg_pieces}} with \code{derivs2}.
#' @param k,l Positions among the structure's free values.
#'
#' @return A \eqn{p \times p} numeric matrix.
#'
#' @keywords internal
mvg_a2 <- function(pc, k, l) {
  nm <- pc$s@free_names
  ij <- sort(c(k, l))
  pc$a2[[paste(nm[ij], collapse = ":")]]
}

#' @title Mean of a Multivariate Gaussian
#' @name mean.MvGaussianDistrib
#' @description The mean vector, which is a parameter of the family.
#' @param x A \code{\link{MvGaussianDistrib}} object.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return A numeric vector of length \eqn{p}.
#' @keywords internal
S7::method(mean, MvGaussianDistrib) <- function(x, theta, ...) {
  mv_location(x, theta)
}


#' @title Variance of a Multivariate Gaussian
#' @name variance.MvGaussianDistrib
#' @description
#' The covariance matrix, which the matrix parameter carries. The return is a matrix
#' rather than the numeric vector a univariate distribution gives, since that
#' is what the second moment of a vector is.
#' @param x A \code{\link{MvGaussianDistrib}} object.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return A \eqn{p \times p} numeric matrix.
#' @keywords internal
S7::method(variance, MvGaussianDistrib) <- function(x, theta, ...) {
  mv_sigma(x, theta)
}


#' @title Random Parameters for a Multivariate Gaussian
#' @name generate_random_theta.MvGaussianDistrib
#' @description
#' A random mean near the origin and a structure drawn near the identity. The
#' default of the base class would draw every free value from the same wide
#' range, which for a log-Cholesky diagonal spans four orders of magnitude in
#' the resulting variances and gives a starting covariance no fit recovers
#' from.
#' @param distrib A \code{\link{MvGaussianDistrib}} object.
#' @param ... Unused.
#' @return A named list of scalars.
#' @keywords internal
S7::method(generate_random_theta, MvGaussianDistrib) <- function(distrib, ...) {
  p <- distrib@n_dim
  s <- distrib@param
  as.list(stats::setNames(
    c(stats::runif(p, -1, 1), stats::runif(s@n_free, -0.4, 0.4)),
    distrib@params
  ))
}


#' Precision Derivative Tensors of a Multivariate Gaussian
#'
#' @description
#' The precision's derivative tensors in the matrix parameter's free values, orders 1
#' to 4, keyed by index tuple. For a precision parametrization they are the
#' parameter's own derivatives; for a covariance they follow from repeated
#' differentiation of the inverse, so no expanded formula is transcribed and
#' no term can be dropped.
#' It takes the pieces rather than the distribution, because the multivariate
#' Student t needs the same tensors of its scale matrix and there must be one
#' copy of the expansion: its first draft double counted the mixed terms, which
#' only a comparison against a stencil caught.
#' @param pc The pieces, as returned by \code{\link{mvg_pieces}} or
#'   \code{mvt_pieces()}: a list carrying \code{s}, \code{eta} and
#'   \code{sigma_inv}.
#' @param order The highest order wanted.
#' @param inverted Whether the free values parametrize the precision, in which
#'   case the tensors are the parameter's own.
#' @return A list with the accessor \code{get}, the log-determinant sign and
#'   the pieces.
#' @keywords internal
mvg_ptensors <- function(pc, order, inverted = FALSE) {
  s <- pc$s
  eta <- pc$eta
  a <- list(parameters7::param_d1(s, eta))
  if (order >= 2L) a[[2L]] <- parameters7::param_d2(s, eta)
  if (order >= 3L) a[[3L]] <- parameters7::param_d3(s, eta)
  if (order >= 4L) a[[4L]] <- parameters7::param_d4(s, eta)

  key <- function(t) paste(sort(t), collapse = ",")
  amap <- list()
  for (o in seq_len(order)) {
    idx <- if (o == 1L) {
      lapply(seq_len(s@n_free), identity)
    } else {
      parameters7::param_tuple_indices(s, o)
    }
    for (i in seq_along(idx)) {
      amap[[paste0(o, ":", key(idx[[i]]))]] <- unname(a[[o]][[i]])
    }
  }
  aget <- function(t) amap[[paste0(length(t), ":", key(t))]]

  if (inverted) {
    return(list(get = aget, sign_ld = 0.5, pc = pc))
  }

  P <- pc$sigma_inv
  pmap <- new.env(parent = emptyenv())
  pget <- function(t) {
    # the empty multiset is the matrix itself; the gaussian never asks for it,
    # but the Student t does, a partition block of pure mean indices carrying
    # no matrix index at all
    if (!length(t)) return(P)
    k <- key(t)
    got <- pmap[[k]]
    if (!is.null(got)) return(got)
    # P_t = sum over ORDERED partitions (B1, ..., Bq) of the multiset t into
    # nonempty blocks of (-1)^q P A_{B1} P A_{B2} ... A_{Bq} P: the expansion
    # of the derivative of an inverse, checked at order one (-P A_k P) and
    # order two (P A_k P A_l P + P A_l P A_k P - P A_kl P).
    acc <- matrix(0, nrow(P), ncol(P))
    for (blocks in mv_ordered_partitions(seq_along(t))) {
      term <- P
      for (b in blocks) term <- term %*% aget(t[b]) %*% P
      acc <- acc + (-1)^length(blocks) * term
    }
    pmap[[k]] <- acc
    acc
  }
  list(get = pget, sign_ld = -0.5, pc = pc)
}


#' Ordered Partitions of a Set of Positions
#'
#' @description
#' All ordered partitions of a set of positions into nonempty blocks: the set
#' partitions, each in every ordering of its blocks. This is how the
#' derivative of an inverse distributes its differentiations,
#' \eqn{P_t = \sum (-1)^q P A_{B_1} P \cdots A_{B_q} P}.
#' @param pos An integer vector of positions.
#' @return A list of lists of integer vectors.
#' @keywords internal
mv_ordered_partitions <- function(pos) {
  if (length(pos) == 1L) return(list(list(pos)))
  # set partitions by recursive insertion of the last element
  set_parts <- function(v) {
    if (length(v) == 1L) return(list(list(v)))
    prev <- set_parts(v[-length(v)])
    out <- list()
    x <- v[length(v)]
    for (pp in prev) {
      for (j in seq_along(pp)) {
        q <- pp
        q[[j]] <- c(q[[j]], x)
        out[[length(out) + 1L]] <- q
      }
      out[[length(out) + 1L]] <- c(pp, list(x))
    }
    out
  }
  perms <- function(n) {
    if (n == 1L) return(list(1L))
    out <- list()
    for (pr in perms(n - 1L)) {
      for (j in seq_len(n)) {
        out[[length(out) + 1L]] <- append(pr, n, after = j - 1L)
      }
    }
    out
  }
  out <- list()
  for (pp in set_parts(pos)) {
    q <- length(pp)
    for (pr in perms(q)) {
      out[[length(out) + 1L]] <- pp[pr]
    }
  }
  out
}


#' The Closed-Form Higher Derivatives of a Multivariate Gaussian
#'
#' @description
#' Shared engine for the third and fourth derivatives: enumerates the
#' parameter tuples the way \code{\link{deriv_names}} does, splits each into
#' mean and structure indices, and reads the surviving cases off the
#' gaussian's algebra.
#' @param distrib A \code{\link{MvGaussianDistrib}} object.
#' @param y An \eqn{n \times p} matrix of observations.
#' @param theta A named list of parameters.
#' @param order 3 or 4.
#' @return A named list of derivative component vectors.
#' @keywords internal
mvg_higher <- function(distrib, y, theta, order) {
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
      quad <- rowSums((r %*% pt_mat) * r)
      ldc <- ld[[match(paste(sort(et), collapse = ","), ldkey)]]
      out[[i]] <- rep(pt$sign_ld * ldc, n) - 0.5 * quad
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
#' @description
#' Closed form, built on the matrix parameter's own third derivatives from
#' \pkg{parameters7}. A component with three mean indices vanishes, the
#' quadratic form being quadratic; one mean index gives
#' \eqn{(P_{klm} r)_i}; two give \eqn{-P_{kl}[i, j]}; none gives
#' \eqn{\mp\tfrac{1}{2}\,\partial^3 \log|M| - \tfrac{1}{2} r^\top P_{klm} r}.
#' The precision's derivative tensors \eqn{P_t} come directly from the
#' structure under a precision parametrization, and by the Leibniz recursion
#' on \eqn{P_k = -P A_k P} under a covariance one.
#' @param distrib A \code{\link{MvGaussianDistrib}} object.
#' @param y An \eqn{n \times p} matrix of observations.
#' @param theta A named list of parameters.
#' @param expected Logical; the expectation is approximated by sampling, as
#'   for the Hessian.
#' @param approx Strategy label; sampling is the only multivariate route.
#' @param nsim Monte Carlo sample size.
#' @return A named list of third-derivative component vectors.
#' @keywords internal
S7::method(distrib_deriv3, MvGaussianDistrib) <- function(distrib, y, theta, expected = FALSE, scale = c("parameter", "link"), approx = c("integrate", "bartlett", "mc", "opg"), nsim = 10000, ...) {
  if (expected) {
    big <- distrib_rng(distrib, nsim, theta)
    ob <- mvg_higher(distrib, big, theta, 3L)
    n <- n_obs(distrib, y)
    return(lapply(ob, function(v) rep(mean(v), n)))
  }
  mvg_higher(distrib, y, theta, 3L)
}

#' @title Multivariate Gaussian Fourth Derivatives
#' @name distrib_deriv4.MvGaussianDistrib
#' @description Closed form; the same algebra as
#'   \code{\link{distrib_deriv3.MvGaussianDistrib}} one order up.
#' @param distrib A \code{\link{MvGaussianDistrib}} object.
#' @param y An \eqn{n \times p} matrix of observations.
#' @param theta A named list of parameters.
#' @param expected Logical; the expectation is approximated by sampling.
#' @param approx Strategy label; sampling is the only multivariate route.
#' @param nsim Monte Carlo sample size.
#' @return A named list of fourth-derivative component vectors.
#' @keywords internal
S7::method(distrib_deriv4, MvGaussianDistrib) <- function(distrib, y, theta, expected = FALSE, scale = c("parameter", "link"), approx = c("integrate", "bartlett", "mc", "opg"), nsim = 10000, ...) {
  if (expected) {
    big <- distrib_rng(distrib, nsim, theta)
    ob <- mvg_higher(distrib, big, theta, 4L)
    n <- n_obs(distrib, y)
    return(lapply(ob, function(v) rep(mean(v), n)))
  }
  mvg_higher(distrib, y, theta, 4L)
}
