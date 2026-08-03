#' @include multivariate.R
NULL

#' Multivariate Gaussian Distribution
#'
#' @description
#' The S7 class of multivariate gaussian distributions, parametrised by a mean
#' vector and by a \pkg{covstructs7} structure for the covariance or the
#' precision. Constructed by \code{\link{mvgaussian_distrib}}.
#'
#' @inheritParams multivariate_distrib
#' @param struct The \pkg{covstructs7} structure carrying the matrix.
#' @param inverted Whether the structure parametrises the precision rather than
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
    struct = covstructs7::covstruct,
    inverted = S7::class_logical
  )
)


#' Construct a Multivariate Gaussian Distribution
#'
#' @description
#' The gaussian distribution on \eqn{\mathbb{R}^p}, with the mean a vector of
#' \eqn{p} free parameters and the matrix carried by a structure from
#' \pkg{covstructs7}.
#'
#' @details
#' Exactly one of \code{struct_sigma} and \code{struct_omega} may be given, and
#' the name of the argument decides which side of the model the structure
#' parametrises: the covariance in the first case, the precision in the second.
#' One constructor returns one of two behaviours, in the manner of
#' \code{\link{truncated}}, which chooses between its continuous and discrete
#' classes from the arguments it is handed.
#'
#' The precision form is the cheaper one and is worth preferring where the
#' modelling allows it. Written in \eqn{\Omega}, the log-density, the score and
#' the Hessian are multiplications, and the first term of the score is the
#' structure's own \code{struct_dlogdet()}; written in \eqn{\Sigma} the same
#' quantities need a solve at every step.
#'
#' \strong{Parameters.} The mean contributes \code{mu1}, ..., \code{mup} and
#' the structure contributes its own free values under their own names, so a
#' two-dimensional gaussian on an unstructured covariance has five parameters:
#' \code{mu1}, \code{mu2}, \code{log_L1}, \code{log_L2}, \code{L2.1}. All of
#' them are unconstrained, and their links are therefore the identity: the
#' constraint that makes the matrix positive definite lives inside the
#' structure, which is why it needs no link to express it. A consequence worth
#' knowing is that the parameter scale and the link scale coincide here, so
#' \code{scale = "link"} changes nothing.
#'
#' \strong{Rank.} A rank-deficient structure is refused. A singular covariance
#' gives a law supported on a subspace, with no density against Lebesgue
#' measure, and a singular precision gives a quadratic form that is flat along
#' its null space and does not normalise. The two are different failures and
#' both are failures; a structure of that kind is a legitimate penalty and not
#' a legitimate density.
#'
#' \strong{The response} is an \eqn{n \times p} matrix, one row per
#' observation. A plain vector of length \eqn{p} is read as a single
#' observation.
#'
#' @param n_dim The dimension \eqn{p}.
#' @param struct_sigma A \pkg{covstructs7} structure for the covariance.
#'   Defaults to \code{covstructs7::log_cholesky(n_dim)} when neither structure
#'   is given.
#' @param struct_omega A \pkg{covstructs7} structure for the precision.
#'
#' @return An object of class \code{\link{MvGaussianDistrib}}.
#'
#' @seealso \code{\link{gaussian_distrib}}, \code{\link{fit_distrib}},
#'   \code{\link[covstructs7]{log_cholesky}}
#'
#' @examples
#' d <- mvgaussian_distrib(2)
#' d
#'
#' theta <- list(mu1 = 0, mu2 = 0, log_L1 = 0, log_L2 = 0, L2.1 = 0.5)
#' y <- rbind(c(0, 0), c(1, -1))
#' distrib_pdf(d, y, theta, log = TRUE)
#'
#' # the covariance the free values describe
#' mv_sigma(d, theta)
#'
#' # a diagonal covariance: two variances instead of three free values
#' mvgaussian_distrib(2, struct_sigma = covstructs7::diag_struct(2))@params
#'
#' # or the precision, which is the cheaper parametrisation
#' mvgaussian_distrib(2, struct_omega = covstructs7::log_cholesky(2))@params
#'
#' @export
mvgaussian_distrib <- function(n_dim, struct_sigma = NULL, struct_omega = NULL) {
  if (!is.numeric(n_dim) || length(n_dim) != 1L || !is.finite(n_dim) ||
    n_dim < 1 || n_dim != round(n_dim)) {
    stop("'n_dim' must be a single positive integer.", call. = FALSE)
  }
  p <- as.integer(n_dim)

  if (!is.null(struct_sigma) && !is.null(struct_omega)) {
    stop(paste0(
      "Give at most one of 'struct_sigma' and 'struct_omega'. They name the\n",
      "  two sides of the same model, and a distribution parametrised by both\n",
      "  would be over-determined."
    ), call. = FALSE)
  }
  inverted <- !is.null(struct_omega)
  s <- if (inverted) struct_omega else struct_sigma
  if (is.null(s)) s <- covstructs7::log_cholesky(p)

  if (!S7::S7_inherits(s, covstructs7::covstruct)) {
    stop("The structure must be a covstructs7 'covstruct' object.", call. = FALSE)
  }
  if (s@dimension != p) {
    stop(sprintf(
      "The structure has dimension %d but the distribution has dimension %d.",
      s@dimension, p
    ), call. = FALSE)
  }
  if (s@rank < s@dimension) {
    stop(sprintf(paste0(
      "The structure is rank deficient (%d of %d), so it does not describe a\n",
      "  density. A singular covariance is supported on a subspace and a\n",
      "  singular precision does not normalise; either way the law has no\n",
      "  density against Lebesgue measure. Such a structure is a penalty, not\n",
      "  a distribution."
    ), s@rank, s@dimension), call. = FALSE)
  }

  mu_names <- paste0("mu", seq_len(p))
  clash <- intersect(mu_names, s@free_names)
  if (length(clash)) {
    stop(sprintf(paste0(
      "The structure's free value '%s' has the same name as a mean component.\n",
      "  Parameter names must be unique, since every derivative component is\n",
      "  keyed by them."
    ), clash[1L]), call. = FALSE)
  }

  params <- c(mu_names, s@free_names)
  n_par <- length(params)

  MvGaussianDistrib(
    # No spaces inside the brackets: print() capitalises the first letter after
    # every space in a distribution's name, and "Covariance Log_cholesky" is
    # not what the structure is called. Same convention as truncated().
    distrib_name = sprintf(
      "multivariate gaussian [%dd, %s=%s]", p,
      if (inverted) "omega" else "sigma", s@struct_name
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
    # structure's free values are free by construction. The constraint that
    # makes the matrix positive definite is inside the structure, which is
    # exactly why it does not need a link to carry it.
    params_bounds = stats::setNames(
      rep(list(c(-Inf, Inf)), n_par), params
    ),
    link_params = stats::setNames(
      rep(list(linkfunctions7::identity_link()), n_par), params
    ),
    struct = s,
    inverted = inverted
  )
}


#' The Pieces a Multivariate Gaussian Evaluates From
#'
#' @description
#' Assembles, once per call, the mean, the covariance, its inverse and its
#' log-determinant from a flat parameter vector, together with the structure's
#' derivative matrices when they are asked for.
#'
#' @details
#' Whichever side the structure parametrises, the arithmetic below is written
#' in the covariance, so a precision structure is inverted once here rather
#' than at every use. The log-determinant follows the structure's own, with its
#' sign flipped for a precision, which is the one place the two forms differ in
#' anything but cost.
#'
#' @param distrib A \code{\link{MvGaussianDistrib}} object.
#' @param theta A named list of parameters, already aligned.
#' @param derivs Whether the structure's first derivative matrices are needed.
#' @param derivs2 Whether its second derivatives are needed as well.
#'
#' @return A list with \code{mu}, \code{sigma}, \code{sigma_inv},
#'   \code{logdet}, \code{eta}, and optionally \code{a} and \code{a2}, the
#'   derivatives of the covariance with respect to the free values.
#'
#' @keywords internal
mvg_pieces <- function(distrib, theta, derivs = FALSE, derivs2 = FALSE) {
  p <- distrib@n_dim
  s <- distrib@struct
  v <- mv_flat_theta(distrib, theta)
  mu <- v[seq_len(p)]
  eta <- v[p + seq_len(s@n_free)]

  m <- covstructs7::struct_matrix(s, eta)
  ld <- covstructs7::struct_logdet(s, eta)

  if (distrib@inverted) {
    omega <- m
    sigma <- covstructs7::struct_solve(s, eta)
    sigma_inv <- omega
    logdet <- -ld
  } else {
    sigma <- m
    sigma_inv <- covstructs7::struct_solve(s, eta)
    logdet <- ld
  }

  out <- list(
    mu = unname(mu), sigma = unname(sigma), sigma_inv = unname(sigma_inv),
    logdet = logdet, eta = eta, p = p, s = s
  )

  if (derivs || derivs2) {
    d <- covstructs7::struct_dmatrix(s, eta)
    if (distrib@inverted) {
      # dSigma/deta = -Sigma (dOmega/deta) Sigma, the derivative of an inverse.
      d <- lapply(d, function(ak) -(sigma %*% ak %*% sigma))
    }
    out$a <- lapply(d, unname)
  }
  if (derivs2) {
    d2 <- covstructs7::struct_d2matrix(s, eta)
    if (distrib@inverted) {
      idx <- covstructs7::struct_pair_indices(s)
      om1 <- covstructs7::struct_dmatrix(s, eta)
      d2 <- lapply(seq_along(idx), function(i) {
        k <- idx[[i]][1L]
        l <- idx[[i]][2L]
        # Differentiating -S B_k S once more, with dS/deta_l = -S B_l S.
        sigma %*% (om1[[l]] %*% sigma %*% om1[[k]] +
          om1[[k]] %*% sigma %*% om1[[l]] - d2[[i]]) %*% sigma
      })
      names(d2) <- covstructs7::struct_pair_names(s)
    }
    out$a2 <- lapply(d2, unname)
  }
  out
}


#' The Mean Vector and Covariance a Parameter List Describes
#'
#' @description
#' Assembles the mean vector and the covariance matrix of a multivariate
#' distribution from its flat parameter list.
#'
#' @details
#' The parameters of a multivariate distribution are scalars, so that every
#' generic of the package can index them, and these two functions put them back
#' into the shapes a reader thinks in. \code{mv_sigma()} returns the covariance
#' whichever side the structure parametrises.
#'
#' @param distrib A \code{\link{multivariate_distrib}} object.
#' @param theta A named list or vector of parameters.
#'
#' @return A numeric vector of length \eqn{p} for \code{mv_mu()}, and a
#'   \eqn{p \times p} matrix for \code{mv_sigma()}.
#'
#' @seealso \code{\link{mvgaussian_distrib}}
#'
#' @examples
#' d <- mvgaussian_distrib(2)
#' theta <- list(mu1 = 1, mu2 = -1, log_L1 = 0, log_L2 = 0, L2.1 = 0.5)
#' mv_mu(d, theta)
#' mv_sigma(d, theta)
#'
#' @export
mv_mu <- function(distrib, theta) {
  theta <- align_theta(distrib, theta)
  v <- mv_flat_theta(distrib, theta)
  stats::setNames(
    unname(v[seq_len(distrib@n_dim)]), paste0("v", seq_len(distrib@n_dim))
  )
}

#' @rdname mv_mu
#' @export
mv_sigma <- function(distrib, theta) {
  theta <- align_theta(distrib, theta)
  pc <- mvg_pieces(distrib, theta)
  nm <- paste0("v", seq_len(distrib@n_dim))
  dimnames(pc$sigma) <- list(nm, nm)
  pc$sigma
}


#' Residuals and Whitened Residuals
#'
#' @description
#' The centred response and its image under the inverse covariance, which are
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
#' evaluated row by row. The quadratic form goes through the structure's own
#' factor rather than through an explicit inverse.
#' @param distrib A \code{\link{MvGaussianDistrib}} object.
#' @param y An \eqn{n \times p} matrix of observations, or a vector of length
#'   \eqn{p} for one observation.
#' @param theta A named list of parameters.
#' @param log Logical; if \code{TRUE}, returns the log-density.
#' @return A numeric vector with one value per observation.
#' @keywords internal
S7::method(distrib_pdf, MvGaussianDistrib) <- function(distrib, y, theta,
                                                       log = FALSE) {
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
#' factor taken from the structure where it parametrises the covariance and
#' from a factorisation of the inverse otherwise.
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
#' derivative of \eqn{\Sigma} in the \eqn{k}-th free value of the structure,
#' \deqn{\frac{\partial \ell}{\partial \mu} = w, \qquad
#'   \frac{\partial \ell}{\partial \eta_k} =
#'   -\frac{1}{2}\frac{\partial \log|\Sigma|}{\partial \eta_k}
#'   + \frac{1}{2} w^\top A_k w.}
#' The first term of the second expression is the structure's own
#' \code{struct_dlogdet()}, so no trace is formed here.
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

  dld <- covstructs7::struct_dlogdet(pc$s, pc$eta)
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

  d2ld <- covstructs7::struct_d2logdet(pc$s, pc$eta)
  if (distrib@inverted) d2ld <- -d2ld
  # Sigma^{-1} A_k, formed once: it appears in the mixed block and again in the
  # matrix block.
  sa <- lapply(pc$a, function(ak) si %*% ak)
  spair <- struct_pair_lookup(pc$s)

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
#' corresponding component of \code{struct_d2matrix()}.
#'
#' @details
#' Built from \pkg{covstructs7}'s own enumeration rather than by taking a
#' component key apart, for the reason that package documents: a free value
#' whose label contains the separator splits into the wrong number of pieces.
#'
#' @param s A \pkg{covstructs7} structure.
#'
#' @return A named list of integers, keyed \code{"k:l"} with \eqn{k \le l}.
#'
#' @keywords internal
struct_pair_lookup <- function(s) {
  idx <- covstructs7::struct_pair_indices(s)
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


#' @title Mean of a Multivariate Gaussian
#' @name mean.MvGaussianDistrib
#' @description The mean vector, which is a parameter of the family.
#' @param x A \code{\link{MvGaussianDistrib}} object.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return A numeric vector of length \eqn{p}.
#' @keywords internal
S7::method(mean, MvGaussianDistrib) <- function(x, theta, ...) {
  mv_mu(x, theta)
}


#' @title Variance of a Multivariate Gaussian
#' @name variance.MvGaussianDistrib
#' @description
#' The covariance matrix, which the structure carries. The return is a matrix
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
  s <- distrib@struct
  as.list(stats::setNames(
    c(stats::runif(p, -1, 1), stats::runif(s@n_free, -0.4, 0.4)),
    distrib@params
  ))
}
