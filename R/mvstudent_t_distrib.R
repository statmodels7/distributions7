#' @include multivariate.R mvgaussian_distrib.R mv_plot.R
NULL

#' Multivariate Student's t Distribution
#'
#' @description
#' The S7 class of multivariate Student t distributions: a mean vector, a
#' \pkg{covstructs7} structure for the scale matrix, and a degrees-of-freedom
#' parameter. Constructed by \code{\link{mvstudent_t_distrib}}.
#'
#' @inheritParams multivariate_distrib
#' @param struct The \pkg{covstructs7} structure carrying the scale matrix.
#'
#' @return An object of class \code{MvStudentTDistrib}.
#'
#' @seealso \code{\link{mvstudent_t_distrib}}
#'
#' @examples
#' S7::S7_inherits(mvstudent_t_distrib(2), MvStudentTDistrib)
#'
#' @export
MvStudentTDistrib <- S7::new_class("MvStudentTDistrib",
  parent = multivariate_distrib,
  properties = list(struct = covstructs7::covstruct)
)


#' Construct a Multivariate Student's t Distribution
#'
#' @description
#' The elliptical \eqn{t} on \eqn{\mathbb{R}^p}, with density
#' \deqn{f(y) \propto |\Sigma|^{-1/2}
#'   \left(1 + \frac{(y-\mu)^\top \Sigma^{-1} (y-\mu)}{\nu}\right)^{-(\nu+p)/2}.}
#'
#' @details
#' \eqn{\Sigma} is the \strong{scale} matrix and not the covariance: the
#' covariance is \eqn{\nu\Sigma/(\nu-2)} where it exists at all, and for
#' \eqn{\nu \le 2} it does not while the distribution is perfectly well
#' defined. \code{\link{mv_sigma}} returns the scale matrix, which is the thing
#' the parametrisation carries, and \code{\link{variance}} returns the
#' covariance, which is a moment. Keeping the two apart is what lets a fit run
#' at \eqn{\nu = 1.5}.
#'
#' \strong{Parameters.} The mean contributes \code{mu1}, ..., \code{mup}, the
#' structure contributes its free values under the \code{sigma_} prefix, and
#' \code{nu} is added last. The mean and the structure are unconstrained and
#' carry identity links; \code{nu} is positive and carries a log link by
#' default, so unlike the multivariate gaussian this family's link scale is not
#' its parameter scale.
#'
#' \strong{Reading a fit.} \code{\link{mv_summary}} reports the square roots of
#' the diagonal of the \strong{scale} matrix and the correlations. The
#' correlations are the response's as well, a positive multiple of a matrix
#' leaving them alone, but the diagonal quantities are not standard deviations
#' of the response and are named \code{scale_sd_} to say so.
#'
#' \strong{What it is for.} A gaussian fitted to data with a few outlying rows
#' inflates its covariance to cover them. A \eqn{t} with \eqn{\nu} estimated
#' does not: the observations far from the centre get a weight
#' \eqn{(\nu+p)/(\nu+q)} that falls away with their Mahalanobis distance
#' \eqn{q}, which is what appears in the score below and what makes the fit
#' resistant. The gaussian is the limit \eqn{\nu \to \infty}.
#'
#' \strong{The expected information} has no closed form here and is
#' approximated by sampling. \code{\link{fit_distrib}} therefore accepts
#' \code{approx}, which it refuses for a family that computes it exactly.
#'
#' @param n_dim The dimension \eqn{p}.
#' @param struct_sigma A \pkg{covstructs7} structure for the scale matrix.
#'   Defaults to \code{covstructs7::log_cholesky(n_dim)}.
#' @param link_nu The link for the degrees of freedom. Defaults to
#'   \code{linkfunctions7::log_link()}.
#'
#' @return An object of class \code{\link{MvStudentTDistrib}}.
#'
#' @seealso \code{\link{mvgaussian_distrib}}, \code{\link{student_t_distrib}}
#'
#' @examples
#' d <- mvstudent_t_distrib(2)
#' d@params
#'
#' theta <- list(mu1 = 0, mu2 = 0, sigma_log_L1 = 0, sigma_log_L2 = 0, sigma_L2.1 = 0.4, nu = 5)
#'
#' # the scale matrix is what the structure carries; the covariance is a moment
#' mv_sigma(d, theta)
#' variance(d, theta)
#'
#' # and below two degrees of freedom the covariance does not exist at all,
#' # while the density does
#' theta$nu <- 1.5
#' distrib_pdf(d, c(0, 0), theta)
#' variance(d, theta)
#'
#' @export
mvstudent_t_distrib <- function(n_dim, struct_sigma = NULL,
                                link_nu = linkfunctions7::log_link()) {
  if (!is.numeric(n_dim) || length(n_dim) != 1L || !is.finite(n_dim) ||
    n_dim < 1 || n_dim != round(n_dim)) {
    stop("'n_dim' must be a single positive integer.", call. = FALSE)
  }
  p <- as.integer(n_dim)
  s <- if (is.null(struct_sigma)) covstructs7::log_cholesky(p) else struct_sigma

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
      "The structure's free value '%s' collides with a parameter name.\n",
      "  Every derivative component is keyed by these, so they must be unique."
    ), clash[1L]), call. = FALSE)
  }

  params <- c(mu_names, free_names, "nu")
  n_par <- length(params)
  ident <- linkfunctions7::identity_link()

  MvStudentTDistrib(
    distrib_name = sprintf("multivariate student t [%dd, sigma=%s]", p, s@struct_name),
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
    struct = s
  )
}


#' The Pieces a Multivariate t Evaluates From
#'
#' @description
#' Assembles the location, the scale matrix, its inverse, the log-determinant
#' and the degrees of freedom from a flat parameter vector, with the
#' structure's derivative matrices when they are needed.
#'
#' @param distrib A \code{\link{MvStudentTDistrib}} object.
#' @param theta A named list of parameters, already aligned.
#' @param derivs Whether the first derivative matrices are needed.
#' @param derivs2 Whether the second derivatives are needed as well.
#'
#' @return A list with \code{mu}, \code{sigma}, \code{sigma_inv},
#'   \code{logdet}, \code{nu}, \code{eta}, \code{p}, \code{s}, and optionally
#'   \code{a} and \code{a2}.
#'
#' @keywords internal
mvt_pieces <- function(distrib, theta, derivs = FALSE, derivs2 = FALSE) {
  p <- distrib@n_dim
  s <- distrib@struct
  v <- mv_flat_theta(distrib, theta)
  mu <- unname(v[seq_len(p)])
  eta <- unname(v[p + seq_len(s@n_free)])
  nu <- unname(v[[length(v)]])

  out <- list(
    mu = mu, eta = eta, nu = nu, p = p, s = s,
    sigma = unname(covstructs7::struct_matrix(s, eta)),
    sigma_inv = unname(covstructs7::struct_solve(s, eta)),
    logdet = covstructs7::struct_logdet(s, eta)
  )
  if (derivs || derivs2) {
    out$a <- lapply(covstructs7::struct_dmatrix(s, eta), unname)
  }
  if (derivs2) {
    out$a2 <- lapply(covstructs7::struct_d2matrix(s, eta), unname)
  }
  out
}


#' @title Multivariate Student t Density
#' @name distrib_pdf.MvStudentTDistrib
#' @description
#' \deqn{\ell = \log\Gamma\!\left(\tfrac{\nu+p}{2}\right)
#'   - \log\Gamma\!\left(\tfrac{\nu}{2}\right)
#'   - \tfrac{p}{2}\log(\nu\pi) - \tfrac{1}{2}\log|\Sigma|
#'   - \tfrac{\nu+p}{2}\log\!\left(1 + \tfrac{q}{\nu}\right),}
#' with \eqn{q = (y-\mu)^\top \Sigma^{-1}(y-\mu)}. The logarithm is taken with
#' \code{log1p}, which is the difference between a number and a loss of every
#' significant digit when \eqn{q/\nu} is small.
#' @param distrib A \code{\link{MvStudentTDistrib}} object.
#' @param y An \eqn{n \times p} matrix of observations.
#' @param theta A named list of parameters.
#' @param log Logical; if \code{TRUE}, returns the log-density.
#' @return A numeric vector with one value per observation.
#' @keywords internal
S7::method(distrib_pdf, MvStudentTDistrib) <- function(distrib, y, theta,
                                                       log = FALSE) {
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
#' @description
#' The scale-mixture representation: \eqn{\mu + L z / \sqrt{g/\nu}} with
#' \eqn{z} standard normal, \eqn{g \sim \chi^2_\nu} and \eqn{LL^\top = \Sigma}.
#' A t is a gaussian whose precision has been multiplied by a gamma variate,
#' which is the same fact that makes it robust.
#' @param distrib A \code{\link{MvStudentTDistrib}} object.
#' @param n The number of observations to draw.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return An \eqn{n \times p} numeric matrix.
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


#' The Weight a Multivariate t Gives Each Observation
#'
#' @description
#' \eqn{(\nu + p)/(\nu + q)} at each observation, with the whitened residuals
#' that go with it.
#'
#' @details
#' This weight is the whole of the family's robustness. At \eqn{q = 0} it is
#' \eqn{(\nu+p)/\nu} and it decays like \eqn{1/q}, so an observation far from
#' the centre contributes less to every derivative rather than dragging the
#' fit towards itself; letting \eqn{\nu \to \infty} sends it to one and
#' recovers the gaussian, where nothing is downweighted.
#'
#' @param y An \eqn{n \times p} matrix.
#' @param pc The result of \code{\link{mvt_pieces}}.
#'
#' @return A list with \code{r}, \code{w}, \code{q} and \code{cw}.
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
#' @description
#' Closed form. With \eqn{w = \Sigma^{-1}(y-\mu)}, \eqn{q = (y-\mu)^\top w} and
#' \eqn{c = (\nu+p)/(\nu+q)},
#' \deqn{\partial_\mu \ell = c\,w, \qquad
#'   \partial_{\eta_k}\ell = -\tfrac{1}{2}\partial_{\eta_k}\log|\Sigma|
#'   + \tfrac{c}{2}\, w^\top A_k w,}
#' \deqn{\partial_\nu \ell = \tfrac{1}{2}\left[
#'   \psi\!\left(\tfrac{\nu+p}{2}\right) - \psi\!\left(\tfrac{\nu}{2}\right)
#'   - \tfrac{p}{\nu} - \log\!\left(1+\tfrac{q}{\nu}\right)
#'   + \tfrac{(\nu+p)q}{\nu(\nu+q)}\right].}
#' The gaussian score is the limit \eqn{c \to 1}.
#' @param distrib A \code{\link{MvStudentTDistrib}} object.
#' @param y An \eqn{n \times p} matrix of observations.
#' @param theta A named list of parameters.
#' @param scale Handled by the generic before dispatch.
#' @param ... Unused.
#' @return A named list with one numeric vector per parameter.
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

  dld <- covstructs7::struct_dlogdet(pc$s, pc$eta)
  for (k in seq_along(pc$a)) {
    out[[p + k]] <- -0.5 * dld[[k]] +
      0.5 * z$cw * rowSums((z$w %*% pc$a[[k]]) * z$w)
  }

  out[[distrib@n_params]] <- 0.5 * (
    digamma((nu + p) / 2) - digamma(nu / 2) - p / nu -
      base::log1p(z$q / nu) + (nu + p) * z$q / (nu * (nu + z$q))
  )
  out
}


#' @title Multivariate Student t Observed Hessian
#' @name distrib_hessian.MvStudentTDistrib
#' @description
#' Closed form, obtained by differentiating the score of
#' \code{\link[=distrib_gradient.MvStudentTDistrib]{distrib_gradient()}} once
#' more. Every block picks up a term in \eqn{\partial c/\partial\cdot}, because
#' the weight depends on the observation through \eqn{q}; that dependence is
#' what distinguishes the family from the gaussian, where \eqn{c} is one and
#' those terms are absent.
#' @param distrib A \code{\link{MvStudentTDistrib}} object.
#' @param y An \eqn{n \times p} matrix of observations.
#' @param theta A named list of parameters.
#' @param scale Handled by the generic before dispatch.
#' @param ... Unused.
#' @return A named list keyed as \code{\link{hess_names}(distrib@params)}.
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

  d2ld <- covstructs7::struct_d2logdet(pc$s, pc$eta)
  spair <- struct_pair_lookup(pc$s)
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
      0.5 * (
        0.5 * trigamma((nu + p) / 2) - 0.5 * trigamma(nu / 2) + p / nu^2 +
          z$q / (nu * den) -
          z$q * (nu^2 + 2 * p * nu + p * z$q) / (nu * den)^2
      )
    }
    if (length(out[[i]]) == 1L) out[[i]] <- rep(out[[i]], n)
  }
  out
}


#' @title Multivariate Student t Response Gradient
#' @name distrib_grad_y.MvStudentTDistrib
#' @description
#' \eqn{\partial \ell / \partial y = -c\,\Sigma^{-1}(y-\mu)}, the gaussian
#' expression with the family's weight in front of it.
#' @param distrib A \code{\link{MvStudentTDistrib}} object.
#' @param y An \eqn{n \times p} matrix of observations.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return An \eqn{n \times p} numeric matrix.
#' @keywords internal
S7::method(distrib_grad_y, MvStudentTDistrib) <- function(distrib, y, theta, ...) {
  y <- as_mv_matrix(distrib, y)
  pc <- mvt_pieces(distrib, theta)
  z <- mvt_weights(y, pc)
  -z$cw * z$w
}


#' @title Multivariate Student t Response Hessian
#' @name distrib_hess_y.MvStudentTDistrib
#' @description
#' Closed form. With \eqn{w = \Sigma^{-1}(y-\mu)}, \eqn{q = (y-\mu)^\top w} and
#' \eqn{c = (\nu+p)/(\nu+q)},
#' \deqn{\dfrac{\partial^2 \ell}{\partial y \, \partial y^\top}
#'   = -c\,\Sigma^{-1} + \dfrac{2c}{\nu+q}\, w w^\top,}
#' which depends on the observation through \eqn{c} and \eqn{w} --- unlike the
#' gaussian's, which is \eqn{-\Sigma^{-1}} everywhere --- so one matrix is
#' returned per row.
#' @param distrib A \code{\link{MvStudentTDistrib}} object.
#' @param y An \eqn{n \times p} matrix of observations.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return A \eqn{p \times p \times n} numeric array.
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
#' @description
#' The location vector, which is the mean when \eqn{\nu > 1} and undefined
#' otherwise; \code{NaN} is returned there rather than the location, since the
#' location exists as a parameter while the moment does not.
#' @param x A \code{\link{MvStudentTDistrib}} object.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return A numeric vector of length \eqn{p}.
#' @keywords internal
S7::method(mean, MvStudentTDistrib) <- function(x, theta, ...) {
  mu <- mv_location(x, theta)
  nu <- mv_flat_theta(x, align_theta(x, theta))[[x@n_params]]
  if (nu <= 1) mu[] <- NaN
  mu
}


#' @title Covariance of a Multivariate Student t
#' @name variance.MvStudentTDistrib
#' @description
#' \eqn{\nu\Sigma/(\nu-2)} for \eqn{\nu > 2}, and infinite otherwise. This is
#' the moment; the matrix the parametrisation carries is
#' \code{\link{mv_sigma}}, and the two differ by the factor above.
#' @param x A \code{\link{MvStudentTDistrib}} object.
#' @param theta A named list of parameters.
#' @param ... Unused.
#' @return A \eqn{p \times p} numeric matrix.
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
#' @description
#' A marginal of a \eqn{t} is a \eqn{t} with the same degrees of freedom, the
#' subvector of the location and the corresponding block of the scale matrix.
#' The degrees of freedom do not change with the dimension, which is what makes
#' the family closed under marginalisation at all.
#' @param distrib A \code{\link{MvStudentTDistrib}} object.
#' @param theta A named list of parameters.
#' @param which An integer vector of coordinates.
#' @param ... Unused.
#' @return A list with \code{distrib} and \code{theta}.
#' @keywords internal
S7::method(mv_marginal, MvStudentTDistrib) <- function(distrib, theta, which, ...) {
  mu <- as.numeric(mv_location(distrib, theta))[which]
  sg <- mv_sigma(distrib, theta)[which, which, drop = FALSE]
  nu <- mv_flat_theta(distrib, align_theta(distrib, theta))[[distrib@n_params]]
  md <- mvstudent_t_distrib(length(which))
  eta <- covstructs7::struct_free(md@struct, unname(sg))
  list(
    distrib = md,
    theta = as.list(stats::setNames(c(mu, unname(eta), nu), md@params))
  )
}


#' @title The Scale Matrix of a Multivariate Student t
#' @name mv_sigma.MvStudentTDistrib
#' @description
#' The scale matrix, not the covariance. The two differ by
#' \eqn{\nu/(\nu-2)}, and only the scale matrix exists for every admissible
#' \eqn{\nu}; the covariance is \code{\link{variance}}.
#' @param distrib A \code{\link{MvStudentTDistrib}} object.
#' @param theta A named list of parameters.
#' @return A \eqn{p \times p} numeric matrix.
#' @keywords internal
S7::method(mv_sigma, MvStudentTDistrib) <- function(distrib, theta) {
  pc <- mvt_pieces(distrib, theta)
  nm <- paste0("v", seq_len(distrib@n_dim))
  dimnames(pc$sigma) <- list(nm, nm)
  pc$sigma
}


#' @title Random Parameters for a Multivariate Student t
#' @name generate_random_theta.MvStudentTDistrib
#' @description
#' A location near the origin, a scale near the identity, and degrees of
#' freedom in a range where the family is heavy-tailed but the likelihood is
#' still well conditioned.
#' @param distrib A \code{\link{MvStudentTDistrib}} object.
#' @param ... Unused.
#' @return A named list of scalars.
#' @keywords internal
S7::method(generate_random_theta, MvStudentTDistrib) <- function(distrib, ...) {
  p <- distrib@n_dim
  s <- distrib@struct
  as.list(stats::setNames(
    c(stats::runif(p, -1, 1), stats::runif(s@n_free, -0.4, 0.4),
      stats::runif(1, 3, 12)),
    distrib@params
  ))
}

#' @title Location of a Multivariate Student t
#' @name mv_location.MvStudentTDistrib
#' @description
#' The first \eqn{p} parameters. They are the mean when \eqn{\nu > 1} and the
#' centre of symmetry always, which is why the generic is called a location
#' rather than a mean.
#' @param distrib A \code{\link{MvStudentTDistrib}} object.
#' @param theta A named list of parameters.
#' @return A named numeric vector of length \eqn{p}.
#' @keywords internal
S7::method(mv_location, MvStudentTDistrib) <- mv_leading_location
