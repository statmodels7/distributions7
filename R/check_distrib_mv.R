#' @include check_distrib.R multivariate.R
NULL

#' Validate a Multivariate Distribution
#'
#' @description
#' The battery \code{\link{check_distrib}} runs on a multivariate distribution,
#' where the one-dimensional checks have no counterpart.
#'
#' @details
#' Five of the univariate checks do not survive the move to \eqn{p} dimensions:
#' the distribution function is an integral over an orthant, the quantile
#' function inverts an ordering that does not exist, and the two checks built on
#' them go with them. Running them anyway and reporting the refusals as failures
#' is the mistake a validator makes when it does not know about a case, and it
#' is worse than not checking, because a user validating their own distribution
#' cannot tell a real defect from it.
#'
#' What replaces them are checks that do generalise:
#' \enumerate{
#'   \item \strong{the density is positive and finite} on a sample;
#'   \item \strong{the density integrates to one}, by importance sampling from a
#'     gaussian with the same mean and an inflated covariance. The proposal is
#'     deliberately not the distribution itself, which would make the ratio
#'     identically one and the check vacuous;
#'   \item \strong{the score has mean zero}, the first Bartlett identity, under
#'     the distribution's own generator -- so it is also a check that the
#'     generator and the density describe the same law;
#'   \item \strong{gradient and Hessian} against finite differences of the
#'     summed log-density;
#'   \item \strong{the expected information} against the Monte Carlo average of
#'     the observed Hessian, and against the variance of the score, which is the
#'     second Bartlett identity;
#'   \item \strong{the generator} against the first two moments;
#'   \item \strong{the response derivatives} against finite differences in
#'     \eqn{y}.
#' }
#'
#' @param distrib A \code{\link{multivariate_distrib}} object.
#' @param theta A named list of parameters, already aligned.
#' @param n The number of observations drawn for the derivative checks.
#' @param nsim The Monte Carlo sample size.
#' @param tol The relative tolerance.
#'
#' @return A list of one-row data frames, as \code{\link{check_distrib}} builds.
#'
#' @seealso \code{\link{check_distrib}}
#' @keywords internal
check_distrib_mv <- function(distrib, theta, n, nsim, tol) {
  p <- distrib@n_dim
  params <- distrib@params
  res <- list()
  rel <- function(a, e) max(abs(a - e) / pmax(1, abs(e)))

  y <- distrib_rng(distrib, n, theta)
  v0 <- vapply(theta[seq_len(distrib@n_params)], function(z) z[[1L]], numeric(1))
  ll <- function(v) {
    sum(distrib_pdf(distrib, y, as.list(stats::setNames(v, params)), log = TRUE))
  }

  # --- the density ---------------------------------------------------------
  res[[length(res) + 1L]] <- safe_check("density is positive and finite", {
    f <- distrib_pdf(distrib, y, theta)
    new_check("density is positive and finite", all(is.finite(f) & f > 0), NA_real_)
  })

  res[[length(res) + 1L]] <- safe_check("density integrates to 1", {
    # Importance sampling from a gaussian with the same mean and an inflated
    # covariance: a proposal equal to the distribution itself would make every
    # ratio one and certify nothing.
    mu <- as.numeric(mean(distrib, theta))
    sg <- as.matrix(variance(distrib, theta)) * 2
    l <- t(chol(sg))
    z <- matrix(stats::rnorm(nsim * p), nsim, p)
    x <- sweep(z %*% t(l), 2L, mu, "+")
    # The whitened residuals are L^{-1} r, so the system is LOWER triangular
    # and forwardsolve is what solves it; backsolve on the transpose solves
    # L' x = b, which is a different vector with the same shape.
    q <- colSums(forwardsolve(l, t(sweep(x, 2L, mu, "-")))^2)
    lg <- -0.5 * (p * log(2 * pi) + 2 * sum(log(diag(l))) + q)
    lf <- distrib_pdf(distrib, x, theta, log = TRUE)
    total <- mean(exp(lf - lg))
    new_check("density integrates to 1", abs(total - 1) < 20 / sqrt(nsim),
      abs(total - 1))
  })

  # --- the derivatives -----------------------------------------------------
  res[[length(res) + 1L]] <- safe_check("gradient vs finite differences", {
    a <- vapply(distrib_gradient(distrib, y, theta), sum, numeric(1))
    e <- numDeriv_grad(ll, v0)
    new_check("gradient vs finite differences", rel(a, e) < tol, rel(a, e))
  })

  res[[length(res) + 1L]] <- safe_check("hessian vs finite differences", {
    h <- distrib_hessian(distrib, y, theta)
    a <- vapply(hess_names(params), function(nm) sum(h[[nm]]), numeric(1))
    pr <- hess_pairs(params)
    e <- vapply(hess_names(params), function(nm) {
      k <- match(pr[[nm]], params)
      fd_second(ll, v0, k[1L], k[2L])
    }, numeric(1))
    new_check("hessian vs finite differences", rel(a, e) < tol * 10, rel(a, e))
  })

  # --- the information -----------------------------------------------------
  res[[length(res) + 1L]] <- safe_check("expected information vs Monte Carlo", {
    big <- distrib_rng(distrib, nsim, theta)
    hb <- distrib_hessian(distrib, big, theta)
    eh <- distrib_expected_hessian(distrib, big[1L, , drop = FALSE], theta)
    nm <- hess_names(params)
    a <- vapply(nm, function(k) eh[[k]][1L], numeric(1))
    e <- vapply(nm, function(k) mean(hb[[k]]), numeric(1))
    new_check("expected information vs Monte Carlo", rel(a, e) < 0.05, rel(a, e))
  })

  res[[length(res) + 1L]] <- safe_check("score has mean zero", {
    big <- distrib_rng(distrib, nsim, theta)
    s <- do.call(cbind, distrib_gradient(distrib, big, theta))
    m <- colMeans(s)
    # Against its own Monte Carlo standard error, which is the only scale the
    # comparison has: a score component with a large variance has a noisy mean.
    se <- apply(s, 2L, stats::sd) / sqrt(nsim)
    z <- max(abs(m) / pmax(se, 1e-12))
    new_check("score has mean zero", z < 5, z)
  })

  res[[length(res) + 1L]] <- safe_check("information matches the score variance", {
    big <- distrib_rng(distrib, nsim, theta)
    s <- do.call(cbind, distrib_gradient(distrib, big, theta))
    opg <- crossprod(sweep(s, 2L, colMeans(s))) / nrow(s)
    eh <- distrib_expected_hessian(distrib, big[1L, , drop = FALSE], theta)
    pr <- hess_pairs(params)
    e_mat <- matrix(0, length(params), length(params))
    for (nm in hess_names(params)) {
      k <- match(pr[[nm]], params)
      e_mat[k[1L], k[2L]] <- eh[[nm]][1L]
      e_mat[k[2L], k[1L]] <- eh[[nm]][1L]
    }
    d <- max(abs(opg + e_mat)) / max(abs(e_mat))
    new_check("information matches the score variance", d < 0.05, d)
  })

  # --- the generator -------------------------------------------------------
  res[[length(res) + 1L]] <- safe_check("rng matches the first two moments", {
    big <- distrib_rng(distrib, nsim, theta)
    mu <- as.numeric(mean(distrib, theta))
    sg <- as.matrix(variance(distrib, theta))
    sc <- sqrt(max(diag(sg)))
    d <- max(
      max(abs(colMeans(big) - mu)) / sc,
      max(abs(stats::cov(big) - sg)) / max(abs(sg))
    )
    new_check("rng matches the first two moments", d < 20 / sqrt(nsim), d)
  })

  # --- the response --------------------------------------------------------
  res[[length(res) + 1L]] <- safe_check("response derivatives vs finite differences", {
    a <- distrib_grad_y(distrib, y, theta)
    e <- t(vapply(seq_len(nrow(y)), function(i) {
      numDeriv_grad(
        function(z) distrib_pdf(distrib, matrix(z, 1L), theta, log = TRUE),
        y[i, ]
      )
    }, numeric(p)))
    new_check("response derivatives vs finite differences",
      rel(a, e) < tol, rel(a, e))
  })

  res
}


#' A Central-Difference Gradient Without a Dependency
#'
#' @description
#' One central difference per coordinate, which is all the multivariate checks
#' need and which keeps \pkg{numDeriv} a suggestion rather than a requirement.
#'
#' @param f A function of a numeric vector, returning one number.
#' @param x The point.
#' @param h_rel The relative step.
#'
#' @return A numeric vector the length of \code{x}.
#'
#' @keywords internal
numDeriv_grad <- function(f, x, h_rel = .Machine$double.eps^(1 / 3)) {
  vapply(seq_along(x), function(k) {
    h <- h_rel * max(1, abs(x[k]))
    up <- x
    dn <- x
    up[k] <- x[k] + h
    dn[k] <- x[k] - h
    (f(up) - f(dn)) / (2 * h)
  }, numeric(1))
}


#' A Second Derivative From One Stencil
#'
#' @description
#' The mixed or repeated second derivative of a scalar function, from a single
#' stencil rather than from nested first differences.
#'
#' @param f A function of a numeric vector, returning one number.
#' @param x The point.
#' @param k,l The coordinates to differentiate in.
#' @param h_rel The relative step.
#'
#' @return A single number.
#'
#' @keywords internal
fd_second <- function(f, x, k, l, h_rel = .Machine$double.eps^(1 / 4)) {
  hk <- h_rel * max(1, abs(x[k]))
  hl <- h_rel * max(1, abs(x[l]))
  if (k == l) {
    up <- dn <- x
    up[k] <- x[k] + hk
    dn[k] <- x[k] - hk
    return((f(up) - 2 * f(x) + f(dn)) / hk^2)
  }
  pp <- pm <- mp <- mm <- x
  pp[k] <- pp[k] + hk; pp[l] <- pp[l] + hl
  pm[k] <- pm[k] + hk; pm[l] <- pm[l] - hl
  mp[k] <- mp[k] - hk; mp[l] <- mp[l] + hl
  mm[k] <- mm[k] - hk; mm[l] <- mm[l] - hl
  (f(pp) - f(pm) - f(mp) + f(mm)) / (4 * hk * hl)
}
