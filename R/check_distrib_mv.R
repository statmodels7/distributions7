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
#' them go with them. Running them anyway and reporting the rejections as failures
#' is the mistake a validator makes when it does not know about a case, and it
#' is worse than not checking, because a user validating their own distribution
#' cannot tell a real defect from it.
#'
#' What replaces them are checks that do generalize:
#' \enumerate{
#'   \item \strong{the density is positive and finite} on a sample;
#'   \item \strong{the density integrates to one}. For a family that enumerates
#'     its support through \code{\link{mv_support}} this is an exact sum over
#'     that support; otherwise it is importance sampling from the proposal
#'     \code{\link{mv_reference_draw}} supplies, which by default is a gaussian
#'     with the same mean and an inflated covariance. The proposal is
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
#' The last of these is emitted only when it applies, as the univariate battery
#' already omits the checks that a discrete family has no counterpart for. A
#' family with an enumerable support is discrete and has no derivative in the
#' response, and the multivariate base class rejects
#' \code{\link{distrib_grad_y}} by design, so a family that has not registered
#' one has made a choice rather than left a gap.
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

  # Two facts decide which checks apply. A family that enumerates a support is
  # discrete, so its normalization is an exact sum and it has no derivative in
  # the response; and the base class refuses the response derivatives by
  # design, so a family that has not registered them has made a choice rather
  # than left a gap, and reporting that refusal as a failure would be the
  # mistake this battery exists to avoid.
  enumerable <- has_mv_support(distrib)
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
    if (enumerable) {
      # A finite support makes the normalization an exact sum, so the check is
      # an equality rather than a comparison against Monte Carlo error.
      total <- sum(distrib_pdf(distrib, mv_support(distrib, theta), theta))
      new_check("density integrates to 1", abs(total - 1) < 1e-10,
        abs(total - 1))
    } else {
      prop <- mv_reference_draw(distrib, theta, nsim)
      lf <- distrib_pdf(distrib, prop$y, theta, log = TRUE)
      total <- mean(exp(lf - prop$logd))
      new_check("density integrates to 1", abs(total - 1) < 20 / sqrt(nsim),
        abs(total - 1))
    }
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
  if (!enumerable && has_mv_grad_y(distrib)) {
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
  }

  res
}


#' Whether a Multivariate Family Enumerates Its Support
#'
#' @description
#' \code{TRUE} when the family registers \code{\link{mv_support}}, which is
#' what a discrete multivariate family does and a continuous one cannot.
#'
#' @details
#' The question is asked of the method rather than of the class, because the
#' multivariate branch sits beside \code{\link{continuous_distrib}} and
#' \code{\link{discrete_distrib}} rather than under either, so there is no
#' class to test. The owning class of a method is read through
#' \code{\link{is_class}}, never with \code{identical()}, which is object
#' identity and fails whenever the package's code is re-evaluated rather than
#' loaded.
#'
#' @param x An object inheriting from class \code{\link{multivariate_distrib}}.
#'
#' @return A single logical.
#'
#' @seealso \code{\link{mv_support}}, \code{\link{check_distrib}}
#' @keywords internal
has_mv_support <- function(x) {
  m <- tryCatch(S7::method(mv_support, S7::S7_class(x)), error = function(e) NULL)
  !is.null(m) && !is_class(attr(m, "signature")[[1L]], multivariate_distrib)
}

#' Whether a Multivariate Family Implements Its Response Derivatives
#'
#' @description
#' \code{TRUE} when \code{\link{distrib_grad_y}} comes from the family rather
#' than from the base class, whose method rejects.
#'
#' @param x An object inheriting from class \code{\link{multivariate_distrib}}.
#'
#' @return A single logical.
#'
#' @seealso \code{\link{has_mv_support}}, \code{\link{check_distrib}}
#' @keywords internal
has_mv_grad_y <- function(x) {
  m <- tryCatch(S7::method(distrib_grad_y, S7::S7_class(x)), error = function(e) NULL)
  !is.null(m) && !is_class(attr(m, "signature")[[1L]], multivariate_distrib)
}


#' A Central-Difference Gradient Without a Dependency
#'
#' @description
#' One central difference per coordinate, which is all the multivariate checks
#' need and which keeps \pkg{numDeriv} a suggestion rather than a requirement.
#'
#' @details
#' The nodes, the weights and the step are \pkg{numericals7}'s.
#' \code{\link[numericals7]{fd_derivative}} is not called directly because its
#' \code{f} maps a vector of points to the values at those points, while this
#' one maps a whole vector to a single number.
#'
#' @param f A function of a numeric vector, returning one number.
#' @param x The point.
#' @param h_rel Deprecated and unused; the step is
#'   \code{\link[numericals7]{fd_step}}'s.
#'
#' @return A numeric vector the length of \code{x}.
#'
#' @seealso \code{\link[numericals7]{fd_weights}}
#'
#' @keywords internal
numDeriv_grad <- function(f, x, h_rel = NULL) {
  s <- numericals7::fd_offsets(1L, accuracy = 2L)$central
  w <- numericals7::fd_weights(s, 1L)
  vapply(seq_along(x), function(k) {
    h <- numericals7::fd_step(x[k], 1L, accuracy = 2L)
    acc <- 0
    for (j in seq_along(s)) {
      if (w[j] == 0) next
      z <- x
      z[k] <- x[k] + s[j] * h
      acc <- acc + w[j] * f(z)
    }
    acc / h
  }, numeric(1))
}


#' A Second Derivative From One Stencil
#'
#' @description
#' The mixed or repeated second derivative of a scalar function, from a single
#' stencil rather than from nested first differences.
#'
#' @details
#' The nodes, the weights and the step are \pkg{numericals7}'s. Where the two
#' coordinates differ the stencil is the product of two first-order factors,
#' which is one stencil in two variables and not a difference of a difference:
#' nesting is forbidden along ONE variable, and is what a mixed derivative is
#' along two.
#'
#' @param f A function of a numeric vector, returning one number.
#' @param x The point.
#' @param k,l The coordinates to differentiate in.
#' @param h_rel Deprecated and unused; the step is
#'   \code{\link[numericals7]{fd_step}}'s at order two.
#'
#' @return A single number.
#'
#' @seealso \code{\link[numericals7]{fd_weights}}
#'
#' @keywords internal
fd_second <- function(f, x, k, l, h_rel = NULL) {
  hk <- numericals7::fd_step(x[k], 2L, accuracy = 2L)
  hl <- numericals7::fd_step(x[l], 2L, accuracy = 2L)
  if (k == l) {
    s <- numericals7::fd_offsets(2L, accuracy = 2L)$central
    w <- numericals7::fd_weights(s, 2L)
    acc <- 0
    for (j in seq_along(s)) {
      if (w[j] == 0) next
      z <- x
      z[k] <- x[k] + s[j] * hk
      acc <- acc + w[j] * f(z)
    }
    return(acc / hk^2)
  }
  s <- numericals7::fd_offsets(1L, accuracy = 2L)$central
  w <- numericals7::fd_weights(s, 1L)
  acc <- 0
  for (a in seq_along(s)) {
    if (w[a] == 0) next
    for (b in seq_along(s)) {
      if (w[b] == 0) next
      z <- x
      z[k] <- x[k] + s[a] * hk
      z[l] <- x[l] + s[b] * hl
      acc <- acc + w[a] * w[b] * f(z)
    }
  }
  acc / (hk * hl)
}
