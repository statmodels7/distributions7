#' @include distrib.R generics.R utility_functions.R numerical_derivatives.R higher_derivatives.R y_derivatives.R
NULL

#' Record One Check Result
#'
#' @description
#' Builds the single-row data frame that \code{\link{check_distrib}} accumulates
#' into its report.
#'
#' @param name The check's name, as it appears in the report.
#' @param ok Logical; whether the check passed.
#' @param stat The numeric statistic the check produced, typically a maximum
#'   discrepancy.
#' @param detail An optional message, used to carry the reason for a failure.
#'
#' @return A one-row data frame with columns \code{check}, \code{status},
#'   \code{statistic} and \code{detail}.
#'
#' @seealso \code{\link{check_distrib}}, \code{\link{safe_check}}
#' @keywords internal
new_check <- function(name, ok, stat, detail = NA_character_) {
  data.frame(
    check = name, status = if (isTRUE(ok)) "OK" else "FAIL",
    statistic = stat, detail = detail, stringsAsFactors = FALSE
  )
}

#' Run a Check, Turning an Error Into a Failure
#'
#' @description
#' Evaluates a check expression and converts any error into a failed row rather
#' than letting it abort the report.
#'
#' @details
#' A distribution under validation is by assumption possibly broken, so a check
#' that throws is itself a result. Without this, the first component to raise
#' would end the run and hide every check after it -- the least useful moment to
#' stop being informative.
#'
#' @param name The check's name, used for the row built on failure.
#' @param expr The expression to evaluate; normally returns a row from
#'   \code{\link{new_check}}.
#'
#' @return The value of \code{expr}, or a failed row carrying the error message.
#'
#' @seealso \code{\link{check_distrib}}, \code{\link{new_check}}
#' @keywords internal
safe_check <- function(name, expr) {
  tryCatch(expr, error = function(e) new_check(name, FALSE, NA_real_, conditionMessage(e)))
}

#' Numerically Validate a Distribution
#'
#' @description
#' Runs a battery of numerical self-consistency checks on a \code{distrib} object.
#' This is the recommended way to validate a distribution you have defined
#' yourself: it verifies that the density integrates (or sums) to one, that the
#' CDF, quantile function and random generator agree with each other and with the
#' density, and that every analytical derivative matches its finite-difference
#' counterpart.
#'
#' @param distrib An object inheriting from class \code{"distrib"}.
#' @param theta A named list of parameter values at which to run the checks. If
#'   \code{NULL} (default) a random admissible value is drawn with
#'   \code{\link{generate_random_theta}}.
#' @param n Integer. Number of observations used for the derivative comparisons.
#'   Defaults to 100.
#' @param nsim Integer. Monte Carlo sample size used for the random generator and
#'   expected-information checks. Defaults to 200000.
#' @param orders Integer vector. Which parameter-derivative orders to check.
#'   Defaults to \code{1:4}; use e.g. \code{1:2} for a faster run.
#' @param tol Numeric. Relative tolerance for the finite-difference comparisons.
#'   Defaults to \code{1e-3}.
#' @param verbose Logical. If \code{TRUE} (default) a readable report is printed.
#'
#' @return Invisibly, a \code{data.frame} with one row per check and columns
#'   \code{check}, \code{status} (\code{"OK"} or \code{"FAIL"}), \code{statistic}
#'   and \code{detail}.
#'
#' @details
#' The checks performed are:
#' \itemize{
#'   \item \strong{density}: non-negativity and integration/summation to 1 over the support.
#'   \item \strong{cdf}: values in \eqn{[0,1]} and monotonicity along a grid of quantiles.
#'   \item \strong{quantile}: round-trip against the CDF (\eqn{F(Q(p)) = p} for continuous
#'     distributions, and the generalized-inverse inequalities for discrete ones).
#'   \item \strong{rng}: the sample mean and variance of a large draw agree with
#'     \code{\link{mean}} and \code{\link{variance}} within Monte Carlo error.
#'   \item \strong{gradient, hessian, deriv3, deriv4}: analytical values against
#'     \code{\link{numerical_gradient}}, \code{\link{numerical_hessian}},
#'     \code{\link{numerical_deriv3}} and \code{\link{numerical_deriv4}}.
#'   \item \strong{expected information}: \code{\link{distrib_expected_hessian}} against a
#'     Monte Carlo estimate of \eqn{-\mathbb{E}[\nabla\ell\,\nabla\ell^\top]}. The outer
#'     product of the score is used as reference because it remains valid when the
#'     log-likelihood is not differentiable in a parameter (see \code{\link{laplace_distrib}}).
#'   \item \strong{response derivatives} (continuous only): \code{\link{distrib_grad_y}} and
#'     \code{\link{distrib_hess_y}} against finite differences in \eqn{y}.
#'   \item \strong{link scale}: \code{scale = "link"} derivatives against finite
#'     differences of the log-likelihood in \eqn{\eta}.
#' }
#' Distributions that rely on the numerical fallbacks will trivially pass the
#' corresponding derivative checks, since analytical and numerical values then
#' coincide by construction.
#'
#' Mixed distributions --- a density with point masses on top of it, as produced by
#' \code{\link{zero_adjusted}()} on a continuous parent --- are handled as long as they
#' declare their atoms through \code{\link{distrib_atoms}}. The density is then expected
#' to integrate to one minus the atomic mass, quantiles falling inside a jump of the CDF
#' are checked as generalized inverses rather than exact ones, and finite differences in
#' \eqn{y} are kept away from the atoms, where no derivative exists.
#'
#' @examples
#' \dontrun{
#' check_distrib(gaussian_distrib())
#' check_distrib(laplace_distrib(), theta = list(mu = 1, b = 2))
#' check_distrib(poisson_distrib(), orders = 1:2, nsim = 5e4)
#' }
#'
#' @seealso \code{\link{numerical_gradient}}, \code{\link{numerical_hessian}},
#'   \code{\link{link_scale_derivatives}}
#' @export
check_distrib <- function(distrib, theta = NULL, n = 100, nsim = 2e5,
                          orders = 1:4, tol = 1e-3, verbose = TRUE) {
  if (is.null(theta)) theta <- generate_random_theta(distrib)
  theta <- align_theta(distrib, theta)

  is_cont <- S7::S7_inherits(distrib, continuous_distrib)
  b <- distrib@bounds
  res <- list()
  smooth_all <- all(param_smoothness(distrib))

  rel <- function(a, e) max(abs(a - e) / pmax(1, abs(e)))

  # A mixed distribution -- zero_adjusted() of a continuous parent is the one the
  # package ships -- is not covered by either branch below. Its density
  # integrates to 1 minus the atomic mass, its cdf jumps, and a finite difference
  # straddling an atom returns a number for a derivative that does not exist. Run
  # as if it were purely continuous it fails four checks while being perfectly
  # correct, which is worse than useless: a user validating their own mixed
  # distribution would have no way to tell a real defect from this one. The atoms
  # are asked for instead, and each check is told what to do about them.
  atoms <- distrib_atoms(distrib, theta)
  atom_mass <- if (length(atoms$y)) sum(atoms$p) else 0
  # Points far enough from every atom that differencing does not cross one.
  away_from_atoms <- function(x, rel_margin = 1e-3) {
    if (!length(atoms$y)) return(rep(TRUE, length(x)))
    keep <- rep(TRUE, length(x))
    for (a in atoms$y) keep <- keep & abs(x - a) > rel_margin * max(1, abs(a))
    keep
  }

  # --- density -------------------------------------------------------------
  res[[length(res) + 1L]] <- safe_check("density integrates to 1", {
    total <- if (is_cont) {
      m <- suppressWarnings(distrib_quantile(distrib, 0.5, theta))
      knots <- if (is.finite(m) && m > b[1] && m < b[2]) c(b[1], m, b[2]) else b
      # Integrate up to each atom and away from it again: the panels stay on one
      # side of the jump, and the mass in the jump is added explicitly.
      knots <- sort(unique(c(knots, atoms$y[atoms$y > b[1] & atoms$y < b[2]])))
      atom_mass + sum(vapply(seq_len(length(knots) - 1L), function(k) {
        stats::integrate(function(t) distrib_pdf(distrib, t, theta),
                         knots[k], knots[k + 1L])$value
      }, numeric(1)))
    } else {
      numerical_series(function(k) distrib_pdf(distrib, k, theta), b[1], b[2])
    }
    new_check("density integrates to 1", abs(total - 1) < 1e-5, abs(total - 1))
  })

  res[[length(res) + 1L]] <- safe_check("density is non-negative", {
    grid <- distrib_quantile(distrib, seq(0.01, 0.99, length.out = 25), theta)
    v <- distrib_pdf(distrib, grid, theta)
    new_check("density is non-negative", all(is.finite(v) & v >= 0), min(v))
  })

  # --- cdf -----------------------------------------------------------------
  res[[length(res) + 1L]] <- safe_check("cdf in [0,1] and non-decreasing", {
    grid <- distrib_quantile(distrib, seq(0.02, 0.98, length.out = 40), theta)
    Fv <- distrib_cdf(distrib, sort(grid), theta)
    new_check("cdf in [0,1] and non-decreasing",
              all(Fv >= -1e-10 & Fv <= 1 + 1e-10) && all(diff(Fv) >= -1e-10),
              min(diff(Fv)))
  })

  # Nothing above ties the cdf to the density, and the quantile round-trip cannot:
  # when the quantile function is the numerical fallback it is derived from the
  # cdf, so the two agree with each other however wrong the cdf is. A cdf shifted
  # by a constant passes every other check. Comparing it against the density is
  # what pins it down: F' = f for a continuous distribution, and F(k) - F(k-1) =
  # P(Y = k) for a discrete one.
  res[[length(res) + 1L]] <- safe_check("cdf agrees with the density", {
    if (is_cont) {
      grid <- distrib_quantile(distrib, seq(0.1, 0.9, length.out = 15), theta)
      grid <- unique(grid[away_from_atoms(grid)])
      # Everything below the highest atom collapsed onto it; take the grid from
      # the continuous stretch above instead of reporting a failure.
      if (length(grid) < 3 && length(atoms$y)) {
        p_hi <- max(distrib_cdf(distrib, atoms$y, theta))
        grid <- distrib_quantile(distrib, seq(p_hi + (1 - p_hi) * 0.1, 0.9, length.out = 15), theta)
        grid <- unique(grid[away_from_atoms(grid)])
      }
      if (length(grid) < 3) {
        stop("the atoms leave too little of the support to compare F' against f.",
             call. = FALSE)
      }
      h <- pmax(abs(grid), 1) * 1e-5
      d_num <- (distrib_cdf(distrib, grid + h, theta) -
                distrib_cdf(distrib, grid - h, theta)) / (2 * h)
      err <- rel(d_num, distrib_pdf(distrib, grid, theta))
    } else {
      ks <- distrib_quantile(distrib, seq(0.1, 0.9, length.out = 9), theta)
      ks <- unique(ks[ks > b[1]])
      err <- rel(distrib_cdf(distrib, ks, theta) - distrib_cdf(distrib, ks - 1, theta),
                 distrib_pdf(distrib, ks, theta))
    }
    new_check("cdf agrees with the density", err < 1e-4, err)
  })

  # --- quantile ------------------------------------------------------------
  res[[length(res) + 1L]] <- safe_check("quantile/cdf round-trip", {
    p <- c(0.05, 0.25, 0.5, 0.75, 0.95)
    q <- distrib_quantile(distrib, p, theta)
    if (is_cont) {
      # F(Q(p)) = p only where F is continuous. A probability falling inside a
      # jump maps to the atom, and F there is the top of the jump; the
      # generalized inverse is still right, so those p are checked the discrete
      # way and dropped from the equality.
      inside <- !away_from_atoms(q)
      ok_jump <- !any(inside) ||
        all(distrib_cdf(distrib, q[inside], theta) >= p[inside] - 1e-10)
      err <- if (all(inside)) 0 else max(abs(distrib_cdf(distrib, q[!inside], theta) - p[!inside]))
      new_check("quantile/cdf round-trip", err < 1e-5 && ok_jump, err)
    } else {
      ok <- all(distrib_cdf(distrib, q, theta) >= p - 1e-10) &&
        all(distrib_cdf(distrib, q - 1, theta) < p + 1e-10)
      new_check("quantile/cdf round-trip", ok, NA_real_)
    }
  })

  # --- rng -----------------------------------------------------------------
  # Compared against the cdf rather than against the first two moments: a mean
  # and a variance are not available for every distribution (the Cauchy has
  # neither), whereas the cdf always is, and matching it is the stronger claim.
  # The comparison is valid for discrete distributions too, since both sides are
  # P(Y <= q) at the same points.
  res[[length(res) + 1L]] <- safe_check("rng matches the cdf", {
    ys <- distrib_rng(distrib, nsim, theta)
    probs <- seq(0.05, 0.95, by = 0.05)
    q_th <- distrib_quantile(distrib, probs, theta)
    p_th <- distrib_cdf(distrib, q_th, theta)
    p_emp <- vapply(q_th, function(q) base::mean(ys <= q), numeric(1))
    se <- sqrt(pmax(p_th * (1 - p_th), .Machine$double.eps) / nsim)
    z <- max(abs(p_emp - p_th) / se)
    new_check("rng matches the cdf", z < 6, z)
  })

  # --- parameter-scale derivatives ----------------------------------------
  y <- distrib_rng(distrib, n, theta)

  if (1 %in% orders) {
    res[[length(res) + 1L]] <- safe_check("gradient vs finite differences", {
      a <- distrib_gradient(distrib, y, theta)
      e <- numerical_gradient(distrib, y, theta)
      keep <- fd_is_reliable(function(h) numerical_gradient(distrib, y, theta, h_rel = h),
                             e, .Machine$double.eps^(1 / 3), smooth_all)
      err <- max(vapply(names(a), function(k) rel(a[[k]][keep], e[[k]][keep]), numeric(1)))
      new_check("gradient vs finite differences", err < tol, err)
    })
  }
  if (2 %in% orders) {
    res[[length(res) + 1L]] <- safe_check("hessian vs finite differences", {
      a <- distrib_hessian(distrib, y, theta)
      e <- numerical_hessian(distrib, y, theta)
      keep <- fd_is_reliable(function(h) numerical_hessian(distrib, y, theta, h_rel = h),
                             e, .Machine$double.eps^(1 / 4), smooth_all)
      err <- max(vapply(names(a), function(k) rel(a[[k]][keep], e[[k]][keep]), numeric(1)))
      new_check("hessian vs finite differences", err < tol, err)
    })
  }
  if (3 %in% orders) {
    res[[length(res) + 1L]] <- safe_check("deriv3 vs finite differences", {
      a <- distrib_deriv3(distrib, y, theta)
      e <- numerical_deriv3(distrib, y, theta)
      err <- max(vapply(names(a), function(k) rel(a[[k]], e[[k]]), numeric(1)))
      new_check("deriv3 vs finite differences", err < tol, err)
    })
  }
  if (4 %in% orders) {
    res[[length(res) + 1L]] <- safe_check("deriv4 vs finite differences", {
      a <- distrib_deriv4(distrib, y, theta)
      e <- numerical_deriv4(distrib, y, theta)
      err <- max(vapply(names(a), function(k) rel(a[[k]], e[[k]]), numeric(1)))
      new_check("deriv4 vs finite differences", err < tol, err)
    })
  }

  # --- expected information (outer product of the score) -------------------
  res[[length(res) + 1L]] <- safe_check("expected information vs Monte Carlo", {
    ys <- distrib_rng(distrib, nsim, theta)
    g <- distrib_gradient(distrib, ys, theta)
    a <- distrib_expected_hessian(distrib, 0, theta)
    params <- distrib@params
    worst <- 0
    for (i in seq_along(params)) {
      for (j in i:length(params)) {
        nm <- paste(params[c(i, j)], collapse = "_")
        prod_ij <- g[[i]] * g[[j]]
        mc <- -base::mean(prod_ij)
        se <- stats::sd(prod_ij) / sqrt(nsim)
        # The score product is not always random: for a non-smooth location
        # parameter it is constant across draws (the Laplace score is
        # sign(y-mu)/b, so its square is 1/b^2 everywhere), leaving a standard
        # error that is pure floating-point dust. Standardising by it turns a
        # perfect agreement into an enormous z. Floor the denominator at the
        # precision the two sides can possibly agree to, which keeps a genuinely
        # wrong analytic information easy to detect.
        den <- max(se, 1e-8 * max(1, abs(mc)))
        worst <- max(worst, abs(a[[nm]][1] - mc) / den)
      }
    }
    new_check("expected information vs Monte Carlo", worst < 6, worst)
  })

  # --- response derivatives (continuous only) ------------------------------
  if (is_cont) {
    res[[length(res) + 1L]] <- safe_check("response derivatives vs finite differences", {
      # The log-density jumps at an atom, so no derivative in y exists there and
      # the finite-difference reference is meaningless: those draws are dropped.
      y <- y[away_from_atoms(y)]
      if (length(y) < 5) {
        stop("too few draws land away from the atoms to check the response derivatives.",
             call. = FALSE)
      }
      # Differencing in y crosses the same kink, so the same guard applies.
      g_ref <- numerical_grad_y(distrib, y, theta)
      k1 <- fd_is_reliable(function(h) list(v = numerical_grad_y(distrib, y, theta, h_rel = h)),
                           list(v = g_ref), .Machine$double.eps^(1 / 3), smooth_all)
      h_ref <- numerical_hess_y(distrib, y, theta)
      k2 <- fd_is_reliable(function(h) list(v = numerical_hess_y(distrib, y, theta, h_rel = h)),
                           list(v = h_ref), .Machine$double.eps^(1 / 4), smooth_all)
      e1 <- rel(distrib_grad_y(distrib, y, theta)[k1], g_ref[k1])
      e2 <- rel(distrib_hess_y(distrib, y, theta)[k2], h_ref[k2])
      err <- max(e1, e2)
      new_check("response derivatives vs finite differences", err < tol, err)
    })
  }

  # --- link scale ----------------------------------------------------------
  res[[length(res) + 1L]] <- safe_check("link-scale gradient vs finite differences", {
    params <- distrib@params
    eta <- vapply(seq_along(params), function(i) {
      linkfunctions7::linkfun(distrib@link_params[[params[i]]], theta[[i]])
    }, numeric(1))
    theta_of <- function(e) {
      out <- lapply(seq_along(params), function(i) {
        linkfunctions7::linkinv(distrib@link_params[[params[i]]], e[i])
      })
      names(out) <- params
      out
    }
    a <- distrib_gradient(distrib, y, theta, scale = "link")
    fd_eta <- function(i, h) {
      ep <- em <- eta
      ep[i] <- ep[i] + h
      em[i] <- em[i] - h
      (distrib_pdf(distrib, y, theta_of(ep), log = TRUE) -
         distrib_pdf(distrib, y, theta_of(em), log = TRUE)) / (2 * h)
    }
    h <- 1e-4
    worst <- 0
    for (i in seq_along(params)) {
      num <- fd_eta(i, h)
      # Same guard as on the parameter scale: a step in eta_mu crosses the kink.
      keep <- fd_is_reliable(function(hh) list(v = fd_eta(i, hh)),
                             list(v = num), h, smooth_all)
      worst <- max(worst, rel(a[[i]][keep], num[keep]))
    }
    new_check("link-scale gradient vs finite differences", worst < tol, worst)
  })

  out <- do.call(rbind, res)
  rownames(out) <- NULL

  if (verbose) {
    cat("Distribution: ", distrib@distrib_name, "\n", sep = "")
    cat("Parameters:   ",
        paste(names(theta), vapply(theta, function(v) format(v[1], digits = 4), character(1)),
              sep = " = ", collapse = ", "), "\n", sep = "")
    cat("Observations: ", n, "   Monte Carlo: ", format(nsim, scientific = FALSE), "\n\n", sep = "")
    width <- max(nchar(out$check))
    for (i in seq_len(nrow(out))) {
      stat <- if (is.na(out$statistic[i])) "" else sprintf("  %.2e", out$statistic[i])
      cat(sprintf("  [%-4s] %-*s%s\n", out$status[i], width, out$check[i], stat))
      if (!is.na(out$detail[i])) cat("         ", out$detail[i], "\n", sep = "")
    }
    n_fail <- sum(out$status == "FAIL")
    cat("\n", if (n_fail == 0) {
      sprintf("All %d checks passed.\n", nrow(out))
    } else {
      sprintf("%d of %d checks FAILED.\n", n_fail, nrow(out))
    }, sep = "")
  }

  invisible(out)
}

#' Which Observations the Finite-Difference Reference Can Be Trusted At
#'
#' @description
#' Flags the observations where a central difference has actually converged, so
#' that \code{\link{check_distrib}} compares an analytical derivative only
#' against a reference worth comparing to.
#'
#' @details
#' A log-likelihood with a kink -- the Laplace's location is the example the
#' package ships -- has no derivative exactly at the kink, and a central
#' difference straddling it returns a number that is simply wrong. An observation
#' landing within a step of that point therefore makes the \emph{reference}
#' invalid, not the analytical value being tested, and comparing against it
#' reports a failure for code that is right. Because the draws are random, that
#' happened rarely and unpredictably.
#'
#' Rather than hard-code where a kink is, the reference is recomputed with the
#' step halved: where the two disagree, finite differencing has not converged and
#' that observation is dropped. For a smooth log-likelihood nothing is ever
#' dropped, so the check keeps its full strength -- a gradient made 5\% wrong is
#' still caught.
#'
#' Two details are load-bearing. The two estimates are compared \strong{relative
#' to their own magnitude} rather than against a denominator floored at one: near
#' a kink both are tiny yet differ by a factor of two, which a floor of one
#' flattens into apparent agreement. And if no observation survives, all are kept:
#' a systematic disagreement is a real failure and should be reported, not hidden
#' by the guard meant to protect against a local one.
#'
#' Note the placement. This is defined \emph{after} \code{check_distrib()}, not
#' before it: a roxygen block attaches to whatever object follows it, so a helper
#' slipped in between silently steals the documentation of the function it
#' belongs to -- which had already happened once here, leaving the package with a
#' \code{fd_is_reliable.Rd} and no \code{check_distrib.Rd}.
#'
#' @param fd_at A function of one argument, the relative step, returning the
#'   finite-difference reference computed with that step.
#' @param ref The reference already computed at \code{h_rel}, a named list of
#'   component vectors.
#' @param h_rel The relative step \code{ref} was computed with.
#' @param smooth_all Logical; \code{TRUE} when every parameter is declared
#'   smooth, in which case no observation is ever dropped and the guard does not
#'   run at all.
#'
#' @return A logical vector as long as the components of \code{ref}.
#'
#' @seealso \code{\link{check_distrib}}
#' @keywords internal
fd_is_reliable <- function(fd_at, ref, h_rel, smooth_all) {
  n <- length(ref[[1L]])
  if (isTRUE(smooth_all)) return(rep(TRUE, n))

  half <- tryCatch(fd_at(h_rel / 2), error = function(e) NULL)
  if (is.null(half)) return(rep(TRUE, n))

  ok <- rep(TRUE, n)
  for (k in names(ref)) {
    # Measured against the estimates' own magnitude, not floored at one: near a
    # kink both are small yet differ by a factor of two, which a denominator of
    # one would flatten into apparent agreement.
    d <- abs(ref[[k]] - half[[k]]) / pmax(abs(ref[[k]]), abs(half[[k]]), 1e-8)
    ok <- ok & (is.finite(d) & d < 1e-3)
  }
  # Never discard everything: if no observation survives, the disagreement is
  # systematic and should be reported rather than hidden.
  if (!any(ok)) rep(TRUE, n) else ok
}
