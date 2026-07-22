#' @include distrib.R generics.R utility_functions.R numerical_derivatives.R higher_derivatives.R y_derivatives.R
NULL

# Internal: accumulate one check result.
new_check <- function(name, ok, stat, detail = NA_character_) {
  data.frame(
    check = name, status = if (isTRUE(ok)) "OK" else "FAIL",
    statistic = stat, detail = detail, stringsAsFactors = FALSE
  )
}

# Internal: run an expression, turning an error into a failed check.
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

  rel <- function(a, e) max(abs(a - e) / pmax(1, abs(e)))

  # --- density -------------------------------------------------------------
  res[[length(res) + 1L]] <- safe_check("density integrates to 1", {
    total <- if (is_cont) {
      m <- suppressWarnings(distrib_quantile(distrib, 0.5, theta))
      knots <- if (is.finite(m) && m > b[1] && m < b[2]) c(b[1], m, b[2]) else b
      sum(vapply(seq_len(length(knots) - 1L), function(k) {
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

  # --- quantile ------------------------------------------------------------
  res[[length(res) + 1L]] <- safe_check("quantile/cdf round-trip", {
    p <- c(0.05, 0.25, 0.5, 0.75, 0.95)
    q <- distrib_quantile(distrib, p, theta)
    if (is_cont) {
      err <- max(abs(distrib_cdf(distrib, q, theta) - p))
      new_check("quantile/cdf round-trip", err < 1e-5, err)
    } else {
      ok <- all(distrib_cdf(distrib, q, theta) >= p - 1e-10) &&
        all(distrib_cdf(distrib, q - 1, theta) < p + 1e-10)
      new_check("quantile/cdf round-trip", ok, NA_real_)
    }
  })

  # --- rng -----------------------------------------------------------------
  res[[length(res) + 1L]] <- safe_check("rng matches theoretical moments", {
    ys <- distrib_rng(distrib, nsim, theta)
    m_th <- mean(distrib, theta)
    v_th <- variance(distrib, theta)
    z_m <- abs(base::mean(ys) - m_th) / sqrt(v_th / nsim)
    new_check("rng matches theoretical moments", z_m < 6, z_m)
  })

  # --- parameter-scale derivatives ----------------------------------------
  set.seed(NULL)
  y <- distrib_rng(distrib, n, theta)

  if (1 %in% orders) {
    res[[length(res) + 1L]] <- safe_check("gradient vs finite differences", {
      a <- distrib_gradient(distrib, y, theta)
      e <- numerical_gradient(distrib, y, theta)
      err <- max(vapply(names(a), function(k) rel(a[[k]], e[[k]]), numeric(1)))
      new_check("gradient vs finite differences", err < tol, err)
    })
  }
  if (2 %in% orders) {
    res[[length(res) + 1L]] <- safe_check("hessian vs finite differences", {
      a <- distrib_hessian(distrib, y, theta)
      e <- numerical_hessian(distrib, y, theta)
      err <- max(vapply(names(a), function(k) rel(a[[k]], e[[k]]), numeric(1)))
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
      e1 <- rel(distrib_grad_y(distrib, y, theta), numerical_grad_y(distrib, y, theta))
      e2 <- rel(distrib_hess_y(distrib, y, theta), numerical_hess_y(distrib, y, theta))
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
    h <- 1e-4
    worst <- 0
    for (i in seq_along(params)) {
      ep <- em <- eta
      ep[i] <- ep[i] + h
      em[i] <- em[i] - h
      num <- (distrib_pdf(distrib, y, theta_of(ep), log = TRUE) -
                distrib_pdf(distrib, y, theta_of(em), log = TRUE)) / (2 * h)
      worst <- max(worst, rel(a[[i]], num))
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
