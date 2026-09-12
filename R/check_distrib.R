#' @include distrib.R generics.R utility_functions.R numerical_derivatives.R higher_derivatives.R y_derivatives.R
NULL

#' Record One Check Result
#'
#' @description
#' Builds the single-row data frame that [check_distrib()] accumulates
#' into its report.
#'
#' @param name The check's name, as it appears in the report.
#' @param ok Logical; whether the check passed.
#' @param stat The numeric statistic the check produced, typically a maximum
#'   discrepancy.
#' @param detail An optional message, used to carry the reason for a failure.
#'
#' @return A one-row data frame with columns `check`, `status`,
#'   `statistic` and `detail`.
#'
#' @seealso [check_distrib()], [safe_check()]
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
#' would end the run and hide every check after it, which is the least useful
#' moment to
#' stop being informative.
#'
#' The failed row is marked with the attribute `from_error`, which
#' [check_distrib()] reads before assembling the table: a check whose
#' computation raised has failed, and it is never mistaken for one whose
#' statistic merely came out non-finite, which is reported as not run.
#'
#' @param name The check's name, used for the row built on failure.
#' @param expr The expression to evaluate; normally returns a row from
#'   [new_check()].
#'
#' @return The value of `expr`, or a failed row carrying the error message and
#'   the attribute `from_error = TRUE`.
#'
#' @seealso [check_distrib()], [new_check()]
#' @keywords internal
safe_check <- function(name, expr) {
  tryCatch(expr, error = function(e) {
    row <- new_check(name, FALSE, NA_real_, conditionMessage(e))
    attr(row, "from_error") <- TRUE
    row
  })
}

#' Numerically Validate a Distribution
#'
#' @description
#' Runs a battery of numerical self-consistency checks on a `distrib` object.
#' Validates a user-defined distribution: it verifies that the density integrates (or sums) to one, that the
#' CDF, quantile function and random generator agree with each other and with the
#' density, and that every analytical derivative matches its finite-difference
#' counterpart.
#'
#' @param distrib An object inheriting from class `"distrib"`.
#' @param theta A named list of parameter values at which to run the checks. If
#'   `NULL` (default) a random admissible value is drawn with
#'   [generate_random_theta()].
#' @param n Integer. Number of observations used for the derivative comparisons.
#'   Defaults to 100.
#' @param nsim Integer. Monte Carlo sample size used for the random generator and
#'   expected-information checks. Defaults to 200000.
#' @param orders Integer vector. Which parameter-derivative orders to check.
#'   Defaults to `1:4`, the orders every family implements analytically; use
#'   e.g. `1:2` for a faster run. Add `5` to check the numerical fifth as
#'   well; see the bullet below for what that row compares and when it is
#'   emitted.
#' @param tol Numeric. Relative tolerance for the finite-difference comparisons.
#'   Defaults to `1e-3`.
#' @param verbose Logical. If `TRUE` (default) a readable report is printed.
#'
#' @return Invisibly, a `data.frame` with one row per check and columns
#'   `check`, `status` (`"OK"` or `"FAIL"`), `statistic`
#'   and `detail`. A check whose statistic came out `NaN` or `NA` is not a row:
#'   it is listed in the attribute `"skipped"`, a data frame with columns
#'   `check` and `reason`, which is present only when some check was skipped.
#'
#' @details
#' The checks performed are:
#'
#' - **density**: non-negativity and integration/summation to 1 over the support.
#' - **cdf**: values in \eqn{[0,1]} and monotonicity along a grid of quantiles,
#'   and its agreement with the density: for a continuous family a central
#'   difference of the cdf against the density, with the step the response
#'   derivatives use, and for a discrete one \eqn{F(k) - F(k-1)} against the
#'   mass.
#' - **quantile**: round-trip against the CDF (\eqn{F(Q(p)) = p} for continuous
#'   distributions, and the generalized-inverse inequalities for discrete ones).
#' - **rng**: the sample mean and variance of a large draw agree with
#'   [mean()] and [variance()] within Monte Carlo error.
#' - **gradient, hessian, deriv3, deriv4**: analytical values against
#'   [numerical_gradient()], [numerical_hessian()],
#'   [numerical_deriv3()] and [numerical_deriv4()].
#' - **deriv5**, when `5` is among `orders`: [distrib_deriv5()] against
#'   [numerical_deriv5()] at a higher accuracy, five stencil nodes instead of
#'   three. No family writes the fifth order out, so there is no analytic value
#'   to compare against and this checks the differencing rather than a family's
#'   algebra. It is emitted only where [has_exact_deriv4()] is `TRUE`: where the
#'   fourth order is itself a fallback the fifth is a difference of a difference
#'   and no verdict on it would mean anything. A family that owns its
#'   fourth-order method while building part of it from stencils passes that
#'   test and may still fail this row: [skewt_distrib()] fails it at `nu = 3`
#'   and passes from `nu = 8` upward, the noise read here being absolute and
#'   falling as `nu` grows. It is why `orders` defaults to `1:4`.
#' - **expected information**: [distrib_expected_hessian()] against a
#'   Monte Carlo estimate of \eqn{-\mathbb{E}[\nabla\ell\,\nabla\ell^\top]}. The outer
#'   product of the score is used as reference because it remains valid when the
#'   log-likelihood is not differentiable in a parameter (see [laplace_distrib()]).
#'   A draw landing exactly on a finite bound of a continuous support, where the
#'   score of a family singular at that bound is infinite, is left out of the
#'   estimate and counted in the row's `detail`; a score that is not finite at
#'   an interior point still leaves the row without a statistic.
#' - **response derivatives** (continuous only): [distrib_grad_y()] and
#'   [distrib_hess_y()] against central differences in \eqn{y} whose step is
#'   chosen observation by observation by [fd_stable_quotient()], between a
#'   step cut to under half the distance to a bound and one scaled on that
#'   distance, each divided by the steps its evaluation points actually lie
#'   at. A draw closer than \eqn{128\lvert b\rvert\varepsilon} to a finite
#'   bound \eqn{b \ne 0}, where the spacing of doubles is absolute and no
#'   central difference compares a derivative, is left out and counted in the
#'   row's `detail`.
#' - **link scale**: `scale = "link"` derivatives against finite
#'   differences of the log-likelihood in \eqn{\eta}.
#'
#' Distributions that rely on the numerical fallbacks pass the corresponding
#' parameter-derivative checks trivially, since analytical and numerical values
#' then coincide by construction. The response fallbacks take the cut step
#' alone, so near a bound they can differ from the reference, which is what the
#' row then reports.
#'
#' A check whose statistic comes out `NaN` or `NA` has nothing to judge, as the
#' expected information of a family where it does not exist, and is not a row
#' of the table: it is listed with its reason in the attribute `"skipped"`, and
#' every row that is in the table reads `"OK"` or `"FAIL"`. An infinite
#' statistic is a failure, being a component that overflows where its reference
#' does not; so is a check whose computation raised an error, and so is a
#' density or a distribution function that is not finite, those three checks
#' being defined on the values themselves.
#'
#' Mixed distributions --- a density with point masses on top of it, as produced by
#' [zero_adjusted()] on a continuous parent --- are handled as long as they
#' declare their atoms through [distrib_atoms()]. The density is then expected
#' to integrate to one minus the atomic mass, quantiles falling inside a jump of the CDF
#' are checked as generalized inverses rather than exact ones, and finite differences in
#' \eqn{y} are kept away from the atoms, where no derivative exists.
#'
#' @examples
#' \dontrun{
#' check_distrib(gaussian1_distrib())
#' check_distrib(laplace_distrib(), theta = list(mu = 1, sigma = 2))
#' check_distrib(poisson_distrib(), orders = 1:2, nsim = 5e4)
#' }
#'
#' @seealso [numerical_gradient()], [numerical_hessian()],
#'   [link_scale_derivatives()]
#' @export
check_distrib <- function(distrib, theta = NULL, n = 100, nsim = 2e5,
                          orders = 1:4, tol = 1e-3, verbose = TRUE) {
  if (is.null(theta)) theta <- generate_random_theta(distrib)
  theta <- align_theta(distrib, theta)

  # A multivariate distribution has no distribution function, no quantile
  # function and no one-dimensional quadrature, so five of the checks below
  # have no counterpart. Running them anyway reports refusals as failures,
  # which is the mistake this function already learned not to make with mixed
  # distributions: a user validating their own code could not tell a real
  # defect from it. The battery that does generalize is a separate one.
  if (S7::S7_inherits(distrib, multivariate_distrib)) {
    out <- do.call(rbind, check_distrib_mv(distrib, theta, n, nsim, tol))
    if (verbose) print_check_table(distrib, out, theta, n, nsim)
    return(invisible(out))
  }

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
      discrete_support_sum(function(k, i) distrib_pdf(distrib, k, theta),
                           b[1], b[2], 1L)
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
      # A GRID POINT IN THE LAST PLACES OF A NON-ZERO BOUND IS LEFT OUT AND
      # COUNTED, the rule the response row already applies to a draw, for the
      # same reason: near zero the spacing of doubles is relative, near any
      # other bound it is absolute, and within a few hundred units of it the
      # reference compares nothing. This grid is the deciles, which on almost
      # every family sit nowhere near a bound -- measured, beta2 is 4.9e+04
      # units away at worst and the two von Mises 3.2e+14. beta1 is the
      # exception, because it carries shape2 = (1 - mu) phi, whose box reaches
      # 0.1 x 0.1 = 0.01, and its upper decile is then 1 - O(exp(-c/shape2)),
      # which stops being representable below about shape2 = 0.06: over 20000
      # draws of generate_random_theta() the row failed 59 times, with
      # statistics up to 8.006 on analytic code, every one of them between 0.5
      # and 85 units of the bound, and a further 185 draws had the grid
      # collapse onto the bound outright.
      near <- rep(FALSE, length(grid))
      if (is.finite(b[1]) && b[1] != 0) {
        near <- near | (grid - b[1] < 128 * abs(b[1]) * .Machine$double.eps)
      }
      if (is.finite(b[2]) && b[2] != 0) {
        near <- near | (b[2] - grid < 128 * abs(b[2]) * .Machine$double.eps)
      }
      n_near <- sum(near)
      n_all <- length(grid)
      if (n_near) grid <- grid[!near]
      if (length(grid) < 3) {
        stop("too much of the grid lies within 128 |b| eps of a non-zero bound to compare F' against f.",
             call. = FALSE)
      }
      # THE SAME STEP AS THE RESPONSE ROW. It was 1e-5 max(1, |x|), not kept
      # inside the support, and it failed where the grid comes close to a
      # bound: measured over this grid for 32 continuous families and three
      # parameter values, 3 rows (gamma1 twice and gamma2 once, the worst at
      # 1.1e-02 on the grid point nearest the bound), and on a sweep toward the
      # bounds 66 of 144 points of the families smooth at the bound and 85 of
      # 264 of the singular ones. fd_stable_quotient() fails none of either,
      # its worst at 9.2e-05. Dividing by the steps actually taken moves no
      # statistic of the census's grid by more than 1.4e-10; near a non-zero
      # bound it moves more, beta1 at mu = 0.5 and phi = 0.5 reading 1.7e-10
      # where the nominal quotient read 5.3e-09.
      d_num <- fd_stable_quotient(
        function(h) fd_first_taken(function(v) distrib_cdf(distrib, v, theta), grid, h),
        grid, b, .Machine$double.eps^(1 / 3))
      err <- rel(d_num, distrib_pdf(distrib, grid, theta))
    } else {
      ks <- distrib_quantile(distrib, seq(0.1, 0.9, length.out = 9), theta)
      ks <- unique(ks[ks > b[1]])
      err <- rel(distrib_cdf(distrib, ks, theta) - distrib_cdf(distrib, ks - 1, theta),
                 distrib_pdf(distrib, ks, theta))
      # a discrete family compares an exact difference of the distribution
      # function against the mass, with no step, so a point on a bound costs
      # it nothing and nothing is left out
      n_near <- 0L
      n_all <- length(ks)
    }
    new_check("cdf agrees with the density", err < 1e-4, err,
              if (n_near) {
                sprintf("%d of %d grid points lay within 128 |b| eps of a non-zero bound and were left out",
                        n_near, n_all)
              } else NA_character_)
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
  # The fifth order has no analytic implementation to compare against: every
  # family reaches numerical_deriv5(), one central difference of the fourth.
  # What can be checked is that difference against a HIGHER-ACCURACY rule from
  # the same library -- five nodes instead of three, different weights and a
  # different step -- which catches a mis-keyed component, a wrong grouping or
  # a wrong scale, and is the only in-package reference there is. The
  # independent one, Richardson on the analytic fourth, lives in the tests,
  # where numDeriv may be named.
  #
  # The row is emitted only where the fourth order is the family's own. Where
  # it is not, the fifth is a difference of a difference -- measured on a
  # density-only gaussian it is 3e+04 relative, which is not a derivative of
  # anything -- and reporting a verdict on it would be reporting a pass nobody
  # earned. This follows the convention the multivariate battery already uses:
  # a check that does not apply is not emitted, rather than emitted with a
  # third status that every consumer reading `status != "OK"` would misread.
  #
  # A family may own that method and still build part of it from single
  # stencils, in which case the row IS emitted and its verdict is a
  # measurement rather than a promise. Note that rel() above floors the
  # denominator at 1, so what this row reads for a small component is an
  # ABSOLUTE error: the skew t reads 3.7e-03 to 4.6e-03 at nu = 3, where it
  # fails, and 4.4e-05 to 4.3e-04 at nu = 8, where it passes. `orders` defaults
  # to 1:4, so the row is reached only by a caller who asked for it.
  if (5 %in% orders && has_exact_deriv4(distrib)) {
    res[[length(res) + 1L]] <- safe_check("deriv5 vs a higher-accuracy rule", {
      a <- distrib_deriv5(distrib, y, theta)
      e <- numerical_deriv5(distrib, y, theta, accuracy = 4L)
      err <- max(vapply(names(a), function(k) rel(a[[k]], e[[k]]), numeric(1)))
      new_check("deriv5 vs a higher-accuracy rule", err < tol, err)
    })
  }

  # --- expected information (outer product of the score) -------------------
  res[[length(res) + 1L]] <- safe_check("expected information vs Monte Carlo", {
    ys <- distrib_rng(distrib, nsim, theta)
    # A draw EXACTLY on a finite bound of a continuous support is left out and
    # counted. There a log-density that is logarithmic at the bound is not
    # finite, nor is the score, and one such draw turned the whole Monte Carlo
    # mean into NaN: measured on beta1 with nsim = 2e5, 5 to 11 draws per seed
    # land exactly on 1 in double precision and the row came back NaN in 3
    # seeds of 3, where without them z reads 0.51, 0.75 and 1.9. A point mass
    # the family declares through distrib_atoms() is not such a draw and stays
    # in. A score that is not finite anywhere else is not left out, and the
    # statistic then stays non-finite.
    on_bound <- rep(FALSE, length(ys))
    if (is_cont) {
      if (is.finite(b[1])) on_bound <- on_bound | ys == b[1]
      if (is.finite(b[2])) on_bound <- on_bound | ys == b[2]
      on_bound <- on_bound & !(ys %in% atoms$y)
    }
    n_bound <- sum(on_bound)
    if (n_bound) ys <- ys[!on_bound]
    g <- distrib_gradient(distrib, ys, theta)
    # approx = "bartlett" IS THE POINT HERE, and must be asked for rather than
    # left to the default. This check wants the actual EXPECTATION, to compare
    # against an independent Monte Carlo estimate of that same expectation;
    # "opg" (the default since 0.44.0) reads the per-observation outer product
    # of the score at the single point y = 0 instead of averaging it, which is
    # a different quantity and disagrees with a Monte Carlo mean almost always.
    a <- distrib_expected_hessian(distrib, 0, theta, approx = "bartlett")
    params <- distrib@params
    worst <- 0
    for (i in seq_along(params)) {
      for (j in i:length(params)) {
        nm <- paste(params[c(i, j)], collapse = "_")
        prod_ij <- g[[i]] * g[[j]]
        mc <- -base::mean(prod_ij)
        se <- stats::sd(prod_ij) / sqrt(length(ys))
        # The score product is not always random: for a non-smooth location
        # parameter it is constant across draws (the Laplace score is
        # sign(y-mu)/b, so its square is 1/b^2 everywhere), leaving a standard
        # error that is pure floating-point dust. Standardizing by it turns a
        # perfect agreement into an enormous z. Floor the denominator at the
        # precision the two sides can possibly agree to, which keeps a genuinely
        # wrong analytic information easy to detect.
        den <- max(se, 1e-8 * max(1, abs(mc)))
        worst <- max(worst, abs(a[[nm]][1] - mc) / den)
      }
    }
    detail <- character(0)
    if (!all(vapply(a, function(v) all(is.finite(v)), logical(1)))) {
      detail <- "the family's expected information is not finite at these parameters"
    }
    if (n_bound) {
      detail <- c(detail, sprintf("%d of %d draws fell exactly on a bound of the support and were left out",
                                  n_bound, nsim))
    }
    new_check("expected information vs Monte Carlo", worst < 6, worst,
              if (length(detail)) paste(detail, collapse = "; ") else NA_character_)
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
      # A DRAW IN THE LAST PLACES OF A NON-ZERO BOUND IS LEFT OUT AND COUNTED,
      # as a draw exactly on a bound is in the Monte Carlo row. Near zero the
      # spacing of doubles is relative; near any other bound it is absolute,
      # the reference's shortest step is two units of it, and its relative
      # error falls as about 5.4/d^2 with d the distance counted in units:
      # measured on beta1 and beta2 singular at 1, 1.3e-03 on the first
      # derivative and 2.0e-03 on the second at 64 units, 8.1e-05 and 1.2e-04
      # at 256. The draws closer than 256 units, 128 |b| eps, compare nothing.
      # On beta1 at mu 0.9 and phi 1 they are 4.4 per cent of the draws, and
      # with them the row failed in 1285 samples of 2000 even on the steps
      # taken; without them, in none. The census below has no such draw.
      near <- rep(FALSE, length(y))
      if (is.finite(b[1]) && b[1] != 0) {
        near <- near | (y - b[1] < 128 * abs(b[1]) * .Machine$double.eps)
      }
      if (is.finite(b[2]) && b[2] != 0) {
        near <- near | (b[2] - y < 128 * abs(b[2]) * .Machine$double.eps)
      }
      n_near <- sum(near)
      n_away <- length(y)
      if (n_near) y <- y[!near]
      if (length(y) < 5) {
        stop("too few draws land away from a non-zero bound to check the response derivatives.",
             call. = FALSE)
      }
      # Differencing in y crosses the same kink, so the same guard applies --
      # and here it applies ALWAYS, not only when a parameter is declared
      # non-smooth. What makes a difference in the RESPONSE unreliable is the
      # support boundary, not the parameters: a gamma is smooth in both of
      # its, so `smooth_all` was TRUE and the guard returned early without
      # examining anything.
      #
      # THE REFERENCE CHOOSES ITS STEP. With the step cut to 0.49 of the
      # distance to a bound, the relative error of a central difference on a
      # log-density carrying (a-1) log y is (h/d)^2/3 once the cut binds, a
      # plateau of 9.4e-02 on the first derivative and 1.4e-01 on the second
      # that halving the step cannot expose, the two steps being cut alike.
      # Measured over 32 continuous families, three parameter values and five
      # seeds, that reference failed 35 rows of 480, all of families singular
      # at a bound (gamma1 14, gamma2 10, chisq 4, beta1 3, beta2 3,
      # lognormal1 1), on analytic derivatives that agree with Richardson to
      # 1e-10. A step scaled on the distance removes the plateau and is worse
      # where the log-density is smooth at the bound, so fd_stable_quotient()
      # reads both and keeps, per observation, the one that agrees with itself
      # at half its step: the same census fails no row, the worst at 3.0e-05.
      # The quotients divide by the steps actually taken, and no step is below
      # two units of the spacing at y (see fd_stable_quotient()): over the
      # census no verdict moves and the worst goes to 5.0e-05, while on beta1
      # at mu 0.5 and phi 0.5, singular at 1, the row failed in 486 samples of
      # 2000 on the nominal quotient, fails in 10 on the steps taken, and in
      # none with the draws above left out.
      lp <- function(v) distrib_pdf(distrib, v, theta, log = TRUE)
      l0 <- lp(y)
      grad_ref <- function(h_rel) {
        fd_stable_quotient(function(h) fd_first_taken(lp, y, h), y, b, h_rel)
      }
      hess_ref <- function(h_rel) {
        fd_stable_quotient(function(h) fd_second_taken(lp, y, h, l0), y, b, h_rel)
      }
      g_ref <- grad_ref(.Machine$double.eps^(1 / 3))
      k1 <- fd_is_reliable(function(h) list(v = grad_ref(h)),
                           list(v = g_ref), .Machine$double.eps^(1 / 3), FALSE)
      h_ref <- hess_ref(.Machine$double.eps^(1 / 4))
      k2 <- fd_is_reliable(function(h) list(v = hess_ref(h)),
                           list(v = h_ref), .Machine$double.eps^(1 / 4), FALSE)
      e1 <- rel(distrib_grad_y(distrib, y, theta)[k1], g_ref[k1])
      e2 <- rel(distrib_hess_y(distrib, y, theta)[k2], h_ref[k2])
      err <- max(e1, e2)
      new_check("response derivatives vs finite differences", err < tol, err,
                if (n_near) {
                  sprintf("%d of %d draws lay within 128 |b| eps of a non-zero bound and were left out",
                          n_near, n_away)
                } else NA_character_)
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

  # A CHECK WITH NO VERDICT IS NOT EMITTED. A row whose statistic came out NaN
  # or NA without its computation raising -- a reference that is not finite
  # somewhere, a sample the statistic cannot be formed from, an expected
  # information the family itself reports as not existing -- compares nothing,
  # and reporting it as FAIL made a correct family fail: the generalized Pareto
  # at xi < -1/2, where the information does not exist, failed in 5 runs of 5.
  # It is left out of the table and recorded with its reason in the attribute
  # "skipped", so every row that is there says "OK" or "FAIL". Only NaN and NA
  # are read that way: an INFINITE statistic is a component that overflows
  # where its reference does not, which is a defect this check exists to catch
  # (gamma1's phi_phi read -Inf against a finite difference until distributions7
  # 0.54.0), so it stays a row and fails. Three rows are excepted, their
  # criterion being defined on the values themselves: a density or a cdf that
  # is not finite fails its own check, and the discrete round-trip carries NA
  # by construction. A row whose computation raised keeps its FAIL, whatever
  # its statistic.
  on_values <- c("density is non-negative", "cdf in [0,1] and non-decreasing",
                 if (!is_cont) "quantile/cdf round-trip")
  skipped <- list()
  kept <- vapply(res, function(r) {
    if (isTRUE(attr(r, "from_error")) || r$check %in% on_values ||
        !is.na(r$statistic)) {
      return(TRUE)
    }
    skipped[[length(skipped) + 1L]] <<- data.frame(
      check = r$check,
      reason = paste(c(sprintf("the statistic is %s", format(r$statistic)),
                       if (!is.na(r$detail)) r$detail), collapse = "; "),
      stringsAsFactors = FALSE)
    FALSE
  }, logical(1))
  out <- if (any(kept)) do.call(rbind, res[kept]) else res[[1L]][0L, , drop = FALSE]
  rownames(out) <- NULL
  if (length(skipped)) {
    sk <- do.call(rbind, skipped)
    rownames(sk) <- NULL
    attr(out, "skipped") <- sk
  }

  if (verbose) print_check_table(distrib, out, theta, n, nsim)

  invisible(out)
}

#' Print a Validation Table
#'
#' @description
#' Renders what [check_distrib()] found, in the same shape for a
#' univariate and a multivariate distribution.
#'
#' @param distrib An object inheriting from class `"distrib"`.
#' @param out The data frame of checks.
#' @param theta The parameters the checks were run at.
#' @param n The number of observations used.
#' @param nsim The Monte Carlo sample size used.
#'
#' @return Invisibly `NULL`.
#'
#' @seealso [check_distrib()]
#' @keywords internal
print_check_table <- function(distrib, out, theta, n, nsim) {
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
  sk <- attr(out, "skipped")
  if (!is.null(sk) && nrow(sk)) {
    wsk <- max(width, nchar(sk$check))
    for (i in seq_len(nrow(sk))) {
      cat(sprintf("  [ -- ] %-*s  not run\n", wsk, sk$check[i]))
      cat("         ", sk$reason[i], "\n", sep = "")
    }
  }
  n_fail <- sum(out$status == "FAIL")
  n_skip <- if (is.null(sk)) 0L else nrow(sk)
  cat("\n", if (n_fail == 0) {
    sprintf("All %d checks passed.\n", nrow(out))
  } else {
    sprintf("%d of %d checks FAILED.\n", n_fail, nrow(out))
  }, sep = "")
  if (n_skip) {
    cat(sprintf("%d check%s not run, the statistic having no value to judge.\n",
                n_skip, if (n_skip == 1L) "" else "s"))
  }
  invisible(NULL)
}

#' Which Observations the Finite-Difference Reference Can Be Trusted At
#'
#' @description
#' Flags the observations where a central difference has actually converged, so
#' that [check_distrib()] compares an analytical derivative only
#' against a reference that is itself reliable.
#'
#' @details
#' A log-likelihood with a kink has no derivative exactly at the kink, the
#' Laplace's location being the example the package ships, and a central
#' difference straddling it returns a number that is simply wrong. An observation
#' landing within a step of that point therefore makes the *reference*
#' invalid, not the analytical value being tested, and comparing against it
#' reports a failure for code that is right. Because the draws are random, that
#' happened rarely and unpredictably.
#'
#' Rather than hard-code where a kink is, the reference is recomputed with the
#' step halved: where the two disagree, finite differencing has not converged and
#' that observation is dropped. For a smooth log-likelihood nothing is ever
#' dropped, so the check keeps its full strength: a gradient made 5\% wrong is
#' still caught.
#'
#' Two details are load-bearing. The two estimates are compared **relative
#' to their own magnitude** rather than against a denominator floored at one: near
#' a kink both are tiny yet differ by a factor of two, which a floor of one
#' flattens into apparent agreement. And if no observation survives, all are kept:
#' a systematic disagreement is a real failure and should be reported, not hidden
#' by the guard meant to protect against a local one.
#'
#' Note the placement. This is defined *after* `check_distrib()`, not
#' before it: a roxygen block attaches to whatever object follows it, so a helper
#' slipped in between silently steals the documentation of the function it
#' belongs to. That had already happened once here, leaving the package with a
#' `fd_is_reliable.Rd` and no `check_distrib.Rd`.
#'
#' @param fd_at A function of one argument, the relative step, returning the
#'   finite-difference reference computed with that step.
#' @param ref The reference already computed at `h_rel`, a named list of
#'   component vectors.
#' @param h_rel The relative step `ref` was computed with.
#' @param smooth_all Logical; `TRUE` when every parameter is declared
#'   smooth, in which case no observation is ever dropped and the guard does not
#'   run at all.
#'
#' @return A logical vector as long as the components of `ref`.
#'
#' @seealso [check_distrib()]
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
