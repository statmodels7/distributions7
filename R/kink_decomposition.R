#' @include generics.R distrib.R
NULL

# THE ORDER OF DIFFERENTIABILITY OF A LOG-DENSITY IN ITS PARAMETERS.
#
# `params_smooth` records a logical per parameter, which is the degenerate
# m in {0, Inf} of a scale that has intermediate rungs on it: a Huber
# likelihood is C^1 and not C^2 in its own cut-off, and the quintic of
# penalties7 is exactly C^3.  What a consumer needs is the integer, because
# each consumer needs a different one -- a score needs one derivative, IWLS
# and the Laplace determinant need two, the exact outer gradient needs three,
# the exact outer Hessian four.
#
# The integer is not declared.  A family declares WHERE its kink is and WHAT
# composition sits on top of it, and the order follows: the order of the
# composition for a parameter the kink's argument moves with, and infinity for
# every other.  One declaration fewer is one declaration fewer that can be
# wrong, which is the reason `deriv_names()` is written in terms of
# `deriv_indices()` rather than as a second list.

#' The Composition That Carries a Kink
#'
#' @description
#' Records the one non-smooth piece of a log-density, as the composition
#' \eqn{c(\theta)\,\phi(v(y, \theta))} in which \eqn{\phi} is not smooth at the
#' origin and everything else is. A family declares one through
#' [kink_decomposition()]; [params_order()] reads the order of
#' differentiability off it.
#'
#' @details
#' The five compositions and the order each one leaves, measured as the first
#' derivative that jumps across \eqn{v = 0}:
#'
#' \tabular{lll}{
#'   `phi`      \tab \eqn{\phi(v)}          \tab order \cr
#'   `"abs"`    \tab \eqn{\lvert v \rvert}  \tab 0 \cr
#'   `"hinge"`  \tab \eqn{(v)_+}            \tab 0 \cr
#'   `"hinge2"` \tab \eqn{(v)_+^2}          \tab 1 \cr
#'   `"hinge3"` \tab \eqn{(v)_+^3}          \tab 2 \cr
#'   `"step"`   \tab \eqn{1\{v > 0\}}       \tab -1
#' }
#'
#' The others reduce to the absolute value, \eqn{(v)_+ = (v + \lvert v
#' \rvert)/2} and \eqn{1\{v>0\} = (1 + \mathrm{sign}(v))/2}, so a smoothing of
#' \eqn{\lvert \cdot \rvert} covers all of them by composition.
#'
#' \eqn{v} is not required to be \eqn{y - \mu}. A Huber likelihood puts its
#' kink at \eqn{v = \lvert y - \mu \rvert - k\sigma}, which moves with the
#' location, the scale and the cut-off together, and the set of non-smooth
#' parameters is read off \eqn{\partial v / \partial \theta_p \ne 0} rather
#' than assumed.
#'
#' @param phi A single string, one of `"abs"`, `"hinge"`, `"hinge2"`,
#'   `"hinge3"` or `"step"`.
#' @param v A function of `(y, theta)` returning the argument of `phi`, one
#'   value per observation.
#' @param dv A function of `(y, theta)` returning a named list with one entry
#'   per parameter, \eqn{\partial v / \partial \theta_p}. An entry that is
#'   identically zero says the kink does not move with that parameter.
#' @param coef A function of `theta` returning the multiplier \eqn{c(\theta)}
#'   standing in front of \eqn{\phi}.
#'
#' @return An object of class `kink_spec` carrying the four fields.
#'
#' @examples
#' # the Laplace: l = -log(2 sigma) - |y - mu| / sigma
#' k <- kink_spec(
#'   phi  = "abs",
#'   v    = function(y, theta) y - theta$mu,
#'   dv   = function(y, theta) list(mu = -1, sigma = 0),
#'   coef = function(theta) -1 / theta$sigma)
#' k
#' kink_order(k@phi)
#'
#' @seealso [kink_decomposition()] for the generic a family answers,
#'   [params_order()] for the order deduced from it, and [check_kink()] for
#'   the validator.
#' @export
kink_spec <- S7::new_class(
  "kink_spec",
  properties = list(
    phi  = S7::class_character,
    v    = S7::class_function,
    dv   = S7::class_function,
    coef = S7::class_function),
  validator = function(self) {
    if (length(self@phi) != 1L || is.na(self@phi) ||
        !self@phi %in% names(.KINK_ORDERS)) {
      return(sprintf("'phi' must be one of %s.",
                     paste(sprintf('"%s"', names(.KINK_ORDERS)), collapse = ", ")))
    }
    if (!identical(names(formals(self@v))[1:2], c("y", "theta")))
      return("'v' must take (y, theta).")
    if (!identical(names(formals(self@dv))[1:2], c("y", "theta")))
      return("'dv' must take (y, theta).")
    if (!identical(names(formals(self@coef))[1L], "theta"))
      return("'coef' must take (theta).")
    NULL
  })

# The order each composition leaves, and the jump of phi' across the origin,
# which is what `check_kink()` holds a declaration to.  Both are measured in
# tests/testthat/test-kink.R rather than asserted here.
.KINK_ORDERS <- c(abs = 0, hinge = 0, hinge2 = 1, hinge3 = 2, step = -1)
.KINK_DPHI   <- c(abs = 2, hinge = 1, hinge2 = 0, hinge3 = 0, step = NA_real_)

#' @export
S7::method(print, kink_spec) <- function(x, ...) {
  cat(sprintf("<kink_spec> c(theta) * phi(v),  phi = %s,  order %s\n",
              x@phi, format(kink_order(x@phi))))
  invisible(x)
}

#' The Order Left by a Composition
#'
#' @description
#' The largest `k` for which a composition has `k` continuous derivatives at
#' the origin: 0 for \eqn{\lvert v \rvert} and \eqn{(v)_+}, whose first
#' derivative jumps; 1 for \eqn{(v)_+^2}; 2 for \eqn{(v)_+^3}; and -1 for
#' \eqn{1\{v>0\}}, which is not continuous at all.
#'
#' @details
#' The order -1 is worse in kind than the others. There a smoothing does not
#' repair a missing derivative but a discontinuity of the log-density itself,
#' and the unsmoothed problem has no well-defined maximum, so it is reported
#' apart rather than treated with the rest.
#'
#' @param phi A single string naming the composition; see [kink_spec()].
#'
#' @return A single number, `-1`, `0`, `1` or `2`.
#'
#' @examples
#' kink_order("abs")
#' kink_order("hinge3")
#'
#' @seealso [kink_spec()], [params_order()].
#' @export
kink_order <- function(phi) {
  if (!is.character(phi) || length(phi) != 1L || is.na(phi) ||
      !phi %in% names(.KINK_ORDERS)) {
    stop(sprintf("'phi' must be one of %s.",
                 paste(sprintf('"%s"', names(.KINK_ORDERS)), collapse = ", ")),
         call. = FALSE)
  }
  unname(.KINK_ORDERS[[phi]])
}

#' The Order of Differentiability of a Family in Each Parameter
#'
#' @description
#' One integer per parameter: the largest `k` such that the log-density has `k`
#' continuous derivatives in that parameter everywhere in `y`. `Inf` says the
#' parameter is smooth, which is the answer for every family that declares no
#' [kink_decomposition()].
#'
#' @details
#' The order is DEDUCED and not declared. Where a family declares a
#' decomposition \eqn{c(\theta)\phi(v(y,\theta))}, a parameter the kink moves
#' with -- one whose \eqn{\partial v/\partial \theta_p} is not identically zero
#' -- inherits the order of \eqn{\phi}, and every other parameter stays `Inf`.
#' Whether \eqn{\partial v/\partial \theta_p} vanishes is measured at probe
#' values rather than assumed, so a Huber cut-off held at zero reports the
#' scale as smooth, which it then is.
#'
#' A family that declares no decomposition but records a parameter as
#' non-smooth through `params_smooth` reports `NA` for it: there is a kink and
#' its order has not been established. A consumer must treat that as unusable
#' rather than as smooth, which is what [check_distrib()] and the fitting layer
#' do. This is the state of every wrapper of a kinked family today.
#'
#' The consumers and the order each one needs:
#'
#' \tabular{ll}{
#'   1 \tab a score continuous in the coefficients \cr
#'   2 \tab IWLS, Newton, `vcov()`, and the determinant of a Laplace criterion \cr
#'   3 \tab the exact gradient of a marginal criterion \cr
#'   4 \tab its exact Hessian \cr
#'   5 \tab the fourth derivative of a filtered predictor
#' }
#'
#' @param distrib An object inheriting from class `"distrib"`.
#' @param ... Passed to [kink_decomposition()].
#'
#' @return A named numeric vector, one entry per element of `distrib@params`.
#'
#' @examples
#' # smooth everywhere
#' params_order(gaussian1_distrib())
#'
#' # the Laplace: the location carries the kink, the scale does not
#' params_order(laplace_distrib())
#'
#' @seealso [kink_decomposition()], [kink_spec()], [param_smoothness()].
#' @export
params_order <- function(distrib, ...) {
  if (!S7::S7_inherits(distrib, distributions7::distrib))
    stop("'distrib' must inherit from class 'distrib'.", call. = FALSE)
  pars <- distrib@params
  out  <- stats::setNames(rep(Inf, length(pars)), pars)
  kd   <- kink_decomposition(distrib, ...)

  if (is.null(kd)) {
    # No decomposition.  A parameter the family itself records as non-smooth
    # has an order, and it has not been established: NA says so.
    sm <- param_smoothness(distrib)
    if (!is.null(sm)) {
      nm <- intersect(names(sm)[!sm], pars)
      out[nm] <- NA_real_
    }
    return(out)
  }

  ord <- kink_order(kd@phi)
  for (p in .kink_touched(distrib, kd)) out[[p]] <- ord
  out
}

#' Which Parameters a Kink Moves With
#'
#' @description
#' The set \eqn{\{p : \partial v/\partial \theta_p \ne 0\}}, evaluated rather
#' than assumed, at probe values inside each parameter's own interval and at a
#' few points of the support.
#'
#' @details
#' The probe is the midpoint of a bounded interval and one unit inside a
#' half-line, which is the rule [reparametrize()] already uses; a component is
#' called non-zero if it is non-zero at any probe, since the set is defined by
#' the derivative not vanishing identically.
#'
#' @param distrib An object inheriting from class `"distrib"`.
#' @param kd A [kink_spec()].
#'
#' @return A character vector of parameter names, possibly empty.
#'
#' @examples
#' d <- laplace_distrib()
#' distributions7:::.kink_touched(d, kink_decomposition(d))
#'
#' @keywords internal
.kink_touched <- function(distrib, kd) {
  th <- .kink_probe_theta(distrib)
  ys <- .kink_probe_y(distrib, th)
  hit <- character(0)
  for (y in ys) {
    dv <- kd@dv(y, th)
    for (p in distrib@params) {
      v <- dv[[p]]
      if (!is.null(v) && any(is.finite(v)) && any(v[is.finite(v)] != 0))
        hit <- c(hit, p)
    }
  }
  intersect(distrib@params, unique(hit))
}

#' A Trial Parameter Value Inside Every Interval
#'
#' @description
#' The midpoint of a bounded interval, one unit inside a half-line, and zero on
#' the whole line. The same rule [reparametrize()] probes a map with.
#'
#' @param distrib An object inheriting from class `"distrib"`.
#'
#' @return A named list, one entry per parameter.
#'
#' @examples
#' distributions7:::.kink_probe_theta(laplace_distrib())
#'
#' @keywords internal
.kink_probe_theta <- function(distrib) {
  b <- distrib@params_bounds
  stats::setNames(lapply(distrib@params, function(k) {
    bk <- b[[k]]
    if (is.null(bk)) return(0)
    if (is.finite(bk[1L]) && is.finite(bk[2L])) return(mean(bk))
    if (is.finite(bk[1L])) return(bk[1L] + 1)
    if (is.finite(bk[2L])) return(bk[2L] - 1)
    0
  }), distrib@params)
}

#' A Few Points of the Support
#'
#' @description
#' Three quantiles of the family at the probe parameters, falling back on the
#' interior of its bounds where the quantile function refuses.
#'
#' @param distrib An object inheriting from class `"distrib"`.
#' @param theta A named list of parameter values.
#'
#' @return A numeric vector.
#'
#' @examples
#' d <- laplace_distrib()
#' distributions7:::.kink_probe_y(d, distributions7:::.kink_probe_theta(d))
#'
#' @keywords internal
.kink_probe_y <- function(distrib, theta) {
  y <- tryCatch(distrib_quantile(distrib, c(0.25, 0.5, 0.75), theta),
                error = function(e) NULL)
  if (!is.null(y) && all(is.finite(y))) return(y)
  b <- distrib@bounds
  lo <- if (is.finite(b[1L])) b[1L] + 1 else -1
  hi <- if (is.finite(b[2L])) b[2L] - 1 else  1
  unique(c(lo, (lo + hi) / 2, hi))
}

#' @name kink_decomposition.distrib
#' @title The Base Method: No Kink
#'
#' @description
#' Returns `NULL`: the family is smooth in every parameter. Every family
#' inherits this except the three that declare an absolute value of the
#' residual in their log-density.
#'
#' @param distrib An object inheriting from class `"distrib"`.
#'
#' @return `NULL`.
#'
#' @examples
#' kink_decomposition(gaussian1_distrib())
#'
#' @seealso [kink_decomposition()] for the generic.
#' @keywords internal
S7::method(kink_decomposition, distrib) <- function(distrib) NULL

#' Check a Declared Kink Against the Family
#'
#' @description
#' Holds a [kink_decomposition()] to the family that declares it. The
#' declaration says where the kink is, what sits on top of it and what stands
#' in front; the family's own density and score are written independently of
#' it, so the two can be compared.
#'
#' @details
#' Three checks, at probe parameters inside every interval:
#'
#' \describe{
#'   \item{`dv`}{each \eqn{\partial v/\partial \theta_p} against a central
#'     difference of \eqn{v} itself. This isolates `dv`, which nothing else
#'     reads.}
#'   \item{`jump`}{the jump of the score across the kink against the one the
#'     declaration predicts. Writing \eqn{\Delta\phi'} for the jump of
#'     \eqn{\phi'} at the origin -- 2 for \eqn{\lvert v \rvert} and 1 for
#'     \eqn{(v)_+} -- the score of \eqn{c(\theta)\phi(v)} jumps by
#'     \deqn{c(\theta)\, \Delta\phi' \, \partial v/\partial \theta_p}
#'     as \eqn{y} crosses the surface \eqn{v = 0}. This is the check that
#'     matters: it reads `phi`, `v`, `dv` and `coef` at once, and a wrong `v`
#'     fails it.}
#'   \item{`smooth`}{the score of every parameter the declaration calls smooth
#'     really is continuous across the same surface. Without this the check
#'     would pass a declaration that named too few parameters.}
#' }
#'
#' A composition whose \eqn{\Delta\phi'} is zero -- \eqn{(v)_+^2} and
#' \eqn{(v)_+^3}, where the second and third derivative jump instead -- is
#' reported as not checked rather than checked against zero, which any
#' declaration would satisfy.
#'
#' @param distrib An object inheriting from class `"distrib"`.
#' @param theta Optional named list of parameter values. Defaults to a probe
#'   inside every interval.
#' @param tol Relative tolerance for the two numerical comparisons.
#' @param verbose Print the table.
#'
#' @return Invisibly, a named logical vector with one entry per check, `NA`
#'   where a check does not apply.
#'
#' @examples
#' check_kink(laplace_distrib())
#' check_kink(laplace2_distrib())
#'
#' @seealso [kink_decomposition()], [kink_spec()], [params_order()],
#'   [check_distrib()] for the family's own battery.
#' @export
check_kink <- function(distrib, theta = NULL, tol = 1e-6, verbose = TRUE) {
  kd <- kink_decomposition(distrib)
  if (is.null(kd)) {
    stop(sprintf("'%s' declares no kink_decomposition(); there is nothing to check.",
                 distrib@distrib_name), call. = FALSE)
  }
  if (is.null(theta)) theta <- .kink_probe_theta(distrib)
  pars    <- distrib@params
  touched <- .kink_touched(distrib, kd)
  ord     <- kink_order(kd@phi)
  dphi    <- unname(.KINK_DPHI[[kd@phi]])

  out <- c(dv = NA, jump = NA, smooth = NA)

  # (1) dv against a central difference of v.
  h <- 1e-5
  gaps <- vapply(pars, function(p) {
    tp <- tm <- theta
    step <- h * max(abs(theta[[p]]), 1)
    tp[[p]] <- theta[[p]] + step; tm[[p]] <- theta[[p]] - step
    y  <- .kink_probe_y(distrib, theta)
    num <- (kd@v(y, tp) - kd@v(y, tm)) / (2 * step)
    ana <- rep_len(kd@dv(y, theta)[[p]], length(y))
    max(abs(num - ana)) / max(1, max(abs(ana)))
  }, numeric(1))
  out[["dv"]] <- all(gaps < tol)

  # (2) and (3): the jump of the score across v = 0.
  ystar <- .kink_locate(distrib, kd, theta)
  if (is.na(ystar) || is.na(dphi) || dphi == 0) {
    if (verbose) .kink_report(out, gaps, NULL, NULL, kd, ord, distrib)
    return(invisible(out))
  }
  e  <- 1e-6 * max(abs(ystar), 1)
  gp <- distrib_gradient(distrib, ystar + e, theta)
  gm <- distrib_gradient(distrib, ystar - e, theta)
  cc <- kd@coef(theta)
  dv <- kd@dv(ystar, theta)

  jumps <- vapply(pars, function(p) unname(gp[[p]] - gm[[p]]), numeric(1))
  pred  <- vapply(pars, function(p) {
    if (!p %in% touched) return(0)
    unname(cc * dphi * rep_len(dv[[p]], 1L))
  }, numeric(1))
  sc <- max(1, max(abs(pred)))
  out[["jump"]]   <- all(abs(jumps[touched] - pred[touched]) < tol * sc)
  rest <- setdiff(pars, touched)
  out[["smooth"]] <- length(rest) == 0L || all(abs(jumps[rest]) < tol * sc)

  if (verbose) .kink_report(out, gaps, jumps, pred, kd, ord, distrib)
  invisible(out)
}

#' Where the Kink Sits on the Response Scale
#'
#' @description
#' Solves \eqn{v(y, \theta) = 0} for `y`, which is where the log-density is not
#' smooth. Returns `NA` when no sign change is bracketed inside the support.
#'
#' @param distrib An object inheriting from class `"distrib"`.
#' @param kd A [kink_spec()].
#' @param theta A named list of parameter values.
#'
#' @return A single number, or `NA_real_`.
#'
#' @examples
#' d <- laplace_distrib()
#' distributions7:::.kink_locate(d, kink_decomposition(d), list(mu = 0.5, sigma = 2))
#'
#' @keywords internal
.kink_locate <- function(distrib, kd, theta) {
  b  <- distrib@bounds
  lo <- if (is.finite(b[1L])) b[1L] else -1e3
  hi <- if (is.finite(b[2L])) b[2L] else  1e3
  f  <- function(y) kd@v(y, theta)
  fl <- tryCatch(f(lo + 1e-8), error = function(e) NA_real_)
  fh <- tryCatch(f(hi - 1e-8), error = function(e) NA_real_)
  if (!is.finite(fl) || !is.finite(fh) || fl * fh > 0) return(NA_real_)
  tryCatch(stats::uniroot(f, c(lo + 1e-8, hi - 1e-8), tol = 1e-12)$root,
           error = function(e) NA_real_)
}

#' Print the Table of check_kink
#'
#' @description Prints one line per check, with the worst gap behind it.
#'
#' @param out The logical vector of results.
#' @param gaps Per-parameter gaps of the `dv` check.
#' @param jumps Measured jumps of the score, or `NULL`.
#' @param pred Predicted jumps, or `NULL`.
#' @param kd A [kink_spec()].
#' @param ord The order the composition leaves.
#' @param distrib The family.
#'
#' @return `NULL`, invisibly; called for the printing.
#'
#' @examples
#' check_kink(laplace_distrib(), verbose = TRUE)
#'
#' @keywords internal
.kink_report <- function(out, gaps, jumps, pred, kd, ord, distrib) {
  cat(sprintf("check_kink: %s,  phi = %s,  order %s\n",
              distrib@distrib_name, kd@phi, format(ord)))
  mark <- function(v) if (is.na(v)) "[not checked]" else if (v) "[PASSED]" else "[FAILED]"
  cat(sprintf("  dv against a difference of v      %s  worst %.2e\n",
              mark(out[["dv"]]), max(gaps)))
  if (is.null(jumps)) {
    cat(sprintf("  the jump of the score             %s  phi' does not jump for \"%s\"\n",
                mark(out[["jump"]]), kd@phi))
    cat(sprintf("  the smooth parameters             %s\n", mark(out[["smooth"]])))
  } else {
    cat(sprintf("  the jump of the score             %s  worst %.2e\n",
                mark(out[["jump"]]), max(abs(jumps - pred))))
    cat(sprintf("  the smooth parameters             %s\n", mark(out[["smooth"]])))
    for (p in names(jumps))
      cat(sprintf("      %-10s measured %12.6g   declared %12.6g\n",
                  p, jumps[[p]], pred[[p]]))
  }
  invisible(NULL)
}
