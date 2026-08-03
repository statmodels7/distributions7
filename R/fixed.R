#' @include distrib.R generics.R utility_functions.R moments.R cross_derivatives.R
NULL

#' @title S7 Class for Distributions With Fixed Parameters (Continuous)
#' @name FixedContinuousDistrib
#'
#' @description
#' A subclass of \code{continuous_distrib} representing a continuous
#' distribution in which some parameters of the wrapped distribution are held
#' at known values. Constructed by \code{\link{fixed}}.
#'
#' @details
#' The free parameters are the parent's minus the fixed ones, in the parent's
#' order. Every method splices the fixed values back into \code{theta} at their
#' positions and delegates to the parent, so the parent's closed forms are used
#' whenever they exist; a derivative method then keeps only the components in
#' which every index is a free parameter. No method of this class computes
#' anything of its own.
#'
#' @inheritParams distrib
#' @param parent_distrib The wrapped \code{continuous_distrib} object.
#' @param fixed_params A named list of the fixed parameter values.
#' @seealso \code{\link{fixed}}
FixedContinuousDistrib <- S7::new_class("FixedContinuousDistrib",
  parent = continuous_distrib,
  properties = list(
    parent_distrib = distrib,
    fixed_params = S7::class_list
  )
)

#' @title S7 Class for Distributions With Fixed Parameters (Discrete)
#' @name FixedDiscreteDistrib
#'
#' @description
#' A subclass of \code{discrete_distrib} representing a discrete distribution
#' in which some parameters of the wrapped distribution are held at known
#' values. Constructed by \code{\link{fixed}}.
#'
#' @details
#' Identical in behaviour to \code{\link{FixedContinuousDistrib}}: every method
#' splices the fixed values into \code{theta} and delegates to the parent, and
#' the derivative methods keep only the components among the free parameters.
#'
#' @inheritParams distrib
#' @param parent_distrib The wrapped \code{discrete_distrib} object.
#' @param fixed_params A named list of the fixed parameter values.
#' @seealso \code{\link{fixed}}
FixedDiscreteDistrib <- S7::new_class("FixedDiscreteDistrib",
  parent = discrete_distrib,
  properties = list(
    parent_distrib = distrib,
    fixed_params = S7::class_list
  )
)

#' Is This a Fixed-Parameter Wrapper?
#'
#' @description
#' \code{TRUE} for a distribution produced by \code{\link{fixed}}, in either of
#' its two forms.
#'
#' @param distrib An object inheriting from class \code{"distrib"}.
#'
#' @return A single logical.
#'
#' @seealso \code{\link{fixed}}
#' @keywords internal
is_fixed <- function(distrib) {
  S7::S7_inherits(distrib, FixedContinuousDistrib) ||
    S7::S7_inherits(distrib, FixedDiscreteDistrib)
}

#' Splice the Fixed Values Back Into a Full Parameter List
#'
#' @description
#' Combines the wrapper's free \code{theta} with its fixed values into the full
#' parameter list of the parent, in the parent's order.
#'
#' @details
#' \code{theta} is aligned against the wrapper first, so the function is safe
#' to call both from generic-dispatched methods, whose \code{theta} is already
#' aligned, and from delegating methods such as \code{mean()}, whose
#' \code{theta} arrives as the caller wrote it. Free values may be vectors --
#' the wrapper is as vectorised in \code{theta} as its parent -- while the
#' fixed values are scalars by construction.
#'
#' @param distrib A fixed-parameter wrapper object.
#' @param theta A named list or vector of the free parameters.
#'
#' @return A named list covering every parameter of the parent.
#'
#' @seealso \code{\link{fixed}}
#' @keywords internal
fixed_full_theta <- function(distrib, theta) {
  theta <- align_theta(distrib, theta)
  free <- distrib@params
  theta <- theta[seq_len(distrib@n_params)]
  names(theta) <- free

  parent <- distrib@parent_distrib
  out <- vector("list", parent@n_params)
  names(out) <- parent@params
  out[free] <- theta
  fp <- distrib@fixed_params
  out[names(fp)] <- fp
  out
}

# ---------------------------------------------------------------------------
# Method registration.
#
# Every method has the same one-line behaviour -- splice and delegate -- so
# they are registered in a loop over the two classes. The functions do not
# read the loop variable, so no closure capture is involved. The derivative
# methods subset the parent's result by the names generated from the free
# parameter set: a combination of free parameters produces the same name
# string under the parent's enumeration as under the wrapper's, because the
# free set preserves the parent's order, so subsetting by name cannot pair a
# component with the wrong indices -- the mistake that re-parsing names by
# splitting on the underscore commits for a parameter whose own name contains
# one.
# ---------------------------------------------------------------------------

for (.fixed_cls in list(FixedContinuousDistrib, FixedDiscreteDistrib)) {
  S7::method(distrib_pdf, .fixed_cls) <- function(distrib, y, theta, log = FALSE) {
    distrib_pdf(distrib@parent_distrib, y, fixed_full_theta(distrib, theta),
      log = log
    )
  }

  S7::method(distrib_cdf, .fixed_cls) <- function(distrib, q, theta, ...) {
    distrib_cdf(distrib@parent_distrib, q, fixed_full_theta(distrib, theta), ...)
  }

  S7::method(distrib_quantile, .fixed_cls) <- function(distrib, p, theta, ...) {
    distrib_quantile(distrib@parent_distrib, p, fixed_full_theta(distrib, theta), ...)
  }

  S7::method(distrib_rng, .fixed_cls) <- function(distrib, n, theta, ...) {
    distrib_rng(distrib@parent_distrib, n, fixed_full_theta(distrib, theta), ...)
  }

  S7::method(distrib_atoms, .fixed_cls) <- function(distrib, theta, ...) {
    distrib_atoms(distrib@parent_distrib, fixed_full_theta(distrib, theta), ...)
  }

  S7::method(distrib_grad_y, .fixed_cls) <- function(distrib, y, theta, ...) {
    distrib_grad_y(distrib@parent_distrib, y, fixed_full_theta(distrib, theta), ...)
  }

  S7::method(distrib_hess_y, .fixed_cls) <- function(distrib, y, theta, ...) {
    distrib_hess_y(distrib@parent_distrib, y, fixed_full_theta(distrib, theta), ...)
  }

  S7::method(distrib_cross_y, .fixed_cls) <- function(distrib, y, theta,
                                                      scale = c("parameter", "link"),
                                                      ...) {
    res <- distrib_cross_y(distrib@parent_distrib, y,
      fixed_full_theta(distrib, theta),
      scale = "parameter", ...
    )
    res[distrib@params]
  }

  S7::method(distrib_gradient, .fixed_cls) <- function(distrib, y, theta,
                                                       scale = c("parameter", "link"),
                                                       ...) {
    res <- distrib_gradient(distrib@parent_distrib, y,
      fixed_full_theta(distrib, theta),
      scale = "parameter", ...
    )
    res[distrib@params]
  }

  S7::method(distrib_hessian, .fixed_cls) <- function(distrib, y, theta,
                                                      scale = c("parameter", "link"),
                                                      ...) {
    res <- distrib_hessian(distrib@parent_distrib, y,
      fixed_full_theta(distrib, theta),
      scale = "parameter", ...
    )
    res[hess_names(distrib@params)]
  }

  S7::method(distrib_expected_hessian, .fixed_cls) <- function(distrib, y, theta,
                                                               scale = c("parameter", "link"),
                                                               approx = c("bartlett", "integrate", "mc", "opg"),
                                                               nsim = 10000, ...) {
    res <- distrib_expected_hessian(distrib@parent_distrib, y,
      fixed_full_theta(distrib, theta),
      scale = "parameter", approx = approx, nsim = nsim, ...
    )
    res[hess_names(distrib@params)]
  }

  S7::method(distrib_deriv3, .fixed_cls) <- function(distrib, y, theta,
                                                     expected = FALSE,
                                                     scale = c("parameter", "link"),
                                                     approx = c("integrate", "bartlett", "mc", "opg"),
                                                     nsim = 10000, ...) {
    res <- distrib_deriv3(distrib@parent_distrib, y,
      fixed_full_theta(distrib, theta),
      expected = expected, scale = "parameter", approx = approx, nsim = nsim, ...
    )
    res[deriv_names(distrib@params, 3L)]
  }

  S7::method(distrib_deriv4, .fixed_cls) <- function(distrib, y, theta,
                                                     expected = FALSE,
                                                     scale = c("parameter", "link"),
                                                     approx = c("integrate", "bartlett", "mc", "opg"),
                                                     nsim = 10000, ...) {
    res <- distrib_deriv4(distrib@parent_distrib, y,
      fixed_full_theta(distrib, theta),
      expected = expected, scale = "parameter", approx = approx, nsim = nsim, ...
    )
    res[deriv_names(distrib@params, 4L)]
  }

  S7::method(distrib_grad_cdf, .fixed_cls) <- function(distrib, q, theta,
                                                       lower.tail = TRUE, log = TRUE,
                                                       ...) {
    res <- distrib_grad_cdf(distrib@parent_distrib, q,
      fixed_full_theta(distrib, theta),
      lower.tail = lower.tail, log = log, ...
    )
    res[distrib@params]
  }

  S7::method(distrib_hess_cdf, .fixed_cls) <- function(distrib, q, theta,
                                                       lower.tail = TRUE, log = TRUE,
                                                       ...) {
    res <- distrib_hess_cdf(distrib@parent_distrib, q,
      fixed_full_theta(distrib, theta),
      lower.tail = lower.tail, log = log, ...
    )
    res[hess_names(distrib@params)]
  }

  # The moments delegate so that a parent with a closed form keeps it; the
  # law is the parent's law at the full parameter vector, so nothing changes
  # but where theta comes from.
  S7::method(mean, .fixed_cls) <- function(x, theta, ...) {
    mean(x@parent_distrib, fixed_full_theta(x, theta), ...)
  }

  S7::method(variance, .fixed_cls) <- function(x, theta, ...) {
    variance(x@parent_distrib, fixed_full_theta(x, theta), ...)
  }

  S7::method(std_dev, .fixed_cls) <- function(x, theta, ...) {
    std_dev(x@parent_distrib, fixed_full_theta(x, theta), ...)
  }

  S7::method(skewness, .fixed_cls) <- function(x, theta, ...) {
    skewness(x@parent_distrib, fixed_full_theta(x, theta), ...)
  }

  S7::method(kurtosis, .fixed_cls) <- function(x, theta, ...) {
    kurtosis(x@parent_distrib, fixed_full_theta(x, theta), ...)
  }

  S7::method(print, .fixed_cls) <- function(x, ...) {
    # The base print iterates over the free parameters, which works down to
    # one of them; with none it computes the width of an empty name set, so
    # that case prints its own header instead. print() is a base generic and a
    # method registered on it is an S3 method, so the parent class's print is
    # reached with NextMethod(); super() only works inside S7 generics, and
    # S7::method() only retrieves from them.
    if (x@n_params > 0L) {
      NextMethod()
    } else {
      d_name <- gsub("(^|[[:space:]])([[:alpha:]])", "\\1\\U\\2",
        x@distrib_name,
        perl = TRUE
      )
      d_type <- if (S7::S7_inherits(x, continuous_distrib)) "Continuous" else "Discrete"
      cat(sprintf("Distribution: %s\n", d_name))
      cat(sprintf("Type:         %s\n", d_type))
      cat(sprintf("Dimensions:   %s\n", x@dimension))
      cat("\nParameters:   none free\n")
    }
    fp <- x@fixed_params
    cat("\nFixed:\n")
    for (nm in names(fp)) {
      cat(sprintf("  %s = %s\n", nm, format(fp[[nm]])))
    }
    invisible(x)
  }
}
rm(.fixed_cls)

#' Fix Parameters of a Distribution at Known Values
#'
#' @description
#' Returns the distribution obtained by holding some parameters of
#' \code{distrib} at known values, leaving only the others to be supplied and
#' estimated.
#'
#' @details
#' The result is the same law with a smaller parameter set: \code{theta}
#' carries only the free parameters, every generic answers as the parent does
#' at the full vector, and the derivative components are the parent's
#' restricted to the free indices -- a subvector of the score, a submatrix of
#' the Hessian, sub-arrays at orders three and four. Nothing is recomputed and
#' no normalising constant changes, so the parent's closed forms are used
#' throughout, and \code{\link{fit_distrib}} estimates the free parameters
#' with standard errors and intervals for them alone.
#'
#' Fixed values are single numbers, strictly inside the open domain of their
#' parameter. Fixing a parameter of a distribution that is already a
#' fixed-parameter wrapper collapses the two into one wrapper around the
#' original parent. Fixing every parameter is allowed and gives a fully known
#' distribution with an empty parameter set.
#'
#' The per-parameter smoothness declaration travels with the free parameters,
#' so fixing the location of a Laplace distribution leaves a distribution
#' whose remaining parameter is smooth.
#'
#' @param distrib The distribution whose parameters are to be fixed.
#' @param ... The fixed values, named after the parameters they fix, as in
#'   \code{fixed(gaussian_distrib(), mu = 0)}.
#'
#' @return An object of class \code{FixedContinuousDistrib} or
#'   \code{FixedDiscreteDistrib}, matching the parent.
#'
#' @seealso \code{\link{zero_inflated}}, \code{\link{truncated}},
#'   \code{\link{transformation}}
#'
#' @examples
#' # a gaussian with known mean: only sigma remains
#' d <- fixed(gaussian_distrib(), mu = 0)
#' d@params
#'
#' theta <- list(sigma = 2)
#' distrib_pdf(d, c(-1, 0, 1), theta)
#' distrib_gradient(d, c(-1, 0, 1), theta)
#'
#' # the score is the corresponding component of the parent's
#' full <- distrib_gradient(gaussian_distrib(), c(-1, 0, 1), list(mu = 0, sigma = 2))
#' all.equal(distrib_gradient(d, c(-1, 0, 1), theta)$sigma, full$sigma)
#'
#' # fixing everything gives a fully known distribution
#' d0 <- fixed(gaussian_distrib(), mu = 0, sigma = 1)
#' distrib_pdf(d0, 0, list())
#'
#' @export
fixed <- function(distrib, ...) {
  if (!S7::S7_inherits(distrib, continuous_distrib) &&
    !S7::S7_inherits(distrib, discrete_distrib)) {
    stop("Input must inherit from 'discrete_distrib' or 'continuous_distrib'.",
      call. = FALSE
    )
  }

  fix <- list(...)
  if (length(fix) == 0L) {
    stop(paste0(
      "fixed() needs at least one named value, as in fixed(d, mu = 0). With\n",
      "  none, the result would be the parent distribution unchanged; returning\n",
      "  it silently would hide a missing argument rather than report it."
    ), call. = FALSE)
  }
  nms <- names(fix)
  if (is.null(nms) || any(!nzchar(nms))) {
    stop("Every fixed value must be named after the parameter it fixes.",
      call. = FALSE
    )
  }
  if (anyDuplicated(nms)) {
    stop(sprintf(
      "Parameter '%s' is fixed more than once in the same call.",
      nms[duplicated(nms)][1L]
    ), call. = FALSE)
  }

  # Collapse a fixed() of a fixed(): one wrapper around the original parent,
  # carrying both sets of values. A parameter already fixed is caught below
  # by the membership check, since it is no longer among the free parameters.
  inherited <- list()
  if (is_fixed(distrib)) {
    unknown <- setdiff(nms, distrib@params)
    if (length(unknown)) {
      stop(sprintf(
        "'%s' is not a free parameter of '%s'. Free parameters: %s.",
        unknown[1L], distrib@distrib_name,
        paste(distrib@params, collapse = ", ")
      ), call. = FALSE)
    }
    inherited <- distrib@fixed_params
    distrib <- distrib@parent_distrib
  }

  unknown <- setdiff(nms, distrib@params)
  if (length(unknown)) {
    stop(sprintf(
      "'%s' is not a parameter of '%s'. Parameters: %s.",
      unknown[1L], distrib@distrib_name,
      paste(distrib@params, collapse = ", ")
    ), call. = FALSE)
  }

  for (nm in nms) {
    v <- fix[[nm]]
    if (!is.numeric(v) || length(v) != 1L || !is.finite(v)) {
      stop(sprintf(
        "The value fixing '%s' must be a single finite number.", nm
      ), call. = FALSE)
    }
    b <- distrib@params_bounds[[nm]]
    # The domains are open intervals, the convention align_theta() enforces
    # for a free parameter; a fixed one obeys the same rule.
    if (v <= b[1] || v >= b[2]) {
      stop(sprintf(
        "The value fixing '%s' (%s) is outside its open domain (%s, %s).",
        nm, format(v), format(b[1]), format(b[2])
      ), call. = FALSE)
    }
    fix[[nm]] <- as.numeric(v)
  }

  all_fixed <- c(inherited, fix)
  # Keep the parent's ordering in the record, so the printed name is stable
  # whatever order the calls arrived in.
  all_fixed <- all_fixed[intersect(distrib@params, names(all_fixed))]
  free <- setdiff(distrib@params, names(all_fixed))

  smooth <- param_smoothness(distrib)[free]
  if (length(free) == 0L) smooth <- logical(0)

  # No spaces in the label: the print method capitalises the first letter
  # after every space in the distribution's name, and a parameter name must
  # not come out as "Sigma". Same convention as truncated's "[lower=0]".
  label <- paste(
    vapply(
      names(all_fixed),
      function(nm) sprintf("%s=%s", nm, format(all_fixed[[nm]])),
      character(1)
    ),
    collapse = ","
  )

  common <- list(
    parent_distrib = distrib,
    fixed_params = all_fixed,
    distrib_name = sprintf("fixed %s [%s]", distrib@distrib_name, label),
    dimension = distrib@dimension,
    bounds = distrib@bounds,
    params = free,
    params_interpretation = distrib@params_interpretation[free],
    n_params = length(free),
    params_bounds = distrib@params_bounds[free],
    link_params = distrib@link_params[free],
    params_smooth = smooth
  )

  if (S7::S7_inherits(distrib, discrete_distrib)) {
    do.call(FixedDiscreteDistrib, common)
  } else {
    do.call(FixedContinuousDistrib, common)
  }
}
