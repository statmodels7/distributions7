#' @include wrapper_derivatives.R partition_sums.R link_scale.R folded.R
#' @include fixed.R transformations.R truncated.R zero_inflated.R
#' @include zero_adjusted.R
NULL

#' @title The Wrappers' Fifth Derivatives
#' @name distrib_deriv5.wrappers
#' @aliases distrib_deriv5.ZeroInflatedDistrib
#'   distrib_deriv5.ZeroAdjustedDiscreteDistrib
#'   distrib_deriv5.ZeroAdjustedContinuousDistrib
#'   distrib_deriv5.TruncatedContinuousDistrib
#'   distrib_deriv5.TruncatedDiscreteDistrib distrib_deriv5.FoldedDistrib
#'   distrib_deriv5.TransformedDistrib distrib_deriv5.FixedContinuousDistrib
#'   distrib_deriv5.FixedDiscreteDistrib
#'
#' @description
#' Computes the fifth derivatives of the log-likelihood of [zero_inflated()],
#' [zero_adjusted()], [truncated()], [folded()], [transformation()] and
#' [fixed()] families by the partition sums that give their third and fourth
#' derivatives, one order up, and carries them to the link scale with
#' [deriv5_scale()].
#'
#' @details
#' The zero wrappers, the truncated and the folded families assemble the
#' fifth derivative of \eqn{\log L} from the block ratios \eqn{d^B L / L} by
#' the moment-to-cumulant relation over the 52 partitions of five indices,
#' each ratio a complete Bell polynomial in the parent's derivatives to order
#' five. A transformed family's fifth derivatives are the parent's at the
#' preimage, and a fixed family's are the parent's at the full parameter
#' vector, subset to the free parameters. The result is analytic where the
#' parent's fifth derivative is.
#'
#' @param distrib A family from one of the wrappers above.
#' @param y A numeric vector of observations.
#' @param theta A named list of parameters.
#' @param scale `"parameter"` or `"link"`.
#' @param ... Unused.
#'
#' @return A named list of fifth-derivative components keyed by
#'   [deriv_names()] at order 5, each of length `length(y)`.
#'
#' @seealso [distrib_deriv5()], [log_deriv()], [bell_f_ratio()]
#'
#' @keywords internal
NULL

#' @rdname distrib_deriv5.wrappers
#' @details
#' The partition sums are used only where the parent's fifth derivative is
#' analytic, as [has_exact_deriv5()] reports. Where it is the stencil of
#' [numerical_deriv5()], the products of the Bell polynomials and, for a
#' truncated family, the quadrature of its ratios amplify the stencil's error
#' (3.4e-10 to 1.8e-9 for a folded gaussian, 3.5e-10 to 3.1e-8 for a
#' truncated one, against symbolic derivatives), and the wrapper's fifth
#' derivative is instead one stencil on its own analytic fourth.
#' @param builder One of the order-generic builders of
#'   `R/wrapper_derivatives.R` or [fold_deriv_k()].
#' @usage NULL
wrapper_deriv5_method <- function(builder) {
  kern <- builder(5L)
  function(distrib, y, theta, scale = c("parameter", "link"), ...) {
    scale <- match.arg(scale)
    if (!has_exact_deriv5(distrib@parent_distrib)) {
      return(numerical_deriv5(distrib, y, theta, scale = scale))
    }
    deriv5_scale(distrib, y, theta, kern(distrib, y, theta), scale)
  }
}

#' Is a Family's Fifth Derivative Analytic?
#'
#' @description
#' `TRUE` where the [distrib_deriv5()] method that dispatches is registered on
#' a class other than the base class, whose method is the stencil of
#' [numerical_deriv5()]. A wrapper registered in this file answers for the
#' family it wraps, since its fifth derivative is analytic exactly when that
#' family's is.
#'
#' @param x A distribution object.
#'
#' @return A single logical.
#'
#' @keywords internal
has_exact_deriv5 <- function(x) {
  passing <- list(ZeroInflatedDistrib, ZeroAdjustedDiscreteDistrib,
                  ZeroAdjustedContinuousDistrib, TruncatedContinuousDistrib,
                  TruncatedDiscreteDistrib, FoldedDistrib, TransformedDistrib,
                  FixedContinuousDistrib, FixedDiscreteDistrib)
  if (any(vapply(passing, function(cl) S7::S7_inherits(x, cl), logical(1)))) {
    return(has_exact_deriv5(x@parent_distrib))
  }
  m <- tryCatch(S7::method(distrib_deriv5, S7::S7_class(x)),
                error = function(e) NULL)
  if (is.null(m)) return(FALSE)
  !is_class(attr(m, "signature")[[1]], distrib)
}

S7::method(distrib_deriv5, ZeroInflatedDistrib) <- wrapper_deriv5_method(zi_deriv_k)
S7::method(distrib_deriv5, ZeroAdjustedDiscreteDistrib) <- wrapper_deriv5_method(za_disc_deriv_k)
S7::method(distrib_deriv5, ZeroAdjustedContinuousDistrib) <- wrapper_deriv5_method(za_cont_deriv_k)
S7::method(distrib_deriv5, TruncatedContinuousDistrib) <- wrapper_deriv5_method(trunc_deriv_k)
S7::method(distrib_deriv5, TruncatedDiscreteDistrib) <- wrapper_deriv5_method(trunc_deriv_k)
S7::method(distrib_deriv5, FoldedDistrib) <- wrapper_deriv5_method(fold_deriv_k)

S7::method(distrib_deriv5, TransformedDistrib) <- function(
    distrib, y, theta, scale = c("parameter", "link"), ...) {
  res <- distrib_deriv5(distrib@parent_distrib,
                        distrib@transformer@trans_inv(y), theta)
  deriv5_scale(distrib, y, theta, res, match.arg(scale))
}

for (.fixed_cls in list(FixedContinuousDistrib, FixedDiscreteDistrib)) {
  S7::method(distrib_deriv5, .fixed_cls) <- function(
      distrib, y, theta, scale = c("parameter", "link"), ...) {
    res <- distrib_deriv5(distrib@parent_distrib, y,
                          fixed_full_theta(distrib, theta))
    deriv5_scale(distrib, y, theta, res[deriv_names(distrib@params, 5L)],
                 match.arg(scale))
  }
}
rm(.fixed_cls)
