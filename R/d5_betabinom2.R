#' @include betabinom2_distrib.R link_scale.R higher_derivatives.R
NULL

# Fifth derivatives of the BB2 log-mass in (alpha, beta), sums over the support above alpha + beta 10 and polygamma pairs below.
# The kernel is betabinom2_deriv5_cpp(), with the link scale applied by deriv5_scale().
S7::method(distrib_deriv5, BetaBinom2Distrib) <- function(distrib, y, theta, scale = c("parameter", "link"), ...) {
  deriv5_scale(distrib, y, theta, betabinom2_deriv5_cpp(y, theta[[1]], theta[[2]], distrib@size), match.arg(scale))
}
