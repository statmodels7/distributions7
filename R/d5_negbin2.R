#' @include negbin2_distrib.R link_scale.R higher_derivatives.R
NULL

# Fifth derivatives of the NB2 log-mass in (mu, theta), sums over the support for large theta.
# The kernel is negbin2_deriv5_cpp(), with the link scale applied by deriv5_scale().
S7::method(distrib_deriv5, NegBin2Distrib) <- function(distrib, y, theta, scale = c("parameter", "link"), threads = 1L, ...) {
  deriv5_scale(distrib, y, theta, negbin2_deriv5_cpp(y, theta[[1]], theta[[2]], threads), match.arg(scale))
}
