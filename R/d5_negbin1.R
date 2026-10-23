#' @include negbin1_distrib.R link_scale.R higher_derivatives.R
NULL

# Fifth derivatives of the NB1 log-mass in (mu, theta), size route above theta 1 and mu + i theta route below.
# The kernel is negbin1_deriv5_cpp(), with the link scale applied by deriv5_scale().
S7::method(distrib_deriv5, NegBin1Distrib) <- function(distrib, y, theta, scale = c("parameter", "link"), ...) {
  deriv5_scale(distrib, y, theta, negbin1_deriv5_cpp(y, theta[[1]], theta[[2]]), match.arg(scale))
}
