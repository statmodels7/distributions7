#' @include betabinom1_distrib.R link_scale.R higher_derivatives.R
NULL

# Fifth derivatives of the BB1 log-mass in (mu, sigma), sums over the support below sigma 0.1 and polygamma pairs above.
# The kernel is betabinom1_deriv5_cpp(), with the link scale applied by deriv5_scale().
S7::method(distrib_deriv5, BetaBinom1Distrib) <- function(distrib, y, theta, scale = c("parameter", "link"), ...) {
  deriv5_scale(distrib, y, theta, betabinom1_deriv5_cpp(y, theta[[1]], theta[[2]], distrib@size), match.arg(scale))
}
