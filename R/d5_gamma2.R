#' @include gamma2_distrib.R link_scale.R higher_derivatives.R
NULL

# Fifth derivative of the gamma log-density in (mu, sigma2), in closed form (gamma2_deriv5_cpp).
# The kernel returns parameter-scale components; deriv5_scale() applies the link.
S7::method(distrib_deriv5, Gamma2Distrib) <- function(distrib, y, theta, scale = c("parameter", "link"), ..., threads = 1L) {
  deriv5_scale(distrib, y, theta, gamma2_deriv5_cpp(y, theta[[1]], theta[[2]], threads), match.arg(scale))
}
