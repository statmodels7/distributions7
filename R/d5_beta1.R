#' @include beta1_distrib.R link_scale.R higher_derivatives.R
NULL

# Fifth derivatives of the beta1 log-density in the parameters, in closed
# form, from beta1_deriv5_cpp() (src/d5_beta1.cpp).
S7::method(distrib_deriv5, Beta1Distrib) <- function(distrib, y, theta, scale = c("parameter", "link"), ..., threads = 1L) {
  deriv5_scale(distrib, y, theta, beta1_deriv5_cpp(y, theta[[1]], theta[[2]], threads), match.arg(scale))
}
