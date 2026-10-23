#' @include beta2_distrib.R link_scale.R higher_derivatives.R
NULL

# Fifth derivatives of the beta2 log-density in the parameters, in closed
# form, from beta2_deriv5_cpp() (src/d5_beta2.cpp).
S7::method(distrib_deriv5, Beta2Distrib) <- function(distrib, y, theta, scale = c("parameter", "link"), ...) {
  deriv5_scale(distrib, y, theta, beta2_deriv5_cpp(y, theta[[1]], theta[[2]]), match.arg(scale))
}
