#' @include lognormal1_distrib.R link_scale.R higher_derivatives.R
NULL

# Fifth derivatives of the lognormal1 log-density in the parameters, in closed
# form, from lognormal1_deriv5_cpp() (src/d5_lognormal1.cpp).
S7::method(distrib_deriv5, Lognormal1Distrib) <- function(distrib, y, theta, scale = c("parameter", "link"), ..., threads = 1L) {
  deriv5_scale(distrib, y, theta, lognormal1_deriv5_cpp(y, theta[[1]], theta[[2]], threads), match.arg(scale))
}
