#' @include weibull1_distrib.R link_scale.R higher_derivatives.R
NULL

# Fifth derivatives of the weibull1 log-density in the parameters, in closed
# form, from weibull1_deriv5_cpp() (src/d5_weibull1.cpp).
S7::method(distrib_deriv5, Weibull1Distrib) <- function(distrib, y, theta, scale = c("parameter", "link"), ..., threads = 1L) {
  deriv5_scale(distrib, y, theta, weibull1_deriv5_cpp(y, theta[[1]], theta[[2]], threads), match.arg(scale))
}
