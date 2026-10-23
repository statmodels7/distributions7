#' @include pig2_distrib.R link_scale.R higher_derivatives.R
NULL

# Fifth derivatives of the pig2 log-mass in closed form, from pig2_deriv5_cpp()
# (src/d5_pig2.cpp), which returns the parameter scale; the link scale is
# applied by deriv5_scale().
S7::method(distrib_deriv5, Pig2Distrib) <- function(distrib, y, theta,
    scale = c("parameter", "link"), ..., threads = 1L) {
  deriv5_scale(distrib, y, theta, pig2_deriv5_cpp(y, theta[[1]], theta[[2]], threads),
               match.arg(scale))
}
