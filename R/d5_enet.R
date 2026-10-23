#' @include enet_distrib.R link_scale.R higher_derivatives.R
NULL

# The fifth derivatives of the log-density in mu, lambda and alpha, from enet_deriv5_cpp();
# the kernel returns parameter-scale components and deriv5_scale() applies the
# link scale.
S7::method(distrib_deriv5, EnetDistrib) <- function(distrib, y, theta,
                                                         scale = c("parameter", "link"),
                                                         ...) {
  deriv5_scale(distrib, y, theta, enet_deriv5_cpp(y, theta[[1]], theta[[2]], theta[[3]]),
               match.arg(scale))
}
