#' @include vonmises1_distrib.R link_scale.R higher_derivatives.R
NULL

# The fifth derivatives of the von Mises log-density in (mu, kappa), from
# vonmises1_deriv5_cpp(); the kernel returns parameter-scale components and
# deriv5_scale() applies the link scale.
S7::method(distrib_deriv5, VonMises1Distrib) <- function(distrib, y, theta,
                                                         scale = c("parameter", "link"),
                                                         ..., threads = 1L) {
  deriv5_scale(distrib, y, theta, vonmises1_deriv5_cpp(y, theta[[1]], theta[[2]], threads),
               match.arg(scale))
}
