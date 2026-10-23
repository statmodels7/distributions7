#' @include bernoulli_distrib.R link_scale.R higher_derivatives.R
NULL

# Fifth derivative of the Bernoulli log-mass in mu, in closed form (bernoulli_deriv5_cpp).
# The kernel returns parameter-scale components; deriv5_scale() applies the link.
S7::method(distrib_deriv5, BernoulliDistrib) <- function(distrib, y, theta, scale = c("parameter", "link"), ..., threads = 1L) {
  deriv5_scale(distrib, y, theta, bernoulli_deriv5_cpp(y, theta[[1]], threads), match.arg(scale))
}
