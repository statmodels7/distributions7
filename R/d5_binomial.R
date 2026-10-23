#' @include binomial_distrib.R link_scale.R higher_derivatives.R
NULL

# Fifth derivative of the Binomial log-mass in mu, in closed form (binomial_deriv5_cpp).
# The size is recycled as in the fourth order; the kernel returns parameter-scale components.
S7::method(distrib_deriv5, BinomialDistrib) <- function(distrib, y, theta, scale = c("parameter", "link"), ..., threads = 1L) {
  a <- binomial_recycle(y, theta[[1]], distrib@size)
  deriv5_scale(distrib, y, theta, binomial_deriv5_cpp(a$y, a$mu, a$size, threads), match.arg(scale))
}
