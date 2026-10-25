#' @include poisson_distrib.R d34expected_analytic.R
NULL

# The third and fourth derivatives of poisson's expected information,
# compiled in src/d34e_poisson.cpp.
register_d34expected(PoissonDistrib, function(d, y, th, k, t)
  (if (k == 3L) poisson_dexpected3_cpp(y, th[[1]], t)
    else poisson_dexpected4_cpp(y, th[[1]], t)))
