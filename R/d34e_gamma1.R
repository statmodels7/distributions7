#' @include gamma1_distrib.R d34expected_analytic.R
NULL

# The third and fourth derivatives of gamma1's expected information,
# compiled in src/d34e_gamma1.cpp.
register_d34expected(Gamma1Distrib, function(d, y, th, k, t)
  (if (k == 3L) gamma1_dexpected3_cpp(y, th[[1]], th[[2]], t)
    else gamma1_dexpected4_cpp(y, th[[1]], th[[2]], t)))
