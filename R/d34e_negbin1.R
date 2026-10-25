#' @include negbin1_distrib.R dexpected_support.R d34expected_analytic.R
NULL

# The third and fourth derivatives of negbin1's expected information,
# compiled in src/d34e_negbin1.cpp.
register_d34expected(NegBin1Distrib, function(d, y, th, k, t)
  (if (k == 3L) negbin1_dexpected3_cpp(y, th[[1]], th[[2]], t)
    else negbin1_dexpected4_cpp(y, th[[1]], th[[2]], t)))
