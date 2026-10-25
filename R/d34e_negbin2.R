#' @include negbin2_distrib.R d34expected_analytic.R
NULL

# The third and fourth derivatives of negbin2's expected information,
# compiled in src/d34e_negbin2.cpp.
register_d34expected(NegBin2Distrib, function(d, y, th, k, t)
  (if (k == 3L) negbin2_dexpected3_cpp(y, th[[1]], th[[2]], t)
    else negbin2_dexpected4_cpp(y, th[[1]], th[[2]], t)))
