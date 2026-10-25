#' @include gaussian3_distrib.R dexpected_families.R d34expected_analytic.R
NULL

# The third and fourth derivatives of gaussian3's expected information,
# compiled in src/d34e_gaussian3.cpp.
register_d34expected(Gaussian3Distrib, function(d, y, th, k, t)
  (if (k == 3L) gaussian3_dexpected3_cpp(y, th[[1]], th[[2]], t)
    else gaussian3_dexpected4_cpp(y, th[[1]], th[[2]], t)))
