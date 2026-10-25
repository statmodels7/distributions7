#' @include invgauss1_distrib.R dexpected_families.R d34expected_analytic.R
NULL

# The third and fourth derivatives of invgauss1's expected information,
# compiled in src/d34e_invgauss1.cpp.
register_d34expected(InvGauss1Distrib, function(d, y, th, k, t)
  (if (k == 3L) invgauss1_dexpected3_cpp(y, th[[1]], th[[2]], t)
    else invgauss1_dexpected4_cpp(y, th[[1]], th[[2]], t)))
