#' @include invgauss2_distrib.R dexpected_families.R d34expected_analytic.R
NULL

# The third and fourth derivatives of invgauss2's expected information,
# compiled in src/d34e_invgauss2.cpp.
register_d34expected(InvGauss2Distrib, function(d, y, th, k, t)
  (if (k == 3L) invgauss2_dexpected3_cpp(y, th[[1]], th[[2]], t)
    else invgauss2_dexpected4_cpp(y, th[[1]], th[[2]], t)))
