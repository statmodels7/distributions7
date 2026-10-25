#' @include gaussian2_distrib.R dexpected_families.R d34expected_analytic.R
NULL

# The third and fourth derivatives of gaussian2's expected information,
# compiled in src/d34e_gaussian2.cpp.
register_d34expected(Gaussian2Distrib, function(d, y, th, k, t)
  (if (k == 3L) gaussian2_dexpected3_cpp(y, th[[1]], th[[2]], t)
    else gaussian2_dexpected4_cpp(y, th[[1]], th[[2]], t)))
