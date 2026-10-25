#' @include pig1_distrib.R pig_expected.R d34expected_analytic.R
NULL

# The third and fourth derivatives of pig1's expected information,
# compiled in src/d34e_pig1.cpp.
register_d34expected(Pig1Distrib, function(d, y, th, k, t)
  (if (k == 3L) pig1_dexpected3_cpp(y, th[[1]], th[[2]], t)
    else pig1_dexpected4_cpp(y, th[[1]], th[[2]], t)))
