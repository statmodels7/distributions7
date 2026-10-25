#' @include pig2_distrib.R pig_expected.R d34expected_analytic.R
NULL

# The third and fourth derivatives of pig2's expected information,
# compiled in src/d34e_pig2.cpp.
register_d34expected(Pig2Distrib, function(d, y, th, k, t)
  (if (k == 3L) pig2_dexpected3_cpp(y, th[[1]], th[[2]], t)
    else pig2_dexpected4_cpp(y, th[[1]], th[[2]], t)))
