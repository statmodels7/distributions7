#' @include beta1_distrib.R d34expected_analytic.R
NULL

# The third and fourth derivatives of beta1's expected information,
# compiled in src/d34e_beta1.cpp.
register_d34expected(Beta1Distrib, function(d, y, th, k, t)
  (if (k == 3L) beta1_dexpected3_cpp(y, th[[1]], th[[2]], t)
    else beta1_dexpected4_cpp(y, th[[1]], th[[2]], t)))
