#' @include lognormal2_distrib.R d34expected_analytic.R
NULL

# The third and fourth derivatives of lognormal2's expected information,
# compiled in src/d34e_lognormal2.cpp.
register_d34expected(Lognormal2Distrib, function(d, y, th, k, t)
  (if (k == 3L) lognormal2_dexpected3_cpp(y, th[[1]], th[[2]], t)
    else lognormal2_dexpected4_cpp(y, th[[1]], th[[2]], t)))
