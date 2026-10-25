#' @include gengamma1_distrib.R d34expected_analytic.R
NULL

# The third and fourth derivatives of gengamma1's expected information,
# compiled in src/d34e_gengamma1.cpp.
register_d34expected(GenGamma1Distrib, function(d, y, th, k, t)
  (if (k == 3L) gengamma1_dexpected3_cpp(y, th[[1]], th[[2]], th[[3]], t)
    else gengamma1_dexpected4_cpp(y, th[[1]], th[[2]], th[[3]], t)))
