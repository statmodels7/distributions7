#' @include gengamma2_distrib.R d34expected_analytic.R
NULL

# The third and fourth derivatives of gengamma2's expected information,
# compiled in src/d34e_gengamma2.cpp.
register_d34expected(GenGamma2Distrib, function(d, y, th, k, t)
  (if (k == 3L) gengamma2_dexpected3_cpp(y, th[[1]], th[[2]], th[[3]], t)
    else gengamma2_dexpected4_cpp(y, th[[1]], th[[2]], th[[3]], t)))
