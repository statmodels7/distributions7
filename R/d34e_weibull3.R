#' @include weibull3_distrib.R dexpected_families.R d34expected_analytic.R
NULL

# The third and fourth derivatives of weibull3's expected information,
# compiled in src/d34e_weibull3.cpp.
register_d34expected(Weibull3Distrib, function(d, y, th, k, t)
  (if (k == 3L) weibull3_dexpected3_cpp(y, th[[1]], th[[2]], t)
    else weibull3_dexpected4_cpp(y, th[[1]], th[[2]], t)))
