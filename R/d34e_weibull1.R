#' @include weibull1_distrib.R dexpected_families.R d34expected_analytic.R
NULL

# The third and fourth derivatives of weibull1's expected information,
# compiled in src/d34e_weibull1.cpp.
register_d34expected(Weibull1Distrib, function(d, y, th, k, t)
  (if (k == 3L) weibull1_dexpected3_cpp(y, th[[1]], th[[2]], t)
    else weibull1_dexpected4_cpp(y, th[[1]], th[[2]], t)))
