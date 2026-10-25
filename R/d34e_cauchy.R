#' @include cauchy_distrib.R dexpected_families.R d34expected_analytic.R
NULL

# The third and fourth derivatives of cauchy's expected information,
# compiled in src/d34e_cauchy.cpp.
register_d34expected(CauchyDistrib, function(d, y, th, k, t)
  (if (k == 3L) cauchy_dexpected3_cpp(y, th[[1]], th[[2]], t)
    else cauchy_dexpected4_cpp(y, th[[1]], th[[2]], t)))
