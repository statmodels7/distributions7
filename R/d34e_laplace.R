#' @include laplace_distrib.R dexpected_families.R d34expected_analytic.R
NULL

# The third and fourth derivatives of laplace's expected information,
# compiled in src/d34e_laplace.cpp.
register_d34expected(LaplaceDistrib, function(d, y, th, k, t)
  (if (k == 3L) laplace_dexpected3_cpp(y, th[[1]], th[[2]], t)
    else laplace_dexpected4_cpp(y, th[[1]], th[[2]], t)))
