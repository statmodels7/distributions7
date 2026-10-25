#' @include chisq_distrib.R dexpected_families.R d34expected_analytic.R
NULL

# The third and fourth derivatives of chisq's expected information,
# compiled in src/d34e_chisq.cpp.
register_d34expected(ChisqDistrib, function(d, y, th, k, t)
  (if (k == 3L) chisq_dexpected3_cpp(y, th[[1]], t)
    else chisq_dexpected4_cpp(y, th[[1]], t)))
