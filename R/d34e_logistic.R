#' @include logistic_distrib.R dexpected_families.R d34expected_analytic.R
NULL

# The third and fourth derivatives of logistic's expected information,
# compiled in src/d34e_logistic.cpp.
register_d34expected(LogisticDistrib, function(d, y, th, k, t)
  (if (k == 3L) logistic_dexpected3_cpp(y, th[[1]], th[[2]], t)
    else logistic_dexpected4_cpp(y, th[[1]], th[[2]], t)))
