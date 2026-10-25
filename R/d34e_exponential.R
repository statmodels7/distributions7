#' @include exponential_distrib.R dexpected_families.R d34expected_analytic.R
NULL

# The third and fourth derivatives of exponential's expected information,
# compiled in src/d34e_exponential.cpp.
register_d34expected(ExponentialDistrib, function(d, y, th, k, t)
  (if (k == 3L) exponential_dexpected3_cpp(y, th[[1]], t)
    else exponential_dexpected4_cpp(y, th[[1]], t)))
