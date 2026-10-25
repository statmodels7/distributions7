#' @include binomial_distrib.R dexpected_families.R d34expected_analytic.R
NULL

# The third and fourth derivatives of binomial's expected information,
# compiled in src/d34e_binomial.cpp.
register_d34expected(BinomialDistrib, function(d, y, th, k, t)
  (if (k == 3L) binomial_dexpected3_cpp(y, th[[1]], d@size, t)
    else binomial_dexpected4_cpp(y, th[[1]], d@size, t)))
