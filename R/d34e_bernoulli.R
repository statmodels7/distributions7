#' @include bernoulli_distrib.R dexpected_families.R d34expected_analytic.R
NULL

# The third and fourth derivatives of bernoulli's expected information,
# compiled in src/d34e_bernoulli.cpp.
register_d34expected(BernoulliDistrib, function(d, y, th, k, t)
  (if (k == 3L) bernoulli_dexpected3_cpp(y, th[[1]], t)
    else bernoulli_dexpected4_cpp(y, th[[1]], t)))
