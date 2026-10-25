#' @include gumbel_distrib.R dexpected_families.R d34expected_analytic.R
NULL

# The third and fourth derivatives of gumbel's expected information,
# compiled in src/d34e_gumbel.cpp.
register_d34expected(GumbelDistrib, function(d, y, th, k, t)
  (if (k == 3L) gumbel_dexpected3_cpp(y, th[[1]], th[[2]], t)
    else gumbel_dexpected4_cpp(y, th[[1]], th[[2]], t)))
