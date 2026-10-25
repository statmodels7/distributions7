#' @include gpd_distrib.R dexpected_families.R d34expected_analytic.R
NULL

# The third and fourth derivatives of gpd's expected information,
# compiled in src/d34e_gpd.cpp.
register_d34expected(GPDDistrib, function(d, y, th, k, t)
  (if (k == 3L) gpd_dexpected3_cpp(y, th[[1]], th[[2]], t)
    else gpd_dexpected4_cpp(y, th[[1]], th[[2]], t)))
