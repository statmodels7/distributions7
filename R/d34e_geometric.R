#' @include geometric_distrib.R dexpected_families.R d34expected_analytic.R
NULL

# The third and fourth derivatives of geometric's expected information,
# compiled in src/d34e_geometric.cpp.
register_d34expected(GeometricDistrib, function(d, y, th, k, t)
  (if (k == 3L) geometric_dexpected3_cpp(y, th[[1]], t)
    else geometric_dexpected4_cpp(y, th[[1]], t)))
