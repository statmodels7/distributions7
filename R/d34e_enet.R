#' @include enet_distrib.R dexpected_families.R d34expected_analytic.R
NULL

# The third and fourth derivatives of enet's expected information,
# compiled in src/d34e_enet.cpp.
register_d34expected(EnetDistrib, function(d, y, th, k, t)
  (if (k == 3L) enet_dexpected3_cpp(y, th[[1]], th[[2]], th[[3]], t)
    else enet_dexpected4_cpp(y, th[[1]], th[[2]], th[[3]], t)))
