#' @include betabinom2_distrib.R dexpected_support.R d34expected_analytic.R
NULL

# The third and fourth derivatives of betabinom2's expected information,
# compiled in src/d34e_betabinom2.cpp.
register_d34expected(BetaBinom2Distrib, function(d, y, th, k, t)
  (if (k == 3L) betabinom2_dexpected3_cpp(y, th[[1]], th[[2]], d@size, t)
    else betabinom2_dexpected4_cpp(y, th[[1]], th[[2]], d@size, t)))
