#' @include betabinom1_distrib.R dexpected_support.R d34expected_analytic.R
NULL

# The third and fourth derivatives of betabinom1's expected information,
# compiled in src/d34e_betabinom1.cpp (direct sums at every sigma).
register_d34expected(BetaBinom1Distrib, function(d, y, th, k, t)
  (if (k == 3L) betabinom1_dexpected3_cpp(y, th[[1]], th[[2]], d@size, t)
    else betabinom1_dexpected4_cpp(y, th[[1]], th[[2]], d@size, t)))
