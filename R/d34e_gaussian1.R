#' @include gaussian1_distrib.R d34expected_analytic.R
NULL

# The third and fourth derivatives of gaussian1's expected information,
# compiled in src/d34e_gaussian1.cpp.
register_d34expected(Gaussian1Distrib, function(d, y, th, k, t)
  (if (k == 3L) gaussian1_dexpected3_cpp(y, th[[1]], th[[2]], t)
    else gaussian1_dexpected4_cpp(y, th[[1]], th[[2]], t)))
