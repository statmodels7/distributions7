#' @include vonmises2_distrib.R dexpected_vonmises.R d34expected_analytic.R
NULL

# The third and fourth derivatives of vonmises2's expected information,
# compiled in src/d34e_vonmises2.cpp.
register_d34expected(VonMises2Distrib, function(d, y, th, k, t)
  (if (k == 3L) vonmises2_dexpected3_cpp(y, th[[1]], th[[2]], t)
    else vonmises2_dexpected4_cpp(y, th[[1]], th[[2]], t)))
