#' @include vonmises1_distrib.R dexpected_vonmises.R d34expected_analytic.R
NULL

# The third and fourth derivatives of vonmises1's expected information,
# compiled in src/d34e_vonmises1.cpp.
register_d34expected(VonMises1Distrib, function(d, y, th, k, t)
  (if (k == 3L) vonmises1_dexpected3_cpp(y, th[[1]], th[[2]], t)
    else vonmises1_dexpected4_cpp(y, th[[1]], th[[2]], t)))
