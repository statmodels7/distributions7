#' @include lognormal1_distrib.R dexpected_families.R d34expected_analytic.R
NULL

# The third and fourth derivatives of lognormal1's expected information, from
# the gaussian2 kernels: its (mu, sigma2) information is gaussian2's.
register_d34expected(Lognormal1Distrib, function(d, y, th, k, t)
  (if (k == 3L) gaussian2_dexpected3_cpp(y, th[[1]], th[[2]], t)
    else gaussian2_dexpected4_cpp(y, th[[1]], th[[2]], t)))
