#' @include cdf_higher.R cdf_mapped_higher.R cdf_compiled.R gamma1_distrib.R
#' @include gamma2_distrib.R chisq_distrib.R gengamma1_distrib.R
#' @include gengamma2_distrib.R beta1_distrib.R beta2_distrib.R
#' @include vonmises1_distrib.R vonmises2_distrib.R
NULL

#' Higher Log-CDF Derivatives by the Compiled Quadrature
#'
#' @description
#' Builds the [distrib_deriv3_cdf()] or [distrib_deriv4_cdf()] body of the
#' families whose compiled distribution function integrates its shape
#' derivatives and which are not location-scale: the derivatives of \eqn{F}
#' of every order up to `order` from [cdf_tables()], whose third and fourth
#' orders are the exact integrals of [compiled_cdf_deriv_k()], put on the
#' requested tail and scale by [cdf_scale_k()].
#'
#' @param order The derivative order, 3 or 4.
#'
#' @return A function of `(distrib, q, theta, lower.tail, log, ...)` suitable
#'   for registering as an S7 method on either generic.
#'
#' @aliases distrib_deriv3_cdf.Gamma1Distrib distrib_deriv4_cdf.Gamma1Distrib
#'   distrib_deriv3_cdf.Gamma2Distrib distrib_deriv4_cdf.Gamma2Distrib
#'   distrib_deriv3_cdf.ChisqDistrib distrib_deriv4_cdf.ChisqDistrib
#'   distrib_deriv3_cdf.GenGamma1Distrib distrib_deriv4_cdf.GenGamma1Distrib
#'   distrib_deriv3_cdf.GenGamma2Distrib distrib_deriv4_cdf.GenGamma2Distrib
#'   distrib_deriv3_cdf.Beta1Distrib distrib_deriv4_cdf.Beta1Distrib
#'   distrib_deriv3_cdf.Beta2Distrib distrib_deriv4_cdf.Beta2Distrib
#'   distrib_deriv3_cdf.VonMises1Distrib distrib_deriv4_cdf.VonMises1Distrib
#'   distrib_deriv3_cdf.VonMises2Distrib distrib_deriv4_cdf.VonMises2Distrib
#'
#' @keywords internal
compiled_deriv_cdf_k <- function(order) {
  function(distrib, q, theta, lower.tail = TRUE, log = TRUE, ...) {
    cdf_scale_k(distrib, distrib_cdf(distrib, q, theta),
                cdf_tables(distrib, q, theta, order), order, lower.tail, log)
  }
}

for (.cls in list(Gamma1Distrib, Gamma2Distrib, ChisqDistrib, GenGamma1Distrib,
                  GenGamma2Distrib, Beta1Distrib, Beta2Distrib,
                  VonMises1Distrib, VonMises2Distrib)) {
  S7::method(distrib_deriv3_cdf, .cls) <- compiled_deriv_cdf_k(3L)
  S7::method(distrib_deriv4_cdf, .cls) <- compiled_deriv_cdf_k(4L)
}
rm(.cls)
