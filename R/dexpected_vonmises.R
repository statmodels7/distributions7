#' @include dexpected_hessian.R dexpected_families.R vonmises1_distrib.R
#' @include vonmises2_distrib.R
NULL

#' @title Derivatives of the Expected Information, von Mises
#' @name distrib_dexpected_hessian.VonMises1Distrib
#' @aliases distrib_d2expected_hessian.VonMises1Distrib
#'   distrib_dexpected_hessian.VonMises2Distrib
#'   distrib_d2expected_hessian.VonMises2Distrib
#' @description
#' The first and second derivatives of the expected information in the
#' parameters, from the derivatives of the Bessel ratio
#' \eqn{A(\kappa) = I_1(\kappa)/I_0(\kappa)} that
#' [numericals7::bessel_i_ratio_derivs()] and
#' [numericals7::bessel_i_ratio_inverse()] return.
#'
#' @details
#' By the concentration, \eqn{\mathbb{E}[\ell_{\mu\mu}] = -\kappa A(\kappa)} and
#' \eqn{\mathbb{E}[\ell_{\kappa\kappa}] = -A'(\kappa)}, so the derivatives in
#' \eqn{\kappa} read \eqn{A'}, \eqn{A''} and \eqn{A'''}. By the mean resultant
#' length \eqn{\rho = A(\kappa)}, \eqn{\mathbb{E}[\ell_{\mu\mu}] = -\rho\,\kappa(\rho)}
#' and \eqn{\mathbb{E}[\ell_{\rho\rho}] = -\kappa'(\rho)}, so they read the
#' inverse's derivatives \eqn{\kappa'}, \eqn{\kappa''} and \eqn{\kappa'''}.
#' Nothing moves with the direction \eqn{\mu}. The expected information of
#' either family is itself computed in R from the same functions, which is why
#' these are too.
#'
#' @param distrib A von Mises distribution object.
#' @param y A numeric vector of observations, read for its length.
#' @param theta A named list of parameters.
#' @param scale `"parameter"` or `"link"`.
#' @param approx,nsim Unused.
#' @param ... Unused.
#' @param threads Unused; accepted for the shared signature.
#' @return A named list keyed as [dexpected_names()] or [d2expected_names()].
#' @seealso [distrib_dexpected_hessian()], [distrib_d2expected_hessian()]
#' @keywords internal
NULL

#' The von Mises Families' Derivative Tables
#'
#' @description
#' Assembles the keyed components for a two-parameter family in which nothing
#' moves with the first parameter, \eqn{E_{11} = f(c)},
#' \eqn{E_{22} = g(c)} and \eqn{E_{12} = 0} in the second parameter \eqn{c}.
#'
#' @param P The two parameter names.
#' @param n The length of the result.
#' @param order `1L` or `2L`.
#' @param f1,f2 The first and second derivatives of \eqn{f}.
#' @param g1,g2 The first and second derivatives of \eqn{g}.
#'
#' @return A named list keyed as [dexpected_names()] or [d2expected_names()].
#'
#' @keywords internal
vm_dexpected <- function(P, n, order, f1, f2, g1, g2) {
  r <- function(v) rep_len(v, n)
  z <- r(0)
  if (order == 1L) {
    out <- list(z, r(f1), z, r(g1), z, z)
    names(out) <- c(dexpected_key(P, 1, 1, 1), dexpected_key(P, 1, 1, 2),
                    dexpected_key(P, 2, 2, 1), dexpected_key(P, 2, 2, 2),
                    dexpected_key(P, 1, 2, 1), dexpected_key(P, 1, 2, 2))
    return(out)
  }
  pr <- list(c(1, 1), c(2, 2), c(1, 2))
  out <- list()
  for (ab in pr) for (cd in pr) {
    v <- if (identical(cd, c(2, 2)) && identical(ab, c(1, 1))) f2 else
      if (identical(cd, c(2, 2)) && identical(ab, c(2, 2))) g2 else 0
    out[[d2expected_key(P, ab[1], ab[2], cd[1], cd[2])]] <- r(v)
  }
  out
}

register_dexpected(VonMises1Distrib, function(d, y, th, k, t) {
  a <- numericals7::bessel_i_ratio_derivs(th[[2]])
  kap <- th[[2]]
  # E_mm = -k A: (-(A + k A'), -(2 A' + k A'')); E_kk = -A': (-A'', -A''')
  vm_dexpected(d@params, length(y), k,
               -(a$A + kap * a$d1), -(2 * a$d1 + kap * a$d2), -a$d2, -a$d3)
})

register_dexpected(VonMises2Distrib, function(d, y, th, k, t) {
  kd <- numericals7::bessel_i_ratio_inverse(th[[2]])
  rho <- th[[2]]
  # E_mm = -rho kappa(rho): (-(kappa + rho kappa'), -(2 kappa' + rho kappa''));
  # E_rr = -kappa'(rho): (-kappa'', -kappa''')
  vm_dexpected(d@params, length(y), k,
               -(kd$kappa + rho * kd$d1), -(2 * kd$d1 + rho * kd$d2),
               -kd$d2, -kd$d3)
})
