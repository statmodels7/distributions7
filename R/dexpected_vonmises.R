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
#' \eqn{A(\kappa) = I_1(\kappa)/I_0(\kappa)} and of its inverse, computed by
#' one compiled kernel per order from numericals7's compiled Bessel ratio.
#'
#' @details
#' By the concentration, \eqn{\mathbb{E}[\ell_{\mu\mu}] = -\kappa A(\kappa)} and
#' \eqn{\mathbb{E}[\ell_{\kappa\kappa}] = -A'(\kappa)}, so the derivatives in
#' \eqn{\kappa} read \eqn{A'}, \eqn{A''} and \eqn{A'''}. By the mean resultant
#' length \eqn{\rho = A(\kappa)}, \eqn{\mathbb{E}[\ell_{\mu\mu}] = -\rho\,\kappa(\rho)}
#' and \eqn{\mathbb{E}[\ell_{\rho\rho}] = -\kappa'(\rho)}, so they read the
#' inverse's derivatives \eqn{\kappa'}, \eqn{\kappa''} and \eqn{\kappa'''}.
#' Nothing moves with the direction \eqn{\mu}. The first-order kernel reads
#' \eqn{A} to \eqn{A''} (or \eqn{\kappa} to \eqn{\kappa''}), the second-order
#' one \eqn{A'} to \eqn{A'''} (or \eqn{\kappa'} to \eqn{\kappa'''}).
#'
#' @param distrib A von Mises distribution object.
#' @param y A numeric vector of observations, read for its length.
#' @param theta A named list of parameters.
#' @param scale `"parameter"` or `"link"`.
#' @param approx,nsim Unused.
#' @param ... Unused.
#' @param threads A single positive integer, how many threads the kernel may
#'   use.
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
#' @param f The derivative of order `order` of \eqn{f}.
#' @param g The derivative of order `order` of \eqn{g}.
#'
#' @return A named list keyed as [dexpected_names()] or [d2expected_names()].
#'
#' @keywords internal
vm_dexpected <- function(P, n, order, f, g) {
  r <- function(v) rep_len(v, n)
  z <- r(0)
  if (order == 1L) {
    out <- list(z, r(f), z, r(g), z, z)
    names(out) <- c(dexpected_key(P, 1, 1, 1), dexpected_key(P, 1, 1, 2),
                    dexpected_key(P, 2, 2, 1), dexpected_key(P, 2, 2, 2),
                    dexpected_key(P, 1, 2, 1), dexpected_key(P, 1, 2, 2))
    return(out)
  }
  pr <- list(c(1, 1), c(2, 2), c(1, 2))
  out <- list()
  for (ab in pr) for (cd in pr) {
    v <- if (identical(cd, c(2, 2)) && identical(ab, c(1, 1))) f else
      if (identical(cd, c(2, 2)) && identical(ab, c(2, 2))) g else 0
    out[[d2expected_key(P, ab[1], ab[2], cd[1], cd[2])]] <- r(v)
  }
  out
}

register_dexpected(VonMises1Distrib, function(d, y, th, k, t) {
  # E_mm = -k A: (-(A + k A'), -(2 A' + k A'')); E_kk = -A': (-A'', -A''')
  fg <- if (k == 1L) vonmises1_dexpected1_cpp(th[[2]], t) else
    vonmises1_dexpected2_cpp(th[[2]], t)
  vm_dexpected(d@params, length(y), k, fg$f, fg$g)
})

register_dexpected(VonMises2Distrib, function(d, y, th, k, t) {
  # E_mm = -rho kappa(rho): (-(kappa + rho kappa'), -(2 kappa' + rho kappa''));
  # E_rr = -kappa'(rho): (-kappa'', -kappa''')
  fg <- if (k == 1L) vonmises2_dexpected1_cpp(th[[2]], t) else
    vonmises2_dexpected2_cpp(th[[2]], t)
  vm_dexpected(d@params, length(y), k, fg$f, fg$g)
})
