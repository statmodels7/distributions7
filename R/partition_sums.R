#' @include generics.R
NULL

# The two partition sums the toolkit is built on, and the small helpers they
# need. They live here rather than beside the wrappers that first used them
# because they are general: the wrappers assemble a parent's derivatives with
# them, and so do the higher derivatives of the distribution function.
#
#   bell_f_ratio()  d^I f / f from the derivatives of log f -- the complete
#                   Bell polynomial, the Bartlett lemma read backwards
#   log_deriv()     d^I log L from the ratios d^B L / L alone -- the
#                   moment-to-cumulant relation

#' Set Partitions of a Multi-Index
#'
#' @description
#' Every way of splitting a multi-index into blocks, with the parameter names
#' carried through rather than the positions.
#'
#' @param idx A character vector of parameter names, with repetition.
#'
#' @return A list of partitions, each a list of character-vector blocks.
#'
#' @seealso \code{\link[numericals7]{set_partitions}}
#' @keywords internal
index_partitions <- function(idx) {
  lapply(numericals7::set_partitions(length(idx)), function(p) lapply(p, function(b) idx[b]))
}

#' Canonical Component Name of a Block
#'
#' @description
#' Names a block of parameters the way \code{\link{deriv_names}} does: in the
#' order the distribution declares them, joined by an underscore.
#'
#' @param block A character vector of parameter names.
#' @param params The distribution's parameter names, in declaration order.
#'
#' @return A single string.
#'
#' @keywords internal
canon_key <- function(block, params) {
  paste(block[order(match(block, params))], collapse = "_")
}

#' Complete Bell Polynomial in the Parent's Log-Derivatives
#'
#' @description
#' Computes \eqn{d^I f / f} from the derivatives of \eqn{\log f}, as the sum over
#' set partitions \eqn{\sum_\pi \prod_{B \in \pi} \ell^{(B)}}.
#'
#' @details
#' This is the Bartlett lemma read backwards: instead of using the identity to
#' eliminate a derivative, it uses it to build one. Every wrapper needs it,
#' because each of their log-likelihoods is the parent's log-density plus, or
#' instead of, \eqn{\log L} for some \eqn{\theta}-dependent \eqn{L}.
#'
#' @param idx A character vector of parameter names, with repetition.
#' @param ell A function returning the parent's log-derivative for a block.
#'
#' @return A numeric vector.
#'
#' @seealso \code{\link{log_deriv}}, the companion identity.
#' @keywords internal
bell_f_ratio <- function(idx, ell) {
  total <- 0
  for (p in index_partitions(idx)) {
    term <- 1
    for (b in p) term <- term * ell(b)
    total <- total + term
  }
  total
}

#' Derivatives of a Logarithm From the Ratios Alone
#'
#' @description
#' Computes \eqn{d^I \log L} as
#' \eqn{\sum_\pi (-1)^{|\pi|-1}(|\pi|-1)! \prod_{B \in \pi} (d^B L / L)}.
#'
#' @details
#' The moment-to-cumulant relation. What makes it the right tool here is that
#' only the \strong{ratios} \eqn{d^B L / L} are needed, never \eqn{L}'s own
#' derivatives -- and the ratios are exactly what each wrapper can supply
#' cheaply, a truncated expectation for truncation and an affine expression for
#' the zero wrappers.
#'
#' @param idx A character vector of parameter names, with repetition.
#' @param ratio A function returning \eqn{d^B L / L} for a block.
#'
#' @return A numeric vector.
#'
#' @seealso \code{\link{bell_f_ratio}}, the companion identity.
#' @keywords internal
log_deriv <- function(idx, ratio) {
  total <- 0
  for (p in index_partitions(idx)) {
    k <- length(p)
    term <- (-1)^(k - 1) * gamma(k)          # (-1)^{k-1} (k-1)!
    for (b in p) term <- term * ratio(b)
    total <- total + term
  }
  total
}

#' Look Up the Parent's Derivative Components by Block
#'
#' @description
#' Fetches the parent's derivatives up to \code{max_order} and returns a function
#' giving the component belonging to any block.
#'
#' @details
#' The orders are fetched once per call rather than once per block. A partition
#' sum at fourth order asks for the same handful of components repeatedly, and
#' the parent's derivative may itself be expensive.
#'
#' @param parent The parent distribution.
#' @param y A numeric vector of observations.
#' @param theta A named list of the parent's parameters.
#' @param max_order The highest order needed, 1 to 4.
#' @param params The parent's parameter names, in declaration order.
#'
#' @return A function of one block, returning that component's vector.
#'
#' @seealso \code{\link{canon_key}}
#' @keywords internal
parent_ell <- function(parent, y, theta, max_order, params) {
  d <- vector("list", max_order)
  d[[1]] <- distrib_gradient(parent, y, theta)
  if (max_order >= 2) d[[2]] <- distrib_hessian(parent, y, theta)
  if (max_order >= 3) d[[3]] <- distrib_deriv3(parent, y, theta)
  if (max_order >= 4) d[[4]] <- distrib_deriv4(parent, y, theta)
  function(block) d[[length(block)]][[canon_key(block, params)]]
}
