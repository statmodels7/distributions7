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
#' @seealso [numericals7::set_partitions()]
#' @keywords internal
index_partitions <- function(idx) {
  lapply(numericals7::set_partitions(length(idx)), function(p) lapply(p, function(b) idx[b]))
}

#' Canonical Component Name of a Block
#'
#' @description
#' Names a block of parameters the way [deriv_names()] does: in the
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

#' @title Complete Bell Polynomial in the Parent's Log-Derivatives
#'
#' @description
#' Computes \eqn{d^I f / f} from the derivatives of \eqn{\log f}, as the sum
#' over set partitions
#' \deqn{\frac{d^I f}{f} = \sum_{\pi} \prod_{B \in \pi} \ell^{(B)},}
#' where \eqn{\pi} runs over the partitions of the multi-index \eqn{I} and
#' \eqn{\ell^{(B)}} is the derivative of the log-density in the parameters that
#' block names.
#'
#' @details
#' This is the Bartlett lemma read backwards. The identity is normally used to
#' eliminate a derivative; here it builds one, and that is what every wrapper
#' needs: each wrapper's log-likelihood is the parent's log-density plus, or in
#' place of, \eqn{\log L} for some \eqn{\theta}-dependent \eqn{L}, and the
#' orders above two are assembled from this sum together with [log_deriv()].
#'
#' At order 1 it returns \eqn{\ell^{(i)}} and at order 2
#' \eqn{\ell^{(ij)} + \ell^{(i)}\ell^{(j)}}, the ordinary relation between the
#' second derivative of a density and of its logarithm. Those two cases
#' reproduce the hand-written closed forms exactly, and that agreement is the
#' license for the orders that have nothing to compare against.
#'
#' @param idx A character vector of parameter names, with repetition, naming
#'   the multi-index \eqn{I}: `c("mu", "mu", "sigma")` is
#'   \eqn{\partial^3/\partial\mu^2\partial\sigma}. Its length is the order.
#' @param ell A function of one block, returning \eqn{\ell^{(B)}} for the
#'   parameters that block names, as a numeric vector. Called once per block of
#'   every partition, so a caller with an expensive parent memoizes it.
#'
#' @return A numeric vector, the length of whatever `ell` returns: the ratio
#'   \eqn{d^I f / f} at each observation.
#'
#' @seealso [log_deriv()], the companion identity for a logarithm;
#'   [index_partitions()], which supplies the partitions;
#'   [numericals7::set_partitions()] for the enumeration itself.
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

#' @title Derivatives of a Logarithm From the Ratios Alone
#'
#' @description
#' Computes \eqn{d^I \log L} from the ratios \eqn{d^B L / L}, as
#' \deqn{d^I \log L = \sum_{\pi} (-1)^{|\pi|-1}(|\pi|-1)!
#'       \prod_{B \in \pi} \frac{d^B L}{L},}
#' the moment-to-cumulant relation. Only the **ratios** are needed, never
#' \eqn{L}'s own derivatives and never \eqn{L} itself, so a wrapper can
#' differentiate a normalizing constant it can only evaluate up to scale.
#'
#' @details
#' It is the inverse of [bell_f_ratio()]: that sum carries derivatives of a
#' logarithm to derivatives of the function, this one carries them back. At
#' order 1 it returns \eqn{d_i L / L} and at order 2
#' \eqn{d_{ij}L/L - (d_i L/L)(d_j L/L)}, the familiar relation for the second
#' derivative of a logarithm.
#'
#' @param idx A character vector of parameter names, with repetition, naming
#'   the multi-index \eqn{I}. Its length is the order.
#' @param ratio A function of one block, returning \eqn{d^B L / L} for that
#'   block as a numeric vector. Called once per block of every partition.
#'
#' @return A numeric vector, the length of whatever `ratio` returns:
#'   \eqn{d^I \log L} at each observation.
#'
#' @seealso [bell_f_ratio()], the companion identity;
#'   [index_partitions()], which supplies the partitions;
#'   [truncated()] and [zero_inflated()], the wrappers that consume it.
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
#' Fetches the parent's derivatives up to `max_order` and returns a function
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
#' @seealso [canon_key()]
#' @keywords internal
parent_ell <- function(parent, y, theta, max_order, params) {
  d <- vector("list", max_order)
  d[[1]] <- distrib_gradient(parent, y, theta)
  if (max_order >= 2) d[[2]] <- distrib_hessian(parent, y, theta)
  if (max_order >= 3) d[[3]] <- distrib_deriv3(parent, y, theta)
  if (max_order >= 4) d[[4]] <- distrib_deriv4(parent, y, theta)
  function(block) d[[length(block)]][[canon_key(block, params)]]
}
