#' Mixed Derivatives by One Stencil on the Highest Analytical Quantity
#'
#' @description
#' Computes the derivatives of the log-density of order `ay` in the response
#' and `ctheta` in the parameters, for a continuous family that has no method
#' of its own for them. Among the quantities that the family implements itself
#' (the log-density, its response derivatives, its score and Hessian in the
#' parameters, and the mixed derivatives of lower order), the one that leaves
#' the smallest order to difference is taken, and each component is a single
#' tensor-product central difference of it, never a difference of a
#' difference.
#'
#' @details
#' The computation is [numericals7::fd_tensor()] over the coordinates
#' \eqn{(y, \theta_1, \dots, \theta_p)}. A candidate quantity of order
#' \eqn{(a, c)} in the response and the parameters is declared to it where
#' \eqn{a \le} `ay`, \eqn{c \le} `ctheta`, \eqn{(a, c)} is not the target
#' itself, and the family's class registers its method, as [owns_method()]
#' reads it, so that a numerical fallback is never differenced again. The
#' log-density is always available. The candidate with the smallest residual
#' order \eqn{r = (\mathrm{ay} - a) + (\mathrm{ctheta} - c)} is used, and
#' among equal ones the candidate with more derivatives in the parameters,
#' which declares more multi-indices, is used; the candidates are declared in
#' decreasing \eqn{c}, so that the order of declaration gives the same
#' preference when there is a single parameter.
#'
#' The step is \eqn{h = \varepsilon^{1/(r+a)}\max(1, |x|)} in every
#' direction, with \eqn{a} the accuracy, reduced so that every node stays
#' inside the support in the response and inside the bounds of each
#' parameter. Values shared by several components are computed once.
#'
#' For [transformation()], which carries the parent's score and Hessian in
#' the parameters and no response derivative, the fourth-order
#' \eqn{\partial^4\ell/\partial y^2\,\partial\theta^2} is a second difference
#' in the response of the analytical Hessian, so \eqn{r = 2}.
#'
#' @param distrib A `continuous_distrib` object.
#' @param y A numeric vector of observations.
#' @param theta A named list of parameters, aligned by the generic.
#' @param ay The order in the response, 1 to 3.
#' @param ctheta The order in the parameters, 1 or 2.
#' @param accuracy The order of accuracy of the one-dimensional stencils,
#'   passed to [numericals7::fd_tensor()].
#'
#' @return A named list: one vector per parameter, keyed by `distrib@params`,
#'   where `ctheta` is 1; one vector per unordered pair, keyed as
#'   [`hess_names(distrib@params)`][hess_names], where it is 2.
#'
#' @seealso [tensor_derivatives()] for the derivatives in the parameters only,
#'   [distrib_cross_y()], [distrib_cross2_y()], [distrib_cross3_y()],
#'   [distrib_grad_y_hess()] and [distrib_hess_y_hess()], whose continuous
#'   fallbacks call this.
#' @keywords internal
#'
#' @examples
#' # transformation() has no response derivative of its own: the mixed
#' # fourth order of a log-gamma is a stencil on its analytical Hessian.
#' d <- fixed(transformation(gamma2_distrib(), log_transform()), mu = 1)
#' b <- c(-0.4, 0.3)
#' k <- 1 / 0.26
#' got <- distributions7:::mixed_tensor_derivatives(d, b, list(sigma2 = 0.26), 2L, 2L)
#' rbind(got$sigma2_sigma2, -2 * k^3 * exp(b))
mixed_tensor_derivatives <- function(distrib, y, theta, ay, ctheta,
                                     accuracy = 4L) {
  ay <- as.integer(ay)
  ctheta <- as.integer(ctheta)
  params <- distrib@params
  p <- length(params)
  # the parameter multi-indices of order 1 and 2, in the order the generics
  # key their results: params, then hess_names()
  pidx <- function(c) {
    if (c == 0L) return(matrix(0L, 1L, p))
    t(vapply(deriv_index_list(p, c), tabulate, integer(p), nbins = p))
  }
  rows_of <- function(a, c) {
    m <- pidx(c)
    if (p == 1L && c > 0L) m <- matrix(m, ncol = 1L)
    cbind(a, m, deparse.level = 0L)
  }
  keys_of <- function(c) if (c == 1L) params else hess_names(params)

  # the candidates, in decreasing order in theta (see the details)
  cands <- list(
    list(gen = distrib_grad_y_hess, a = 1L, c = 2L),
    list(gen = distrib_hessian, a = 0L, c = 2L),
    list(gen = distrib_cross3_y, a = 3L, c = 1L),
    list(gen = distrib_cross2_y, a = 2L, c = 1L),
    list(gen = distrib_cross_y, a = 1L, c = 1L),
    list(gen = distrib_gradient, a = 0L, c = 1L),
    list(gen = distrib_deriv4_y, a = 4L, c = 0L),
    list(gen = distrib_deriv3_y, a = 3L, c = 0L),
    list(gen = distrib_hess_y, a = 2L, c = 0L),
    list(gen = distrib_grad_y, a = 1L, c = 0L))
  usable <- Filter(function(b) {
    b$a <= ay && b$c <= ctheta && !(b$a == ay && b$c == ctheta) &&
      owns_method(distrib, b$gen)
  }, cands)

  cs <- tensor_coords(distrib, theta, "parameter")
  x <- c(list(.y = y), cs$x)
  bounds <- c(list(.y = distrib@bounds), cs$bounds)
  at <- function(x) cs$to_theta(x[-1L])
  known <- lapply(usable, function(b) {
    gen <- b$gen
    fun <- if (b$c == 0L) {
      function(x) {
        v <- gen(distrib, x$.y, at(x))
        if (is.list(v)) v[[1L]] else v
      }
    } else {
      keys <- keys_of(b$c)
      function(x) unname(gen(distrib, x$.y, at(x))[keys])
    }
    numericals7::fd_known(fun, indices = rows_of(b$a, b$c))
  })
  f <- function(x) distrib_pdf(distrib, x$.y, at(x), log = TRUE)
  got <- numericals7::fd_tensor(f, x, components = rows_of(ay, ctheta),
                                known = known, bounds = bounds,
                                accuracy = accuracy)
  out <- unname(got)
  names(out) <- keys_of(ctheta)
  attr(out, "fd_source") <- NULL
  out
}
