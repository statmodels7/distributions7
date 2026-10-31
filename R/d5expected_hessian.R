#' @include distrib.R generics.R utility_functions.R numerical_derivatives.R d34expected_hessian.R
NULL

# The fifth derivative of the expected information in the parameters.
#
# The consumer is a score-driven filter whose score is scaled by a power of the
# expected information: the fifth derivative of the filtered predictor reads
# the fifth derivative of the information. No family writes this order out.
# The default below is one tensor stencil (numericals7::fd_tensor()) on the
# highest order of the expected information that the family computes itself:
# the fourth for every shipped univariate family and for a wrapper of one, a
# first difference; the third or the second for a family defined by the user,
# a stencil of order two or three. A family without an analytic second
# derivative of its information has no fifth, as it has no third or fourth.

#' The Fifth Derivative of the Expected Information
#'
#' @description
#' Returns \eqn{\partial^5\,\mathbb{E}[\ell_{ab}]/\partial\theta_c\,
#' \partial\theta_d\,\partial\theta_e\,\partial\theta_f\,\partial\theta_g},
#' with one component per pair \eqn{(a,b)} and per unordered quintuple
#' \eqn{(c,d,e,f,g)}.
#'
#' @details
#' The components are symmetric in \eqn{(a,b)} and in \eqn{(c,d,e,f,g)}
#' separately, and are keyed by [d5expected_names()].
#'
#' No family computes this order in closed form. The default method applies
#' one tensor-product central stencil of [numericals7::fd_tensor()] to the
#' highest derivative of the expected information that the family computes
#' itself: a first difference of [distrib_d4expected_hessian()] for the
#' univariate families shipped with the package, which compute the fourth
#' order in closed form, and for a wrapper of such a family; a stencil of
#' order two on [distrib_d3expected_hessian()] or of order three on
#' [distrib_d2expected_hessian()] for a family that registers only those. A
#' wrapper never differences from an order above the one of the family that it
#' wraps, since its own orders are built from the parent's. A family without an
#' analytic second derivative of its expected information signals an error,
#' as it does for the third and fourth orders. The step of a stencil of order
#' \eqn{r} is \eqn{\varepsilon^{1/(r+4)}\max(1, |x|)} in the coordinate
#' \eqn{x}, reduced so that every node stays inside the domain.
#'
#' The method returns the parameter scale. On `scale = "link"` the generic
#' composes the result with the expected information and its derivatives of
#' orders one to four on the parameter scale, by Leibniz's rule and the
#' diagonal Faa di Bruno formula of [dexpected_link_higher()], which reads the
#' derivatives of the inverse links up to the sixth.
#'
#' With an object as `method`, the result is the derivative of the matrix that
#' the method defines, computed as for the lower orders by Leibniz's formula
#' for [expect_opg()] and with the moved density for [expect_integrate()] and
#' [expect_mc()]. It reads observed derivatives up to the sixth
#' ([distrib_deriv6()]), so only `identity = "score"` is available; the
#' Hessian identity would need the seventh order and signals an error.
#'
#' @param distrib A distribution object inheriting from `distrib`.
#' @param y A numeric vector of observations.
#' @param theta A named list of parameters; each has length 1 or
#'   `length(y)`.
#' @param scale `"parameter"` or `"link"`.
#' @param method `NULL` (the default) or an object inheriting from
#'   [ExpectMethod], as for [distrib_d4expected_hessian()]. See
#'   [expected_derivative_methods].
#' @param ... Passed to methods; the default method accepts `accuracy`, the
#'   order of accuracy of the stencil passed to [numericals7::fd_tensor()]
#'   (4 by default).
#'
#' @return A named list of numeric vectors, keyed as
#'   [`d5expected_names(distrib@params)`][d5expected_names].
#'
#' @examples
#' d <- poisson_distrib()
#' str(distrib_d5expected_hessian(d, 2, list(mu = 3), scale = "link"))
#'
#' g <- gaussian1_distrib()
#' r <- distrib_d5expected_hessian(g, 0, list(mu = 0, sigma = 2))
#' r[["sigma_sigma_sigma_sigma_sigma_sigma_sigma"]]
#'
#' @seealso [distrib_d4expected_hessian()] for the order below,
#'   [d5expected_names()], [distrib_deriv6()]
#'
#' @export
distrib_d5expected_hessian <- S7::new_generic(
  "distrib_d5expected_hessian", "distrib",
  function(distrib, y, theta, scale = c("parameter", "link"),
           method = NULL,
           ...) {
    args <- check_derivative_args(distrib, y, theta)
    check_method(method)
    y <- args$y
    theta <- args$theta
    scale <- match.arg(scale)
    if (!is.null(method)) {
      return(dexpected_by_method(distrib, y, theta, scale, 5L, method))
    }
    res <- S7::S7_dispatch()
    if (scale == "link") d5expected_link(distrib, y, theta, res) else res
  })


#' @title Default Fifth Derivative of the Expected Information
#' @name distrib_d5expected_hessian.distrib
#' @description
#' Applies one tensor stencil to the highest derivative of the expected
#' information that the family computes itself, through
#' [numerical_d5expected()], and signals an error where the family has no
#' analytic second derivative.
#' @param distrib A distribution object.
#' @param y A numeric vector of observations.
#' @param theta A named list of parameters.
#' @param scale Handled by the generic after dispatch; this method always
#'   returns the parameter scale.
#' @param method `NULL`: the generic handles an object before dispatching on
#'   the family. See [expected_derivative_methods].
#' @param ... Passed to [numerical_d5expected()], which accepts `accuracy`.
#' @return A named list keyed as [d5expected_names()].
#' @keywords internal
S7::method(distrib_d5expected_hessian, distrib) <- function(
    distrib, y, theta, scale = c("parameter", "link"),
    method = NULL, ...) {
  numerical_d5expected(distrib, y, theta, "parameter", ...)
}


#' Carry the Fifth Derivative of the Expected Information to the Link Scale
#'
#' @description
#' Composes the parameter-scale fifth derivative with the expected information
#' and its derivatives of orders one to four, through
#' [dexpected_link_higher()] at order five.
#'
#' @param distrib A distribution object.
#' @param y,theta The observations and the parameters, aligned.
#' @param res The fifth derivative on the parameter scale.
#'
#' @return A named list keyed as [d5expected_names()], on the link scale.
#'
#' @keywords internal
d5expected_link <- function(distrib, y, theta, res) {
  nat <- list(
    distrib_expected_hessian(distrib, y, theta, scale = "parameter"),
    distrib_dexpected_hessian(distrib, y, theta, scale = "parameter"),
    distrib_d2expected_hessian(distrib, y, theta, scale = "parameter"),
    distrib_d3expected_hessian(distrib, y, theta, scale = "parameter"),
    distrib_d4expected_hessian(distrib, y, theta, scale = "parameter"),
    res)
  h <- inverse_link_derivs(distrib, theta, 6L)
  dexpected_link_higher(distrib@params, nat, h, 5L)
}


#' The Names of the Expected Information's Fifth Derivative
#'
#' @description
#' Returns one key per pair \eqn{(a,b)}, spelled as [hess_names()] spells it,
#' and per unordered quintuple, spelled as [deriv_names()] spells it at order
#' five, joined `ab` first.
#'
#' @param params A character vector of parameter names, in the family's order.
#'
#' @return A character vector, `length(hess_names(params)) *
#'   length(deriv_names(params, 5))` long.
#'
#' @examples
#' d5expected_names(c("mu", "sigma"))
#'
#' @seealso [d5expected_key()], [d4expected_names()]
#'
#' @export
d5expected_names <- function(params) {
  as.vector(t(outer(hess_names(params), deriv_names(params, 5L), paste,
                    sep = "_")))
}


#' The Key of One Component of the Expected Information's Fifth Derivative
#'
#' @description
#' Returns the name under which [distrib_d5expected_hessian()] returns
#' \eqn{\partial^5\,\mathbb{E}[\ell_{ab}]/\partial\theta_c\,\partial\theta_d\,
#' \partial\theta_e\,\partial\theta_f\,\partial\theta_g}.
#'
#' @param params A character vector of parameter names, in the family's order.
#' @param a,b Indices of the information's pair; their order does not matter.
#' @param c,d,e,f,g Indices of the parameters with respect to which the
#'   derivative is taken; their order does not matter either.
#'
#' @return A single string.
#'
#' @examples
#' d5expected_key(c("mu", "sigma"), 1, 2, 2, 2, 1, 2, 1)
#'
#' @seealso [d5expected_names()]
#'
#' @export
d5expected_key <- function(params, a, b, c, d, e, f, g) {
  paste0(hess_pair_name(params, a, b), "_",
         paste(params[sort(c(c, d, e, f, g))], collapse = "_"))
}


#' The Order from Which the Fifth Derivative of the Information Is Differenced
#'
#' @description
#' Returns the highest order, 4, 3 or 2, at which the family registers its own
#' derivative of the expected information, and `NA` where it registers none of
#' them. For a wrapper (a family with a `parent_distrib` property) the result
#' is at most the one of the parent, since the wrapper's derivatives of each
#' order are built from the parent's of the same and lower orders.
#'
#' @param x A distribution object.
#'
#' @return An integer, 2 to 4, or `NA_integer_`.
#'
#' @seealso [numerical_d5expected()], which differences that order.
#' @keywords internal
d5expected_base <- function(x) {
  own <- if (owns_method(x, distrib_d4expected_hessian)) 4L else
    if (owns_method(x, distrib_d3expected_hessian)) 3L else
      if (has_d2expected_hessian(x)) 2L else NA_integer_
  if (!is.na(own) && S7::prop_exists(x, "parent_distrib")) {
    own <- min(own, d5expected_base(x@parent_distrib))
  }
  own
}


#' The Fifth Derivative of the Expected Information by One Stencil
#'
#' @description
#' Computes the default route behind [distrib_d5expected_hessian()] by one
#' tensor-product central stencil of [numericals7::fd_tensor()] on the highest
#' derivative of the expected information that the family computes itself, as
#' [d5expected_base()] gives it.
#'
#' @details
#' For each pair \eqn{(a,b)}, [numericals7::fd_tensor()] is called with
#' \eqn{\mathbb{E}[\ell_{ab}]} as the function and its derivatives of the base
#' order as the known quantity, and returns every quintuple. The derivatives of
#' the base order are computed once per node of the stencil and shared by the
#' pairs. With base 4 and accuracy 4 the cost is \eqn{4p} evaluations of
#' [distrib_d4expected_hessian()], for \eqn{p} parameters.
#'
#' @param distrib A distribution object.
#' @param y A numeric vector of observations.
#' @param theta A named list of parameters.
#' @param scale `"parameter"` or `"link"`. [distrib_d5expected_hessian()]
#'   uses the parameter scale and applies the link scale itself; on `"link"`
#'   the stencil is taken in the linear predictors on the link-scale
#'   derivatives, a second route kept for comparison.
#' @param accuracy The order of accuracy of the stencils, passed to
#'   [numericals7::fd_tensor()].
#'
#' @return A named list keyed as [d5expected_names()].
#'
#' @keywords internal
numerical_d5expected <- function(distrib, y, theta, scale, accuracy = 4L) {
  base <- d5expected_base(distrib)
  if (is.na(base)) {
    stop(sprintf(paste0(
      "'%s' has no analytic second derivative of its expected information,\n",
      "  and its derivative of order 5 is one stencil on an analytic\n",
      "  derivative of order 2 to 4."), distrib@distrib_name), call. = FALSE)
  }
  params <- distrib@params
  p <- length(params)
  theta <- align_theta(distrib, theta)
  cs <- tensor_coords(distrib, theta, scale)
  gen <- switch(as.character(base),
                "2" = distrib_d2expected_hessian,
                "3" = distrib_d3expected_hessian,
                "4" = distrib_d4expected_hessian)
  # the inner tuples of the base order in the order of
  # numericals7::tuple_indices(), which is hess_names() at order two and the
  # lexicographic deriv_names() above
  inner <- if (base == 2L) hess_names(params) else deriv_names(params, base)
  # the base order at each node of the stencil, computed once and read by
  # every pair
  memo <- list()
  at <- function(x) {
    for (m in memo) if (identical(m$x, x)) return(m$v)
    v <- gen(distrib, y, cs$to_theta(x), scale = scale)
    memo[[length(memo) + 1L]] <<- list(x = x, v = v)
    v
  }
  comp <- t(vapply(deriv_indices(params, 5L), tabulate, integer(p),
                   nbins = p))
  if (p == 1L) comp <- matrix(comp, ncol = 1L)
  pairs <- which(upper.tri(diag(p), diag = TRUE), arr.ind = TRUE)
  pairs <- pairs[order(pairs[, 1L], pairs[, 2L]), , drop = FALSE]
  out <- list()
  for (r in seq_len(nrow(pairs))) {
    a <- pairs[r, 1L]; b <- pairs[r, 2L]
    pair <- hess_pair_name(params, a, b)
    keys <- paste0(pair, "_", inner)
    f <- function(x) {
      distrib_expected_hessian(distrib, y, cs$to_theta(x), scale = scale)[[pair]]
    }
    known <- numericals7::fd_known(function(x) unname(at(x)[keys]),
                                   order = base)
    got <- numericals7::fd_tensor(f, cs$x, components = comp, known = known,
                                  bounds = cs$bounds, accuracy = accuracy)
    names(got) <- paste0(pair, "_", deriv_names(params, 5L))
    out[names(got)] <- unname(got)
  }
  nm <- d5expected_names(params)
  out <- out[nm]
  # a component is as long as the observations, also where the family's
  # quantities are constant in the response
  n <- max(length(y), lengths(theta[seq_len(p)]))
  lapply(stats::setNames(out, nm), function(v) rep_len(as.numeric(v), n))
}
