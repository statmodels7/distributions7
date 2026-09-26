#' @include generics.R dexpected_hessian.R
NULL

#' The Expected Information of a Location-Scale Family
#'
#' @description
#' Computes the expected information of a family whose first two parameters are
#' a location \eqn{\mu} and a scale \eqn{\sigma}, and optionally its first and
#' second derivatives in the parameters, by one quadrature per distinct value
#' of the remaining shape parameters.
#'
#' @details
#' When the log-density is \eqn{\ell(y;\mu,\sigma,s) = -\log\sigma + g(z;s)}
#' with \eqn{z = (y-\mu)/\sigma}, every derivative of \eqn{\ell} carries one
#' factor \eqn{1/\sigma} per index on \eqn{\mu} or \eqn{\sigma} and is otherwise
#' a function of \eqn{z} and \eqn{s} alone. The same holds after taking the
#' expectation, so for any component
#' \deqn{\partial^{I}\,\mathbb{E}[\ell_{ab}](\mu,\sigma,s)
#'   = \sigma^{-k}\,\partial^{I}\,\mathbb{E}[\ell_{ab}](0,1,s),}
#' with \eqn{k} the number of indices among \eqn{a}, \eqn{b} and \eqn{I} that
#' fall on \eqn{\mu} or \eqn{\sigma}. The quantity at \eqn{(0,1,s)} does not
#' depend on the location or on the scale, so it is computed once for each
#' distinct shape and not once for each observation.
#'
#' At \eqn{(0,1,s)} each quantity is an integral over \eqn{z} of the family's
#' own analytic derivatives against its density \eqn{f}. Differentiating under
#' the integral moves the measure as well as the integrand:
#' \deqn{\mathbb{E}[\ell_{ab}] = \int \ell_{ab}\,f,}
#' \deqn{\partial_c\,\mathbb{E}[\ell_{ab}] = \int (\ell_{abc}
#'   + \ell_{ab}\ell_c)\,f,}
#' \deqn{\partial_{cd}\,\mathbb{E}[\ell_{ab}] = \int (\ell_{abcd}
#'   + \ell_{abc}\ell_d + \ell_{abd}\ell_c + \ell_{ab}\ell_{cd}
#'   + \ell_{ab}\ell_c\ell_d)\,f.}
#'
#' The integrals are taken by an exp-sinh rule on each half-line, from
#' [loc_scale_rule()]. Its nodes cluster double-exponentially at \eqn{z = 0},
#' which is where the three families that use it put their one sharp feature
#' (the step of the skewing function at a large shape, and the core of the
#' pseudo-Huber at a small one), and the same transformation turns an algebraic
#' tail into an exponential one, which is what the Student t's heavy tail
#' needs.
#'
#' A node where the density is exactly zero contributes nothing, whatever the
#' derivatives read there. Any other non-finite contribution makes the whole
#' shape row `NA` rather than a sum over the finite part.
#'
#' @param distrib A family whose first two parameters are a location and a
#'   scale, and which carries analytic derivatives to the order needed.
#' @param theta The aligned parameter list, each component of the common
#'   length `n`.
#' @param order `0L` for the expected information, `1L` for its first
#'   derivatives, `2L` for its second. At order 2 the first derivatives are
#'   returned beside the second, from the same evaluation.
#' @param n The number of observations.
#' @param threads The thread count passed to the family's kernels.
#'
#' @return A named list: at order 0 keyed as [hess_names()], at order 1 as
#'   [dexpected_names()], at order 2 as [dexpected_names()] and
#'   [d2expected_names()] together. Each component has length `n` and is on
#'   the parameter scale.
#'
#' @seealso [loc_scale_rule()] for the quadrature,
#'   [distrib_expected_hessian()], [distrib_dexpected_hessian()] and
#'   [distrib_d2expected_hessian()] for the generics this serves.
#'
#' @keywords internal
loc_scale_expected <- function(distrib, theta, order, n, threads = 1L) {
  params <- distrib@params
  p <- length(params)
  rule <- loc_scale_rule()
  N <- length(rule$x)

  # the distinct shapes, and which one each observation carries
  shp <- params[-(1:2)]
  S <- vapply(shp, function(nm) rep_len(theta[[nm]], n), numeric(n))
  if (!is.matrix(S)) S <- matrix(S, nrow = n)
  key <- if (length(shp)) do.call(paste, c(as.data.frame(S), sep = "\r")) else rep("", n)
  first <- !duplicated(key)
  U <- S[first, , drop = FALSE]
  idx <- match(key, key[first])
  m <- nrow(U)

  # every node against every distinct shape, at location 0 and scale 1
  yy <- rep(rule$x, m)
  th <- c(list(rep(0, N * m), rep(1, N * m)),
          lapply(seq_along(shp), function(j) rep(U[, j], each = N)))
  names(th) <- params
  fw <- distrib_pdf(distrib, yy, th) * rep(rule$w, m)
  live <- fw != 0

  H <- distrib_hessian(distrib, yy, th, threads = threads)
  if (order >= 1L) {
    g <- distrib_gradient(distrib, yy, th, threads = threads)
    T3 <- distrib_deriv3(distrib, yy, th, threads = threads)
  }
  if (order >= 2L) T4 <- distrib_deriv4(distrib, yy, th, threads = threads)

  nm_k <- function(v) paste(params[sort(v)], collapse = "_")
  integrate_rows <- function(v) {
    v[!live] <- 0
    M <- matrix(v * fw, N, m)
    out <- colSums(M)
    out[!is.finite(out)] <- NA_real_
    out
  }
  # sigma^-k with k the count of indices on the location or the scale
  sig <- rep_len(theta[[params[2L]]], n)
  scaled <- function(val, idxs) {
    k <- sum(idxs <= 2L)
    val[idx] / sig^k
  }

  pairs <- which(upper.tri(diag(p), diag = TRUE), arr.ind = TRUE)
  out <- list()
  if (order == 0L) {
    for (r in seq_len(nrow(pairs))) {
      a <- pairs[r, 1L]; b <- pairs[r, 2L]
      out[[hess_pair_name(params, a, b)]] <-
        scaled(integrate_rows(H[[hess_pair_name(params, a, b)]]), c(a, b))
    }
    return(out[hess_names(params)])
  }
  for (r in seq_len(nrow(pairs))) {
    a <- pairs[r, 1L]; b <- pairs[r, 2L]
    Hab <- H[[hess_pair_name(params, a, b)]]
    for (c in seq_len(p)) {
      v <- T3[[nm_k(c(a, b, c))]] + Hab * g[[c]]
      out[[dexpected_key(params, a, b, c)]] <- scaled(integrate_rows(v), c(a, b, c))
    }
    if (order == 2L) {
      for (s in seq_len(nrow(pairs))) {
        c <- pairs[s, 1L]; d <- pairs[s, 2L]
        v <- T4[[nm_k(c(a, b, c, d))]] +
          T3[[nm_k(c(a, b, c))]] * g[[d]] + T3[[nm_k(c(a, b, d))]] * g[[c]] +
          Hab * (H[[hess_pair_name(params, c, d)]] + g[[c]] * g[[d]])
        out[[d2expected_key(params, a, b, c, d)]] <-
          scaled(integrate_rows(v), c(a, b, c, d))
      }
    }
  }
  out
}


#' The Quadrature Rule of the Location-Scale Expectations
#'
#' @description
#' Nodes and weights of an exp-sinh rule on each half-line, for integrals over
#' the whole real line of a density at location 0 and scale 1.
#'
#' @details
#' On \eqn{(0,\infty)} the substitution \eqn{z = \exp(\tfrac{\pi}{2}\sinh t)}
#' turns the integral into one over the whole \eqn{t} line whose integrand
#' decays double-exponentially at both ends, which the trapezoidal rule then
#' integrates with an error falling like \eqn{\exp(-c/h)}; the negative half is
#' its mirror image. The step is \eqn{h = 1/32} and the nodes run from
#' \eqn{10^{-16}} to \eqn{10^{60}} in \eqn{|z|}, about 580 in all. Against the
#' same rule at \eqn{h = 1/64}, it reproduces the expected information to
#' \eqn{3\times10^{-14}} and its second derivative to \eqn{7\times10^{-9}} on the
#' skew normal up to \eqn{|\alpha| = 500}, the hardest case measured; at
#' \eqn{h = 1/24} that second derivative is \eqn{2\times10^{-4}} out. Past the
#' upper end a tail decaying as \eqn{|z|^{-1-\nu}} leaves a relative
#' \eqn{10^{-60\nu}/\nu}, negligible for \eqn{\nu \ge 0.3}. The end is not taken
#' further because the Student t's derivatives in \eqn{\nu} are not finite at
#' \eqn{|z| = 10^{130}}, where its density is still positive.
#'
#' @return A list with the nodes `x` and the weights `w`.
#'
#' @seealso [loc_scale_expected()]
#'
#' @keywords internal
loc_scale_rule <- local({
  cache <- NULL
  function() {
    if (is.null(cache)) {
      c0 <- pi / 2
      h <- 1 / 32
      tlo <- -asinh(16 * log(10) / c0)
      thi <- asinh(60 * log(10) / c0)
      t <- seq(ceiling(tlo / h), floor(thi / h)) * h
      x <- exp(c0 * sinh(t))
      w <- h * c0 * cosh(t) * x
      cache <<- list(x = c(-rev(x), x), w = c(rev(w), w))
    }
    cache
  }
})


#' Is a Family's Exact Expected Information Costly to Compute?
#'
#' @description
#' `TRUE` when the family's expected information is exact but costs far more
#' than its observed information, because it is obtained by an integral or a
#' sum for each observation or each distinct shape; `FALSE` otherwise.
#'
#' @details
#' The answer is a statement about cost, and it is read beside
#' [expected_hessian_exact()], which states accuracy. Two kinds of family
#' answer `TRUE`.
#'
#' The location-scale families [skewnormal1_distrib()],
#' [skewnormal2_distrib()], [skewt_distrib()] and [pseudohuber_distrib()]
#' take one integral over the standardized response per distinct shape, from
#' [loc_scale_expected()]: nothing where the shape is the same for every
#' observation, one integral per observation where it is modelled. Measured
#' on smooth regressions at 1000 observations with the shape developed over a
#' covariate, a fit that takes the scoring step on this information costs 7
#' to 21 times one that takes it on the observed information, at the same
#' estimate.
#'
#' The Poisson-inverse Gaussians [pig1_distrib()] and [pig2_distrib()] sum
#' over the support for each observation, over a number of terms that grows
#' as \eqn{\sigma\mu}: about 25 microseconds an observation at
#' \eqn{\mu = 3}, \eqn{\sigma = 0.3} and 2 milliseconds at \eqn{\mu = 30},
#' \eqn{\sigma = 3}, against a few microseconds for the observed information.
#' Measured at 1000 observations, a fit on it costs 1 to 6 times one on the
#' observed information at the same estimate, and far more where the fit
#' passes through a large \eqn{\sigma\mu}.
#'
#' \pkg{statmodels7}'s `iwls(hessian = "auto")` therefore settles on the
#' observed information for these families. The default answers for a
#' wrapper's parent, read off its `parent_distrib` property, and `FALSE` for a
#' family without one.
#'
#' @param x An object inheriting from class `"distrib"`.
#' @param ... Passed to methods.
#'
#' @return A single logical.
#'
#' @examples
#' expected_hessian_costly(skewt_distrib())
#' expected_hessian_costly(pig1_distrib())
#' expected_hessian_costly(gaussian1_distrib())
#' expected_hessian_costly(fixed(skewt_distrib(), nu = 6))
#'
#' @seealso [expected_hessian_exact()], [loc_scale_expected()]
#'
#' @export
expected_hessian_costly <- S7::new_generic(
  "expected_hessian_costly", "x",
  function(x, ...) S7::S7_dispatch())

S7::method(expected_hessian_costly, distrib) <- function(x, ...) {
  if ("parent_distrib" %in% S7::prop_names(x)) {
    return(expected_hessian_costly(x@parent_distrib))
  }
  FALSE
}


#' Register the Location-Scale Expected Information on a Family
#'
#' @description
#' Registers [distrib_expected_hessian()], [distrib_dexpected_hessian()] and
#' [distrib_d2expected_hessian()] on a location-scale class, each computed by
#' [loc_scale_expected()], and [expected_hessian_costly()] answering
#' `TRUE`.
#'
#' @details
#' The `approx` and `nsim` arguments of the three generics are accepted and
#' ignored, the expectation being taken by the rule of [loc_scale_rule()] and
#' not by one of the approximations `approx` selects between.
#'
#' @param cls An S7 class whose first two parameters are a location and a
#'   scale.
#'
#' @return `NULL`, invisibly; called for its side effect.
#'
#' @keywords internal
register_loc_scale_expected <- function(cls) {
  S7::method(distrib_expected_hessian, cls) <- function(
      distrib, y, theta, scale = c("parameter", "link"),
      approx = c("opg", "bartlett", "integrate", "mc"), nsim = 10000, ...,
      threads = 1L) {
    loc_scale_expected(distrib, theta, 0L, length(y), threads)
  }
  S7::method(distrib_dexpected_hessian, cls) <- function(
      distrib, y, theta, scale = c("parameter", "link"),
      approx = c("opg", "bartlett", "integrate", "mc"), nsim = 10000, ...,
      threads = 1L) {
    dexpected_analytic(distrib, y, theta, match.arg(scale), 1L, threads,
                       function(k) loc_scale_expected(distrib, theta, k,
                                                      length(y), threads))
  }
  S7::method(distrib_d2expected_hessian, cls) <- function(
      distrib, y, theta, scale = c("parameter", "link"),
      approx = c("opg", "bartlett", "integrate", "mc"), nsim = 10000, ...,
      threads = 1L) {
    dexpected_analytic(distrib, y, theta, match.arg(scale), 2L, threads,
                       function(k) loc_scale_expected(distrib, theta, k,
                                                      length(y), threads))
  }
  S7::method(expected_hessian_costly, cls) <- function(x, ...) TRUE
  invisible(NULL)
}
