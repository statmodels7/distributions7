#' @include partition_sums.R
NULL

# The derivatives of the logarithm of a tail probability far in that tail.
#
# Composed from the derivatives of F and divided by the tail probability P,
# they lose digits as P falls (the ratios grow like a power of the distance
# and cancel in the moment-to-cumulant sum) and are NaN once P underflows,
# although log P is finite. Two routes take them directly:
#
#   ls_log_tail_derivs()    the location-scale families, log P = g(z) with z
#                           linear in the location, by Faa di Bruno over the
#                           partials of z with g^(k) written out per base
#   quad_log_tail_derivs()  every other family, log P = l(q) + log J with J
#                           the integral of exp(l(y) - l(q)) over the tail,
#                           its partials by the moment-to-cumulant sum over the
#                           Bell polynomials in the differences l(y) - l(q)
#
# log_tail_route() is called from the bodies of the four cdf-derivative
# generics and returns NULL where neither applies.

#' The Log-Tail Route of the CDF Derivatives
#'
#' @description
#' Called from the bodies of [distrib_grad_cdf()], [distrib_hess_cdf()],
#' [distrib_deriv3_cdf()] and [distrib_deriv4_cdf()] before dispatch. On the
#' log scale it returns the derivatives of \eqn{\log P}, \eqn{P} the requested
#' tail, from [ls_log_tail_derivs()] for the eight location-scale families
#' listed in `log_tail_ls`, and from [quad_log_tail_derivs()] for the
#' families listed in `log_tail_quad` (and, for the lower tail,
#' `log_tail_quad_lower`) at the quantiles where \eqn{\log P}
#' is below `log_tail_cut`, the other quantiles going to the family's own
#' method. Elsewhere it returns `NULL` and the generic dispatches as usual.
#'
#' `log_tail_ls` maps each location-scale class to its base and to the map
#' \eqn{w} of its second parameter; `log_tail_quad` lists the classes of the
#' quadrature route for both tails and `log_tail_quad_lower` those for the
#' lower tail, which adds the two Weibull families; `log_tail_cut` is the
#' threshold, \eqn{-5}; `log_tail_off` holds the switch that turns the
#' route off while quantiles are sent back to the family's method.
#'
#' The quantiles are selected on [tail_cdf()], the value the family's method
#' computes first. While the method runs on the quantiles that are not far,
#' `log_tail_stash` holds that value with the class, quantiles, parameters
#' and tail it belongs to, and [tail_cdf()] returns it for a call with the
#' same arguments instead of evaluating the distribution function again. A
#' parameter of length one is passed to the method as a scalar.
#'
#' @param distrib An object inheriting from `distrib`.
#' @param q A numeric vector of quantiles, already recycled by the generic.
#' @param theta A named list of parameters, already aligned by the generic.
#' @param order The derivative order, 1 to 4.
#' @param lower.tail Is the lower tail wanted? A single logical.
#' @param log Are derivatives of the log probability wanted? A single logical.
#' @param generic The generic of that order, called back on the quantiles
#'   that stay with the family's method: those that are not far, and those
#'   whose tail reaches outside the double range, which
#'   [quad_log_tail_derivs()] marks.
#'
#' @return A named list of numeric vectors, keyed as
#'   [`deriv_names(distrib@params, order)`][deriv_names] (in the order of
#'   [hess_names()] at order two), or `NULL`.
#'
#' @examples
#' d <- gaussian1_distrib()
#' # log S = log(1 - Phi(z)) at z = 40, where 1 - Phi(z) underflows
#' distributions7:::log_tail_route(d, 40, list(mu = 0, sigma = 1), 1L,
#'                                 FALSE, TRUE, distrib_grad_cdf)
#'
#' @keywords internal
log_tail_route <- function(distrib, q, theta, order, lower.tail, log, generic) {
  if (!length(q) || log_tail_off$on) return(NULL)
  cls <- attr(S7::S7_class(distrib), "name")
  if (is.null(cls)) return(NULL)
  spec <- log_tail_ls[[cls]]
  if (!is.null(spec)) {
    if (!isTRUE(log)) return(NULL)
    return(ls_log_tail_derivs(distrib, q, theta, order, spec, lower.tail))
  }
  if (!(cls %in% (if (lower.tail) log_tail_quad_lower else log_tail_quad)))
    return(NULL)
  n <- max(length(q), lengths(theta[distrib@params]))
  q <- rep_len(q, n)
  th <- theta[distrib@params]
  # the selection reads the value the family's method computes first, and
  # hands it to that method through log_tail_stash
  Fq <- tail_cdf(distrib, q, th, lower.tail)
  lp <- base::log(tail_probability(Fq, lower.tail))
  far <- !is.na(lp) & lp < log_tail_cut & is.finite(q) &
    q > distrib@bounds[1L] & q < distrib@bounds[2L]
  nm <- deriv_names(distrib@params, order)
  out <- stats::setNames(lapply(nm, function(i) numeric(n)), nm)
  # a parameter of length one stays a scalar, so that a method with a path
  # for a scalar parameter keeps it
  sub <- function(idx) lapply(th, function(v) if (length(v) == 1L) v else v[idx])
  if (any(!far)) {
    qn <- q[!far]
    thn <- sub(!far)
    Fn <- as.numeric(Fq)[!far]
    if (!lower.tail) attr(Fn, "upper") <- attr(Fq, "upper")[!far]
    list2env(list(cls = cls, q = qn, theta = thn, lower = lower.tail, Fq = Fn),
             envir = log_tail_stash)
    # the call back to the generic does not come through the route again
    log_tail_off$on <- TRUE
    on.exit({
      log_tail_off$on <- FALSE
      rm(list = ls(log_tail_stash), envir = log_tail_stash)
    }, add = TRUE)
    near <- generic(distrib, qn, thn, lower.tail = lower.tail, log = log)
    log_tail_off$on <- FALSE
    rm(list = ls(log_tail_stash), envir = log_tail_stash)
    if (!any(far)) return(near)
    for (i in nm) out[[i]][!far] <- near[[i]]
  }
  if (isTRUE(log)) {
    tail <- quad_log_tail_derivs(distrib, q[far], sub(far), order, lower.tail)
    lost <- attr(tail, "lost")
  } else {
    # d^I P = P B_I, the complete Bell polynomial in the partials of log P
    lt <- quad_log_tail_derivs(distrib, q[far], sub(far), order, lower.tail,
                               all_orders = TRUE)
    ell <- function(b) lt[[length(b)]][[canon_key(b, distrib@params)]]
    P <- exp(distrib_cdf(distrib, q[far], sub(far), lower.tail = lower.tail,
                         log.p = TRUE))
    tail <- stats::setNames(lapply(deriv_indices(distrib@params, order),
                                   function(I) P * bell_f_ratio(distrib@params[I], ell)),
                            nm)
    lost <- attr(lt, "lost")
  }
  for (i in nm) out[[i]][far] <- tail[[i]]
  # a tail with mass outside the double range (a shape near zero, whose
  # density reaches below the smallest double) stays with the family's
  # method, which integrates that part in closed form
  if (any(lost)) {
    back <- which(far)[lost]
    log_tail_off$on <- TRUE
    on.exit(log_tail_off$on <- FALSE, add = TRUE)
    redo <- generic(distrib, q[back], sub(back), lower.tail = lower.tail,
                    log = log)
    log_tail_off$on <- FALSE
    for (i in nm) out[[i]][back] <- redo[[i]]
  }
  if (order == 2L) out[hess_names(distrib@params)] else out
}

#' @rdname log_tail_route
#' @usage NULL
#' @format NULL
log_tail_cut <- -5

#' @rdname log_tail_route
#' @usage NULL
#' @format NULL
log_tail_off <- new.env()
log_tail_off$on <- FALSE

#' @rdname log_tail_route
#' @usage NULL
#' @format NULL
log_tail_stash <- new.env()

#' @rdname log_tail_route
#' @usage NULL
#' @format NULL
log_tail_quad <- c(
  "Beta1Distrib", "Beta2Distrib", "ChisqDistrib", "EnetDistrib",
  "Gamma1Distrib", "Gamma2Distrib", "GenGamma1Distrib", "GenGamma2Distrib",
  "InvGauss1Distrib", "InvGauss2Distrib", "Lognormal2Distrib",
  "PseudoHuberDistrib", "PseudoHuber2Distrib", "SkewNormal1Distrib",
  "SkewNormal2Distrib", "SkewTDistrib", "StudentT1Distrib",
  "StudentT2Distrib", "VonMises1Distrib", "VonMises2Distrib")
#' @rdname log_tail_route
#' @usage NULL
#' @format NULL
log_tail_quad_lower <- c(log_tail_quad, "Weibull1Distrib", "Weibull3Distrib")

#' @rdname log_tail_route
#' @usage NULL
#' @format NULL
log_tail_ls <- list(
  Gaussian1Distrib = list(base = "gaussian", w = "inv"),
  Gaussian2Distrib = list(base = "gaussian", w = "invsqrt"),
  Gaussian3Distrib = list(base = "gaussian", w = "sqrt"),
  Lognormal1Distrib = list(base = "gaussian", w = "invsqrt", logq = TRUE),
  LaplaceDistrib = list(base = "laplace", w = "inv"),
  Laplace2Distrib = list(base = "laplace", w = "rate"),
  LogisticDistrib = list(base = "logistic", w = "inv"),
  GumbelDistrib = list(base = "gumbel", w = "inv"))

#' Log-Tail Derivatives of the Location-Scale Families
#'
#' @description
#' Computes, in compiled code (`ls_log_tail_cpp()`), the derivatives of
#' order `o` of \eqn{\log P} in
#' \eqn{(\mu, \theta_2)} for a family whose tail probability is
#' \eqn{P = G(z)}, \eqn{z = (x - \mu)\,w(\theta_2)}, with \eqn{x = q} or
#' \eqn{\log q}. With \eqn{g = \log G}, every partial derivative is the sum
#' over the set partitions \eqn{\pi} of the multi-index of
#' \eqn{g^{(|\pi|)}(z)\prod_{B \in \pi}\partial^B z}; \eqn{z} is linear in
#' \eqn{\mu}, so a block naming \eqn{\mu} twice is zero.
#'
#' @details
#' # The bases
#'
#' The lower tail of a symmetric base is its upper tail at \eqn{-z}. The
#' derivatives \eqn{g^{(k)}}, \eqn{k \le 4}, are written out:
#' \itemize{
#'   \item Laplace: \eqn{g = \log\tfrac12 - z} for \eqn{z > 0} and
#'     \eqn{\log(1 - u)}, \eqn{u = e^{z}/2}, otherwise, whose derivatives are
#'     rational in \eqn{u}.
#'   \item logistic: polynomials in \eqn{p = F} and \eqn{s = 1 - F}, each
#'     taken from its own tail by [stats::plogis()].
#'   \item Gumbel, upper tail: with \eqn{w = e^{-z}}, \eqn{g' = -\beta(w)},
#'     \eqn{\beta(w) = w/(e^{w} - 1)}, and each further order applies
#'     \eqn{-w\,d/dw}; below \eqn{w = 1/2} by the Bernoulli series
#'     \eqn{\sum_n n^k B_n w^n/n!}, above it in closed form. Lower tail:
#'     \eqn{\log F = -e^{-z}}.
#'   \item Gaussian: below \eqn{z = 9} from the hazard
#'     \eqn{m = \phi/(1 - \Phi)} and \eqn{m' = m(m - z)}; from \eqn{z = 9}
#'     by the asymptotic expansion
#'     \eqn{g = -z^2/2 - \log z - \tfrac12\log 2\pi + \log\sum_{n \le 40}
#'     (-1)^n (2n-1)!!\,z^{-2n}}, whose truncation error there is below
#'     \eqn{e^{-40}}. For \eqn{z > 0} the term \eqn{-z^2/2} is
#'     differentiated exactly as \eqn{-(x - \mu)^2 w^2/2} and only the rest
#'     goes through the partition sum, which keeps the cancellation between
#'     the two out of the result.
#' }
#'
#' Outside the support of a lognormal (\eqn{q \le 0}) the derivatives of
#' \eqn{\log S} are zero and those of \eqn{\log F} are `NaN`.
#'
#' @param distrib One of the families of `log_tail_ls`.
#' @param q A numeric vector of quantiles.
#' @param theta A named list with the location and the second parameter.
#' @param o The derivative order, 1 to 4.
#' @param spec The entry of `log_tail_ls` for the family.
#' @param lower.tail Is the lower tail wanted? A single logical.
#'
#' @return A named list, one component per multi-index of order `o`, in the
#'   order of [hess_names()] at order two.
#'
#' @examples
#' d <- gaussian2_distrib()
#' # log S at 400 standard deviations, where the tail underflows
#' distributions7:::ls_log_tail_derivs(d, 400, list(mu = 0, sigma2 = 1), 2L,
#'   list(base = "gaussian", w = "invsqrt"), lower.tail = FALSE)
#'
#' @keywords internal
ls_log_tail_derivs <- function(distrib, q, theta, o, spec, lower.tail) {
  params <- distrib@params
  n <- max(length(q), length(theta[[1]]), length(theta[[2]]))
  x <- rep_len(as.numeric(q), n)
  inside <- rep_len(TRUE, n)
  if (isTRUE(spec$logq)) {
    inside <- x > 0
    x <- log(ifelse(inside, x, 1))
  }
  base <- if (lower.tail && spec$base == "gumbel") 4L else
    match(spec$base, c("gaussian", "laplace", "logistic", "gumbel")) - 1L
  m <- ls_log_tail_cpp(x, as.numeric(theta[[1]]), as.numeric(theta[[2]]), base,
                       match(spec$w, c("inv", "rate", "invsqrt", "sqrt")) - 1L,
                       lower.tail && spec$base != "gumbel",
                       spec$base == "gaussian", o)
  out <- lapply(seq_len(o + 1L), function(c) {
    v <- m[, c]
    if (isTRUE(spec$logq)) v <- ifelse(inside, v, if (lower.tail) NaN else 0)
    v
  })
  out <- stats::setNames(out, deriv_names(params, o))
  if (o == 2L) out[hess_names(params)] else out
}



#' Log-Tail Derivatives by Quadrature of the Log-Density Derivatives
#'
#' @description
#' Computes the derivatives of order `o` of \eqn{\log P}, \eqn{P} the tail
#' beyond `q` on the side that `lower.tail` selects, from
#' \deqn{\log P = \ell(q) + \log J, \qquad
#'   J = \int_{\mathrm{tail}} e^{D(y)}\,dy, \qquad D(y) = \ell(y) - \ell(q).}
#' The derivatives of \eqn{\log J} follow from the ratios
#' \eqn{\partial^B J / J}, which are weighted means of the complete Bell
#' polynomials in the partial derivatives of \eqn{D}, through the
#' moment-to-cumulant sum; both sums are the recursions over the subsets of a
#' multi-index that contain its first position, run in compiled code
#' (`log_tail_cumulants_cpp()`) on the pairs `log_tail_pairs()` enumerates. Every partial of
#' \eqn{D} is a difference of the family's own log-density derivatives at
#' \eqn{y} and at \eqn{q}, which stays small over the tail, so the sums carry
#' no cancellation of the size the ratios of [cdf_scale_k()] carry, and
#' nothing underflows where \eqn{P} does.
#'
#' @details
#' The integral is taken with the 15-point Kronrod rule on 16 panels in
#' \eqn{u \in [0, 128]}. With \eqn{s = 1/|\partial_y \ell(q)|} the scale on
#' which \eqn{\ell} falls away from \eqn{q}, the map is
#' \eqn{y = q \pm s\,(e^{u} - 1)} toward an infinite end of the support,
#' which turns an exponential tail into a double exponential and a polynomial
#' one into an exponential in \eqn{u}, and
#' \eqn{y = b + (q - b)\exp\{-c\,(e^{u} - 1)\}}, \eqn{c = \min(1, s/|q - b|)},
#' toward a finite end \eqn{b}, which reaches \eqn{b} at \eqn{u = \infty} and
#' turns a power of \eqn{y - b} into a double exponential. The route is
#' meant for the far tail, where the density decreases away from `q`;
#' [log_tail_route()] uses it below \eqn{\log P = -5}. Toward an end at zero,
#' a row whose last node with a representable abscissa still carries a
#' relative weight above \eqn{10^{-12}}, the nodes beyond it having
#' underflowed to zero, has mass below the smallest double (a shape near
#' zero); it is marked in the attribute `"lost"` of the result.
#' Toward an infinite end the nodes past the largest double are dropped,
#' which leaves out a relative mass of order \eqn{(q/x_{\max})^{\nu}} for a
#' polynomial tail of index \eqn{\nu}.
#'
#' `quad_tail_parts()` returns the nodes at every quantile, the response
#' vector \eqn{(q, y_1, \dots)} with the parameters repeated to match, the
#' log densities and the weights \eqn{e^{D}} times the Jacobian and the
#' rule's weight; `log_tail_nodes()` the nodes and weights in \eqn{u},
#' computed once; `log_tail_pairs()` the components of orders 1 to `o` in
#' canonical form and, for each, the pairs (subset containing the first
#' position, rest) as indices into the components.
#'
#' @param distrib A continuous distribution object.
#' @param q A numeric vector of quantiles inside the support.
#' @param theta A named list of parameters, each of the length of `q`.
#' @param o The derivative order, 1 to 4.
#' @param lower.tail Is the lower tail wanted? A single logical.
#'
#' @param all_orders Return the orders 1 to `o`, as a list of such lists?
#'
#' @return A named list, one component per multi-index of order `o`, in the
#'   order of [deriv_names()], with the attribute `"lost"`.
#'
#' @examples
#' d <- gamma1_distrib()
#' # log S of a gamma at 3000, where S underflows
#' distributions7:::quad_log_tail_derivs(d, 3000, list(mu = 1.7, phi = 2.3),
#'                                       1L, lower.tail = FALSE)
#'
#' @keywords internal
quad_log_tail_derivs <- function(distrib, q, theta, o, lower.tail,
                                 all_orders = FALSE) {
  params <- distrib@params
  qp <- quad_tail_parts(distrib, q, theta, lower.tail)
  n <- qp$n
  Wt <- qp$Wt / rowSums(qp$Wt)
  gens <- list(distrib_gradient, distrib_hessian, distrib_deriv3,
               distrib_deriv4)
  tabs <- lapply(seq_len(o), function(k) gens[[k]](distrib, qp$y, qp$theta))
  pt <- log_tail_pairs(params, o)
  tl <- lapply(seq_along(pt$keys),
               function(c) as.numeric(tabs[[pt$order[c]]][[pt$keys[c]]]))
  kap <- log_tail_cumulants_cpp(tl, Wt, pt$pstart, pt$ps, pt$pr)
  one <- function(k) {
    out <- lapply(deriv_indices(params, k), function(I) {
      key <- canon_key(params[I], params)
      c <- match(key, pt$keys)
      tl[[c]][seq_len(n)] + kap[, c]
    })
    stats::setNames(out, deriv_names(params, k))
  }
  out <- if (all_orders) lapply(seq_len(o), one) else one(o)
  attr(out, "lost") <- qp$lost
  out
}

#' @rdname quad_log_tail_derivs
#' @usage NULL
log_tail_pairs <- local({
  cache <- list()
  function(params, o) {
    id <- paste(c(params, o), collapse = "|")
    if (!is.null(cache[[id]])) return(cache[[id]])
    comps <- unlist(lapply(seq_len(o), function(k)
      lapply(deriv_indices(params, k), sort)), recursive = FALSE)
    keys <- vapply(comps, function(I) canon_key(params[I], params), "")
    ps <- pr <- integer()
    pstart <- 0L
    for (I in comps) {
      k <- length(I)
      for (m in seq_len(2^(k - 1)) - 1L) {
        sel <- c(TRUE, bitwAnd(m, 2^(seq_len(k - 1) - 1)) > 0)[seq_len(k)]
        ps <- c(ps, match(canon_key(params[I[sel]], params), keys) - 1L)
        pr <- c(pr, if (all(sel)) -1L else
          match(canon_key(params[I[!sel]], params), keys) - 1L)
      }
      pstart <- c(pstart, length(ps))
    }
    out <- list(keys = keys, order = lengths(comps), ps = ps, pr = pr,
                pstart = pstart)
    cache[[id]] <<- out
    out
  }
})

#' @rdname quad_log_tail_derivs
#' @return `quad_log_tail_prob()` returns \eqn{\log P = \ell(q) + \log J},
#'   a numeric vector of the length of `q`.
quad_log_tail_prob <- function(distrib, q, theta, lower.tail) {
  qp <- quad_tail_parts(distrib, q, theta, lower.tail)
  qp$lf[seq_len(qp$n)] + log(rowSums(qp$Wt))
}

#' @rdname quad_log_tail_derivs
#' @usage NULL
quad_tail_parts <- function(distrib, q, theta, lower.tail) {
  params <- distrib@params
  n <- length(q)
  nd <- log_tail_nodes()
  K <- length(nd$u)
  theta <- lapply(theta[params], rep_len, length.out = n)
  U <- matrix(nd$u, n, K, byrow = TRUE)
  b <- if (lower.tail) distrib@bounds[1L] else distrib@bounds[2L]
  sg <- if (lower.tail) -1 else 1
  # the scale on which log f falls away from q into the tail
  gy <- sg * distrib_grad_y(distrib, q, theta)
  s <- ifelse(is.finite(gy) & gy < 0, -1 / gy, Inf)
  if (is.finite(b)) {
    # toward a finite end, y = b + (q - b) exp(-c (e^u - 1)), c = min(1, s/L)
    # with L = |q - b|: a step of s at q, and the end reached at u = Inf
    L <- abs(q - b)
    cc <- pmin(1, s / L)
    E <- exp(-cc * expm1(U))
    Y <- b + (q - b) * E
    jac <- L * E * cc * exp(U)
  } else {
    s[!is.finite(s)] <- 1
    Y <- q + sg * s * expm1(U)
    jac <- s * exp(U)
  }
  ok <- is.finite(Y) & (if (lower.tail) Y < q else Y > q) &
    Y > distrib@bounds[1L] & Y < distrib@bounds[2L]
  Y[!ok] <- q[row(Y)[!ok]]
  # every node and q itself in one call, row i of the parameters repeated
  # K + 1 times
  y <- c(q, as.vector(Y))
  idx <- c(seq_len(n), rep(seq_len(n), K))
  th <- lapply(theta, `[`, idx)
  lf <- distrib_pdf(distrib, y, th, log = TRUE)
  wq <- matrix(nd$w, n, K, byrow = TRUE)
  Wt <- exp(matrix(lf[-seq_len(n)], n, K) - lf[seq_len(n)]) * jac * wq
  Wt[!ok | !is.finite(Wt)] <- 0
  # toward an end at zero, a row whose last usable node still carries weight,
  # the nodes beyond it having underflowed to the end, has mass below the
  # smallest double
  kl <- apply(Wt > 0, 1L, function(r) if (any(r)) max(which(r)) else 0L)
  wl <- Wt[cbind(seq_len(n), pmax(kl, 1L))]
  lost <- isTRUE(b == 0) & kl > 0L & kl < K & wl > 1e-12 * rowSums(Wt)
  list(n = n, K = K, y = y, theta = th, lf = lf, Wt = Wt, lost = lost)
}

#' @rdname quad_log_tail_derivs
#' @usage NULL
log_tail_nodes <- local({
  cache <- NULL
  function() {
    if (is.null(cache)) {
      r <- numericals7::gauss_kronrod15()
      br <- c(0, 0.25, 0.5, 1, 2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 64, 96, 128)
      a <- br[-length(br)]
      h <- diff(br) / 2
      cache <<- list(u = as.vector(outer(r$nodes, h) + rep(a + h, each = 15)),
                     w = as.vector(outer(r$wk, h)))
    }
    cache
  }
})

#' Refine a Quantile Far in a Tail
#'
#' @description
#' Corrects quantiles `x` of the requested tail where its log probability is
#' below `log_tail_cut`, by a safeguarded Newton iteration on
#' \eqn{r = \log P(x) - \log p}, with \eqn{d\log P/dx = \pm f/P} from the
#' family's own log distribution function and log density.
#'
#' @details
#' On a positive support, and where \eqn{x} lies on the side of the tail with
#' \eqn{|x| > 1}, the variable is \eqn{v = \log|x|}, on which a polynomial
#' tail and a power of \eqn{x} near zero are linear; elsewhere it is \eqn{x}.
#' Each step is the Newton step, limited to 10 units of \eqn{v} (or to
#' \eqn{\max(1, |x|)} in \eqn{x}), and replaced by the midpoint of the
#' current bracket of the root when it leaves it. A start that is not finite
#' or not inside the support is replaced by 1 on a positive support and by 0
#' otherwise. The iteration stops when \eqn{|r|} is within four units in the
#' last place of \eqn{\log p}; an element that does not reach it in 200
#' steps is given to the default quantile of `continuous_distrib`.
#'
#' @param distrib A continuous distribution object.
#' @param x A numeric vector of starting quantiles.
#' @param p,lower.tail,log.p As in [distrib_quantile()].
#' @param theta A named list of parameters.
#'
#' @return A numeric vector of quantiles, `x` itself where the tail is not
#'   far.
#'
#' @examples
#' d <- invgauss1_distrib()
#' th <- list(mu = 1.7, phi = 2.3)
#' x0 <- suppressWarnings(statmod::qinvgauss(-800, mean = 1.7, dispersion = 2.3,
#'                                           lower.tail = FALSE, log.p = TRUE))
#' x <- distributions7:::tail_quantile_refine(d, x0, -800, th, FALSE, TRUE)
#' distrib_cdf(d, c(x0, x), th, lower.tail = FALSE, log.p = TRUE)
#'
#' @keywords internal
tail_quantile_refine <- function(distrib, x, p, theta, lower.tail, log.p) {
  params <- distrib@params
  n <- max(length(x), length(p), lengths(theta[params]))
  x <- rep_len(as.numeric(x), n)
  lt <- rep_len(if (log.p) p else log(p), n)
  far <- which(!is.na(lt) & lt < log_tail_cut)
  if (!length(far)) return(x)
  th <- lapply(theta[params], function(v) rep_len(v, n)[far])
  b <- distrib@bounds
  pos <- b[1] >= 0
  sg <- if (lower.tail) 1 else -1
  xx <- x[far]
  tg <- lt[far]
  bad0 <- !is.finite(xx) | xx <= b[1] | xx >= b[2]
  xx[bad0] <- if (pos) 1 else 0
  # the variable: log|x| on a positive support or on the tail's side
  lg <- pos | (abs(xx) > 1 & sign(xx) == -sg)
  sx <- ifelse(xx < 0, -1, 1)
  to_x <- function(v, i) ifelse(lg[i], sx[i] * exp(v), v)
  v <- ifelse(lg, log(abs(xx)), xx)
  vlo <- rep(-Inf, length(v))
  vhi <- rep(Inf, length(v))
  done <- logical(length(v))
  for (it in seq_len(200L)) {
    act <- which(!done)
    if (!length(act)) break
    xa <- to_x(v[act], act)
    tha <- lapply(th, `[`, act)
    lP <- distrib_cdf(distrib, xa, tha, lower.tail = lower.tail, log.p = TRUE)
    r <- lP - tg[act]
    hit <- !is.na(r) & abs(r) <= 4 * .Machine$double.eps * pmax(1, abs(tg[act]))
    done[act[hit]] <- TRUE
    # d log P / dv, and the side of the root given by the sign of r
    slope <- sg * exp(distrib_pdf(distrib, xa, tha, log = TRUE) - lP) *
      ifelse(lg[act], xa, 1)
    up <- (r < 0) == (slope > 0)
    ok <- !hit & !is.na(r) & !is.na(slope) & slope != 0
    vlo[act[ok & up]] <- v[act[ok & up]]
    vhi[act[ok & !up]] <- v[act[ok & !up]]
    lim <- ifelse(lg[act], 10, pmax(1, abs(xa)))
    vn <- v[act] - pmax(-lim, pmin(lim, r / slope))
    inb <- vn > vlo[act] & vn < vhi[act]
    both <- is.finite(vlo[act]) & is.finite(vhi[act])
    vn[!inb & both] <- ((vlo[act] + vhi[act]) / 2)[!inb & both]
    # in x, a step past an end of the support is halved towards it
    xn <- to_x(vn, act)
    out_s <- !lg[act] & (xn <= b[1] | xn >= b[2])
    vn[out_s] <- ((v[act] + pmin(pmax(xn, b[1]), b[2])) / 2)[out_s]
    stall <- ok & vn == v[act]
    done[act[stall]] <- TRUE
    v[act[ok]] <- vn[ok]
    done[act[!ok & !hit]] <- TRUE
    bad0[act[!ok & !hit]] <- TRUE
  }
  xx <- to_x(v, seq_along(v))
  # whatever did not converge goes to the default root finding
  redo <- which(!done | (bad0 & !is.finite(xx)))
  if (length(redo)) {
    def <- S7::method(distrib_quantile, continuous_distrib)
    xx[redo] <- def(distrib, tg[redo], lapply(th, `[`, redo),
                    lower.tail = lower.tail, log.p = TRUE)
  }
  x[far] <- xx
  x
}
