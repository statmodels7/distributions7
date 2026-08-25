#' @include cdf_survival_higher.R
NULL

# The families whose distribution function is a sum of normal tails.
#
# Both the inverse gaussian and the elastic net have
#
#   F = const + sum_k s_k exp(w_k) Phi(x_k),
#
# so every derivative is a Leibniz split between the weight and the tail, with
# a Faa di Bruno on each side. A family supplies, per term, its sign, the
# partial derivatives of the LOG weight and those of the argument, and gets
# four orders. The weight is never formed on its own: it is combined with the
# tail on the log scale, `exp(w + pnorm(x, log.p = TRUE))`, because the
# inverse gaussian's `exp(2/(phi mu))` overflows exactly where its `Phi(b)`
# underflows.

#' Derivatives of the Standard Normal Distribution Function
#'
#' @description
#' Returns the factor \eqn{h_m} in \eqn{\Phi^{(m)}(x) = h_m(x)\,\varphi(x)} for
#' \eqn{m = 1, \ldots, 4}: \eqn{1}, \eqn{-x}, \eqn{x^2 - 1} and
#' \eqn{3x - x^3}. These are the probabilists' Hermite polynomials up to a
#' sign, and factoring the density out of them is what keeps the tails
#' evaluable: \eqn{\varphi} is supplied by the caller and multiplied in once.
#'
#' @section Notation:
#' \eqn{\Phi} is the standard normal distribution function and \eqn{\varphi}
#' its density.
#'
#' @param x A numeric vector of arguments.
#' @param m The order, 1 to 4. Any other value returns `NULL`, `switch()`
#'   falling through; no caller passes one.
#'
#' @return A numeric vector the length of `x`, the polynomial factor alone.
#'   Multiply by `dnorm(x)` for the derivative itself.
#'
#' @seealso [phi_terms_cdf_deriv_k()], the one consumer.
#'
#' @examples
#' # The first derivative of Phi is the density itself, so h1 is 1.
#' all.equal(distributions7:::phi_hermite(0.5, 1) * dnorm(0.5), dnorm(0.5))
#'
#' @keywords internal
phi_hermite <- function(x, m) {
  switch(m,
    rep_len(1, length(x)),
    -x,
    x^2 - 1,
    3 * x - x^3
  )
}

#' CDF Derivatives of a Sum of Weighted Normal Tails
#'
#' @description
#' Evaluates every component of \eqn{\partial^I F} of the requested order for a
#' distribution function of the form
#' \deqn{F = c_0 + \sum_k s_k\,e^{w_k}\,\Phi(x_k).}
#' Two families in the package have that shape, the inverse Gaussian and the
#' elastic net, and both reach all four orders through this one function.
#'
#' @details
#' # The two sums
#'
#' The Leibniz rule splits the positions of the multi-index \eqn{I} between the
#' weight and the tail. The weight side is
#' \eqn{\partial^{S} e^{w} = e^{w} B_{S}(w)}, the complete Bell polynomial in
#' the partials of the **log** weight; the tail side is one Faa di Bruno pass
#' over \eqn{x}, whose inner derivatives are the Hermite factors of
#' [phi_hermite()] times \eqn{\varphi(x)}. Nothing is transcribed: both sums
#' run on the package's own partition enumeration.
#'
#' # Why the log weight
#'
#' The weight is never formed on its own. It is combined with the tail as
#' `exp(w + pnorm(x, log.p = TRUE))`, because for the inverse Gaussian
#' \eqn{e^{2/(\phi\mu)}} overflows exactly where \eqn{\Phi(b)} underflows: at
#' \eqn{\mu = 0.01}, \eqn{\phi = 0.1} the weight is `Inf` while the product is
#' about \eqn{3\times10^{-106}}, and the fourth derivative there comes back at
#' \eqn{3\times10^{-91}}.
#'
#' @section Notation:
#' \eqn{F} is the distribution function, \eqn{\Phi} and \eqn{\varphi} the
#' standard normal distribution and density, \eqn{w} a log weight, \eqn{x} a
#' tail argument and \eqn{B_S} the complete Bell polynomial.
#'
#' @param distrib An object inheriting from `distrib`. Its `params` name and
#'   order the components.
#' @param q A numeric vector of quantiles.
#' @param order The derivative order, 1 to 4.
#' @param terms A list of terms. Each is a list with `sign` (\eqn{\pm 1}),
#'   `logw` (the log weight), `wderiv` (a function of a block of parameter
#'   names returning that partial of the log weight), `x` (the tail argument)
#'   and `xderiv` (the same for `x`). A term whose weight is 1 passes a `logw`
#'   of zero and a `wderiv` returning zero.
#'
#' @return A named list of numeric vectors, derivatives of \eqn{F} itself on
#'   the natural scale, keyed as
#'   [`deriv_names(distrib@params, order)`][deriv_names].
#'
#' @seealso [register_phi_terms_cdf()], which turns a term function into the
#'   four methods; [phi_hermite()]; [separable_deriv()] for the commonest
#'   shape of `wderiv` and `xderiv`.
#'
#' @keywords internal
phi_terms_cdf_deriv_k <- function(distrib, q, order, terms) {
  params <- distrib@params
  n <- length(q)

  one_term <- function(tm, I) {
    # exp(w) times the tail and times the density, both formed on the log
    # scale so that a large weight and a small tail cancel before either is
    # exponentiated
    ewPhi <- exp(tm$logw + stats::pnorm(tm$x, log.p = TRUE))
    ewphi <- exp(tm$logw + stats::dnorm(tm$x, log = TRUE))
    pos <- seq_along(I)
    acc <- rep(0, n)
    for (mask in seq_len(bitwShiftL(1L, length(pos))) - 1L) {
      take <- pos[bitwAnd(mask, bitwShiftL(1L, pos - 1L)) > 0L]
      rest <- setdiff(pos, take)
      wfac <- if (!length(take)) 1 else bell_f_ratio(params[I[take]], tm$wderiv)
      tail <- if (!length(rest)) {
        ewPhi
      } else {
        s <- rep(0, n)
        for (part in index_partitions(params[I[rest]])) {
          term <- ewphi * phi_hermite(tm$x, length(part))
          for (b in part) term <- term * tm$xderiv(b)
          s <- s + term
        }
        s
      }
      acc <- acc + wfac * tail
    }
    tm$sign * acc
  }

  idx <- deriv_indices(params, order)
  out <- lapply(idx, function(I) {
    acc <- rep(0, n)
    for (tm in terms) acc <- acc + one_term(tm, I)
    acc
  })
  stats::setNames(out, deriv_names(params, order))
}

#' Register All Four CDF Orders on a Sum of Normal Tails
#'
#' @description
#' Turns one term function into the four S7 methods a family of the shape
#' \eqn{F = c_0 + \sum_k s_k e^{w_k}\Phi(x_k)} needs, so that the family states
#' its terms once instead of four times. The inverse Gaussian and the elastic
#' net are the two families registered through it.
#'
#' @details
#' All four orders are registered, [distrib_grad_cdf()] included, so these
#' families take the closed route from the first order up and never reach a
#' stencil. `force(o)` inside the factory is what keeps the four registrations
#' from sharing one order.
#'
#' @param cls The S7 class to register on.
#' @param term_fn A function of `(distrib, q, theta)` returning the list of
#'   terms [phi_terms_cdf_deriv_k()] documents.
#'
#' @return Invisibly `NULL`. Called for the registration.
#'
#' @seealso [phi_terms_cdf_deriv_k()], the body it registers;
#'   [distrib_grad_cdf.InvGauss1Distrib()] and
#'   [distrib_grad_cdf.EnetDistrib()], the two families.
#'
#' @keywords internal
register_phi_terms_cdf <- function(cls, term_fn) {
  make <- function(o) {
    force(o)
    function(distrib, q, theta, lower.tail = TRUE, log = TRUE, ...) {
      tms <- term_fn(distrib, q, theta)
      tabs <- lapply(seq_len(o),
                     function(k) phi_terms_cdf_deriv_k(distrib, q, k, tms))
      cdf_scale_k(distrib, distrib_cdf(distrib, q, theta), tabs, o,
                  lower.tail, log)
    }
  }
  S7::method(distrib_grad_cdf, cls) <- make(1L)
  S7::method(distrib_hess_cdf, cls) <- make(2L)
  S7::method(distrib_deriv3_cdf, cls) <- make(3L)
  S7::method(distrib_deriv4_cdf, cls) <- make(4L)
  invisible(NULL)
}

#' Partial Derivatives of a Product of a Location Term and a Scale Term
#'
#' @description
#' Builds the `wderiv` or `xderiv` callback [phi_terms_cdf_deriv_k()] wants,
#' for a quantity that factorizes as \eqn{u(\theta_1)\,v(\theta_2)}. A mixed
#' partial of such a product is the product of two one-variable derivatives,
#' so a block of parameter names is answered by counting how many of each it
#' holds.
#'
#' @details
#' Both families this file serves are separable in exactly this way. The
#' inverse Gaussian's \eqn{a}, \eqn{b} and \eqn{c} are each a function of the
#' mean times a function of the dispersion; the elastic net's \eqn{s} and
#' \eqn{x} are each a function of \eqn{\lambda} times a function of
#' \eqn{\alpha}. Separability is why four orders are cheap here: no
#' multivariate expansion is ever formed.
#'
#' @param nm A character vector of length 2, the two parameter names, the one
#'   `uderiv` differentiates first.
#' @param uderiv A function of a non-negative whole number \eqn{j} returning
#'   \eqn{\partial^j u/\partial\theta_1^j}.
#' @param vderiv The same for \eqn{v} in the second parameter.
#'
#' @return A function of a block, a character vector of parameter names, that
#'   returns the corresponding mixed partial. A block naming a parameter
#'   outside `nm` gives the order-zero factor for both, which is the value
#'   itself; no caller does that.
#'
#' @seealso [phi_terms_cdf_deriv_k()], whose `wderiv` and `xderiv` this builds;
#'   [dpow_affine()] for the other shape.
#'
#' @keywords internal
separable_deriv <- function(nm, uderiv, vderiv) {
  function(block) {
    uderiv(sum(block == nm[1L])) * vderiv(sum(block == nm[2L]))
  }
}

#' @title Inverse Gaussian Log-CDF Derivatives
#' @name distrib_grad_cdf.InvGauss1Distrib
#'
#' @description
#' Closed form at every order from one to four. Unusually for a positive
#' family, the inverse Gaussian's distribution function is elementary,
#' \deqn{F(q) = \Phi(a) + e^{c}\,\Phi(b), \qquad
#'       a = \frac{q/\mu - 1}{\sqrt{\phi q}}, \quad
#'       b = -\frac{q/\mu + 1}{\sqrt{\phi q}}, \quad
#'       c = \frac{2}{\phi\mu},}
#' so it can simply be differentiated. This page's four registrations are made
#' together by [register_phi_terms_cdf()].
#'
#' @details
#' # Why four orders are cheap here
#'
#' Each of \eqn{a}, \eqn{b} and \eqn{c} is a product of a function of the mean
#' and a function of the dispersion, so every mixed partial is a product of two
#' one-variable derivatives and no multivariate expansion is formed. The two
#' terms are then a Leibniz split between the weight and the tail, which
#' [phi_terms_cdf_deriv_k()] runs on the package's own partition enumeration.
#'
#' # The overflow, and how it is avoided
#'
#' \eqn{e^{c}} is `Inf` at ordinary settings: the exponent is 2000 at
#' \eqn{\mu = 0.01}, \eqn{\phi = 0.1}. That is exactly where \eqn{\Phi(b)}
#' underflows, so the product is finite and neither factor is. The weight and
#' the tail are combined as `exp(c + pnorm(b, log.p = TRUE))`, and at that
#' setting the gradient comes back at \eqn{3\times10^{-106}} and the fourth
#' derivative at \eqn{3\times10^{-91}}.
#'
#' # What the closed route is worth
#'
#' Against a product stencil on the same cdf: \eqn{1.0\times10^{-10}} at order
#' 1, \eqn{8.8\times10^{-7}} at order 2, \eqn{1.8\times10^{-5}} at order 3 and
#' \eqn{3.7\times10^{-4}} at order 4. The gap is the stencil's error, and it is
#' the reason the family registers all four.
#'
#' @section Notation:
#' \eqn{\mu > 0} is the mean, \eqn{\phi > 0} the dispersion, \eqn{\Phi} the
#' standard normal distribution function and \eqn{F} the inverse Gaussian's.
#'
#' @param distrib An `InvGauss1Distrib` object, from [invgauss1_distrib()].
#' @param q A numeric vector of quantiles, positive.
#' @param theta A named list with components `mu` (positive) and `phi`
#'   (positive), each a numeric vector of length 1 or `n`.
#' @param lower.tail Is the lower tail wanted? A single logical, `TRUE` by
#'   default.
#' @param log Are derivatives of the log probability wanted? A single logical,
#'   `TRUE` by default.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A named list of numeric vectors of the order the generic asked for,
#'   keyed as [`deriv_names(distrib@params, order)`][deriv_names]: two
#'   components for the gradient, three for the Hessian, four at order 3 and
#'   five at order 4.
#'
#' @seealso [phi_terms_cdf_deriv_k()] for the construction;
#'   [register_phi_terms_cdf()], which makes the four registrations;
#'   [distrib_grad_cdf.EnetDistrib()], the other family of this shape;
#'   [invgauss1_distrib()].
#'
#' @examples
#' d <- invgauss1_distrib()
#' q <- c(0.5, 2, 5)
#' th <- list(mu = 2, phi = 0.5)
#'
#' # Against a central difference of the cdf, which shares no arithmetic.
#' fd <- numerical_cdf_deriv(d, q, th, order = 1)
#' max(abs(unlist(distrib_grad_cdf(d, q, th, log = FALSE)) / unlist(fd) - 1))
#'
#' # Finite where the weight alone is not: exp(2 / (phi mu)) is Inf here.
#' exp(2 / (0.01 * 0.1))
#' distrib_grad_cdf(d, 0.02, list(mu = 0.01, phi = 0.1), log = FALSE)
#'
#' @keywords internal
register_phi_terms_cdf(InvGauss1Distrib, function(distrib, q, theta) {
  mu <- theta[[1]]
  ph <- theta[[2]]
  nm <- distrib@params
  n <- length(q)
  rt <- sqrt(ph * q)
  a <- (q / mu - 1) / rt
  b <- -(q / mu + 1) / rt
  cval <- 2 / (ph * mu)

  # A = q/mu - 1 and C = q/mu + 1 differ only at order zero
  A <- function(j) if (j == 0L) q / mu - 1 else q * (-1)^j * factorial(j) / mu^(1 + j)
  C <- function(j) if (j == 0L) q / mu + 1 else q * (-1)^j * factorial(j) / mu^(1 + j)
  # B = (phi q)^(-1/2)
  B <- function(k) {
    coef <- if (k == 0L) 1 else prod(-0.5 - seq_len(k) + 1)
    coef * q^-0.5 * ph^(-0.5 - k)
  }
  # 1/mu and 1/phi, for c = 2/(phi mu)
  M <- function(j) (-1)^j * factorial(j) / mu^(1 + j)
  P <- function(k) (-1)^k * factorial(k) / ph^(1 + k)

  list(
    list(sign = 1, logw = rep(0, n), wderiv = function(block) rep(0, n),
         x = a, xderiv = separable_deriv(nm, A, B)),
    list(sign = 1, logw = rep_len(cval, n),
         wderiv = separable_deriv(nm, function(j) 2 * M(j), P),
         x = b,
         xderiv = separable_deriv(nm, function(j) -C(j), B))
  )
})


#' Derivatives of a Power of an Affine Function
#'
#' @description
#' Returns \eqn{\partial^k u^p/\partial\theta^k} where \eqn{u} is affine in
#' \eqn{\theta} with slope `inner`, which is
#' \eqn{(\mathrm{inner})^k\,p(p-1)\cdots(p-k+1)\,u^{p-k}}. The elastic net's
#' tail arguments are powers of affine functions of its two hyperparameters, so
#' this supplies their one-variable derivatives to [separable_deriv()].
#'
#' @param u A numeric vector, the affine function evaluated at the parameter.
#' @param p The power, a single number, not necessarily a whole one.
#' @param k The derivative order, a non-negative whole number. Zero returns
#'   \eqn{u^p} itself.
#' @param inner The slope of the affine function, a single number or a numeric
#'   vector recyclable against `u`.
#'
#' @return A numeric vector the length of `u`.
#'
#' @seealso [separable_deriv()], which combines two of these;
#'   [distrib_grad_cdf.EnetDistrib()], the consumer.
#'
#' @examples
#' # The second derivative of (2 theta)^3 at theta = 1 is 8 * 3 * 2 * 2 = 24.
#' distributions7:::dpow_affine(u = 2, p = 3, k = 2, inner = 2)
#'
#' @keywords internal
dpow_affine <- function(u, p, k, inner) {
  if (k == 0L) return(u^p)
  inner^k * prod(p - seq_len(k) + 1) * u^(p - k)
}

#' @title Elastic Net Log-CDF Derivatives
#' @name distrib_grad_cdf.EnetDistrib
#'
#' @description
#' Closed form at every order from one to four. Each half of the distribution
#' function is a truncated Gaussian: with \eqn{z = q - \mu}, \eqn{s = \sqrt c}
#' and \eqn{x = a/\sqrt c} it is \eqn{e^{w}\Phi(X)} below the location and
#' \eqn{1 - e^{w}\Phi(X)} above it, for \eqn{X = \pm sz - x} and a weight
#' \eqn{w = -\log M(x) + x^2/2 + \mathrm{const}} written through the Mills
#' ratio the family already carries. The four registrations are made together
#' by [register_phi_terms_cdf()].
#'
#' @details
#' # Why four orders are cheap here
#'
#' Both \eqn{s} and \eqn{x} are a function of \eqn{\lambda} times a function of
#' \eqn{\alpha}, so their mixed partials are products of one-variable ones
#' through [separable_deriv()], and the one-variable ones are powers of affine
#' functions through [dpow_affine()]. No multivariate expansion is formed at
#' any order.
#'
#' # The kink at the location
#'
#' The location is the non-regular direction, as it is in the Laplace this
#' family contains: the second derivative in \eqn{\mu} carries a point mass at
#' \eqn{q = \mu}, and the formulas hold on either side of it. A numerical check
#' that straddles the location is checking the arithmetic against a reference
#' that is not valid there.
#'
#' # What the closed route is worth
#'
#' Against a product stencil on the same cdf: \eqn{1.6\times10^{-9}} at order
#' 1, \eqn{3.5\times10^{-7}} at order 2, \eqn{3.0\times10^{-5}} at order 3 and
#' \eqn{3.9\times10^{-4}} at order 4.
#'
#' @section Notation:
#' \eqn{\mu} is the location, \eqn{\lambda > 0} and \eqn{\alpha \in (0,1)} the
#' two hyperparameters, \eqn{a = \lambda\alpha} and \eqn{c = \lambda(1-\alpha)}
#' the two rates, \eqn{M} the Mills ratio, \eqn{\Phi} the standard normal
#' distribution function and \eqn{F} the elastic net's.
#'
#' @param distrib An `EnetDistrib` object, from [enet_distrib()].
#' @param q A numeric vector of quantiles.
#' @param theta A named list with components `mu` (any real value), `lambda`
#'   (positive) and `alpha` (strictly between 0 and 1), each a numeric vector
#'   of length 1 or `n`.
#' @param lower.tail Is the lower tail wanted? A single logical, `TRUE` by
#'   default.
#' @param log Are derivatives of the log probability wanted? A single logical,
#'   `TRUE` by default.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A named list of numeric vectors of the order the generic asked for,
#'   keyed as [`deriv_names(distrib@params, order)`][deriv_names]: three
#'   components for the gradient, six for the Hessian, ten at order 3 and
#'   fifteen at order 4.
#'
#' @seealso [phi_terms_cdf_deriv_k()] for the construction;
#'   [separable_deriv()] and [dpow_affine()] for the partials;
#'   [distrib_grad_cdf.InvGauss1Distrib()], the other family of this shape;
#'   [enet_distrib()].
#'
#' @examples
#' d <- enet_distrib()
#' q <- c(-1, 0.5, 2)
#' th <- list(mu = 0, lambda = 1, alpha = 0.5)
#'
#' # Against a central difference of the cdf, which shares no arithmetic.
#' fd <- numerical_cdf_deriv(d, q, th, order = 1)
#' max(abs(unlist(distrib_grad_cdf(d, q, th, log = FALSE)) / unlist(fd) - 1))
#'
#' # Three, six, ten and fifteen components as the order rises.
#' lengths(list(distrib_grad_cdf(d, q, th),
#'              distrib_hess_cdf(d, q, th),
#'              distrib_deriv3_cdf(d, q, th),
#'              distrib_deriv4_cdf(d, q, th)))
#'
#' @keywords internal
register_phi_terms_cdf(EnetDistrib, function(distrib, q, theta) {
  p <- .enet_g_higher(.enet_parts(theta))
  nm <- distrib@params
  n <- length(q)
  lam <- p$lam
  al <- p$al
  x0 <- p$x
  s <- sqrt(p$c)
  z <- q - p$mu
  eps <- ifelse(z <= 0, 1, -1)

  # s = sqrt(lambda) sqrt(1 - alpha) and x = sqrt(lambda) alpha / sqrt(1 - alpha)
  sqrt_lam <- function(j) dpow_affine(lam, 0.5, j, 1)
  v_s <- function(k) dpow_affine(1 - al, 0.5, k, -1)
  v_x <- function(k) {
    # d^k of alpha (1 - alpha)^(-1/2), by the product rule. The second term is
    # guarded rather than multiplied by zero: R evaluates it either way, and
    # dpow_affine has no order below zero.
    out <- al * dpow_affine(1 - al, -0.5, k, -1)
    if (k > 0L) out <- out + k * dpow_affine(1 - al, -0.5, k - 1L, -1)
    out
  }
  # a block naming the location contributes nothing to either
  sep <- function(vfun) function(block) {
    if (any(block == nm[1L])) return(rep(0, n))
    rep_len(sqrt_lam(sum(block == nm[2L])) * vfun(sum(block == nm[3L])), n)
  }
  ds <- sep(v_s)
  dx <- sep(v_x)

  # w(x) = -log M(x) + x^2/2 + const, so its derivatives in x are written
  # through G = dlogM/dx, which the family already computes
  wx <- list(-p$g + x0, -p$dg + 1, -p$d2g, -p$d3g)
  wderiv <- function(block) {
    if (any(block == nm[1L])) return(rep(0, n))
    acc <- rep(0, n)
    for (part in index_partitions(block)) {
      term <- rep_len(wx[[length(part)]], n)
      for (b in part) term <- term * dx(b)
      acc <- acc + term
    }
    acc
  }

  # X = eps s (q - mu) - x, linear in the location
  xderiv <- function(block) {
    nmu <- sum(block == nm[1L])
    if (nmu >= 2L) return(rep(0, n))
    rest <- block[block != nm[1L]]
    if (nmu == 1L) {
      if (!length(rest)) return(-eps * rep_len(s, n))
      return(-eps * ds(rest))
    }
    eps * z * ds(block) - dx(block)
  }

  logw <- -.enet_logQ(x0) - log(2)
  list(list(sign = eps, logw = rep_len(logw, n),
            wderiv = wderiv, x = eps * s * z - x0, xderiv = xderiv))
})
