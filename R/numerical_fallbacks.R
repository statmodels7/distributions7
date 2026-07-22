#' @include distrib.R generics.R utility_functions.R

# =========================== CONTINUOUS FALLBACKS ===========================
#
# Default methods registered on `continuous_distrib`: any continuous subclass
# that implements only distrib_pdf automatically gets a CDF (numerical
# integration), a quantile function (root-finding on the CDF) and an RNG
# (inverse transform sampling). Subclasses with closed forms override these
# through normal S7 dispatch and are both faster and more accurate.

# Internal: locate an interior high-density anchor point (approximate mode)
# used to split integrals, to scale root-finding brackets and to recentre the
# GRoU kernel. The support is compactified to (0, 1) and the maximum is found by
# repeatedly evaluating the log-density on a grid and keeping the two cells
# around the largest value: for a unimodal density that bracket provably still
# contains the mode, and it shrinks by a factor of 64 per pass.
#
# A golden-section search on the compactified scale is not accurate enough here.
# Its tolerance is expressed in t, and the derivative of the compactification can
# be enormous -- with the tangent map dy/dt is of order y^2, so the default
# tolerance of about 1e-4 turns into an error of several hundred units in y for a
# density centred far from the origin. The grid refinement instead stops on the
# width of the bracket measured in y.
find_pdf_anchor <- function(distrib, theta) {
  find_lp_anchor(
    function(y) distrib_pdf(distrib, y, theta, log = TRUE),
    distrib@bounds
  )
}

# The same search expressed on a bare log-density over an interval, so that it
# can also be applied to a reparameterised density that has no distrib object.
find_lp_anchor <- function(lp_raw, b) {
  lp <- function(y) {
    v <- lp_raw(y)
    v[is.na(v)] <- -Inf
    v
  }

  to_y <- if (all(is.finite(b))) {
    function(t) b[1] + t * (b[2] - b[1])
  } else if (is.finite(b[1])) {
    function(t) b[1] + t / (1 - t)
  } else if (is.finite(b[2])) {
    function(t) b[2] - (1 - t) / t
  } else {
    function(t) tan(pi * (t - 0.5))
  }

  lo <- 1e-10
  hi <- 1 - 1e-10
  n_grid <- 129L

  for (i in seq_len(60L)) {
    t <- seq(lo, hi, length.out = n_grid)
    v <- lp(to_y(t))
    if (all(!is.finite(v))) break
    k <- which.max(v)
    lo <- t[max(1L, k - 1L)]
    hi <- t[min(n_grid, k + 1L)]
    y_lo <- to_y(lo)
    y_hi <- to_y(hi)
    if (is.finite(y_hi - y_lo) &&
        y_hi - y_lo <= 1e-10 * max(1, abs(y_lo), abs(y_hi))) {
      break
    }
  }

  to_y((lo + hi) / 2)
}

#' @title Default Numerical CDF for Continuous Distributions
#' @name distrib_cdf.continuous_distrib
#' @description
#' Fallback method: continuous distributions that do not implement an analytical
#' CDF get one by numerical integration of \code{distrib_pdf}. An approximate mode
#' is located first and the integral is taken over the side of the mode containing
#' \eqn{q} (using the complement for the other side), so that the quadrature nodes
#' concentrate where the probability mass is.
#' @param distrib An object inheriting from class \code{"continuous_distrib"}.
#' @param q A numeric vector of quantiles.
#' @param theta A named list of parameters.
#' @param lower.tail Logical; if \code{TRUE} (default), probabilities are \eqn{P(Y \le q)}.
#' @param log.p Logical; if \code{TRUE}, probabilities are returned as logs.
#' @keywords internal
S7::method(distrib_cdf, continuous_distrib) <- function(distrib, q, theta, lower.tail = TRUE, log.p = FALSE) {
  b <- distrib@bounds
  all_params <- expand_params(c(list(.q = q), theta))
  rows <- transpose_params(all_params)

  prev_th <- NULL
  m <- NULL

  res <- vapply(rows, function(r) {
    r <- as.list(r)
    th <- r[distrib@params]
    qi <- r$.q

    if (is.na(qi)) return(NA_real_)
    if (qi <= b[1]) return(0)
    if (qi >= b[2]) return(1)

    if (!identical(th, prev_th)) {
      m <<- find_pdf_anchor(distrib, th)
      prev_th <<- th
    }

    if (qi <= m) {
      stats::integrate(function(t) distrib_pdf(distrib, t, th), b[1], qi)$value
    } else {
      1 - stats::integrate(function(t) distrib_pdf(distrib, t, th), qi, b[2])$value
    }
  }, numeric(1))

  res <- pmin(pmax(res, 0), 1)
  if (!lower.tail) res <- 1 - res
  if (log.p) log(res) else res
}

#' @title Default Numerical Quantile Function for Continuous Distributions
#' @name distrib_quantile.continuous_distrib
#' @description
#' Fallback method: continuous distributions that do not implement an analytical
#' quantile function get one by root-finding on \code{\link{distrib_cdf}} (which may
#' itself be the numerical fallback). Brackets start from an approximate mode and
#' expand geometrically, with the step scaled by the density height at the mode.
#' @param distrib An object inheriting from class \code{"continuous_distrib"}.
#' @param p A numeric vector of probabilities.
#' @param theta A named list of parameters.
#' @param lower.tail Logical; if \code{TRUE} (default), probabilities are \eqn{P(Y \le p)}.
#' @param log.p Logical; if \code{TRUE}, probabilities are given as logs.
#' @keywords internal
S7::method(distrib_quantile, continuous_distrib) <- function(distrib, p, theta, lower.tail = TRUE, log.p = FALSE) {
  if (log.p) p <- exp(p)
  if (!lower.tail) p <- 1 - p

  b <- distrib@bounds
  all_params <- expand_params(c(list(.p = p), theta))
  rows <- transpose_params(all_params)

  prev_th <- NULL
  m <- NULL
  scale_m <- NULL
  Fm <- NULL

  vapply(rows, function(r) {
    r <- as.list(r)
    th <- r[distrib@params]
    pi <- r$.p

    if (is.na(pi) || pi < 0 || pi > 1) return(NaN)
    if (pi == 0) return(b[1])
    if (pi == 1) return(b[2])

    if (!identical(th, prev_th)) {
      m <<- find_pdf_anchor(distrib, th)
      # Typical width of the bulk of a unimodal density ~ 1 / f(mode)
      scale_m <<- 1 / max(distrib_pdf(distrib, m, th), 1e-12)
      Fm <<- distrib_cdf(distrib, m, th)
      prev_th <<- th
    }

    f <- function(q) distrib_cdf(distrib, q, th) - pi

    if (pi <= Fm) {
      hi <- m
      if (is.finite(b[1])) {
        lo <- b[1]
      } else {
        s <- scale_m
        lo <- m - s
        it <- 0L
        while (f(lo) > 0 && it < 80L) {
          s <- 2 * s
          lo <- m - s
          it <- it + 1L
        }
      }
    } else {
      lo <- m
      if (is.finite(b[2])) {
        hi <- b[2]
      } else {
        s <- scale_m
        hi <- m + s
        it <- 0L
        while (f(hi) < 0 && it < 80L) {
          s <- 2 * s
          hi <- m + s
          it <- it + 1L
        }
      }
    }

    stats::uniroot(f, lower = lo, upper = hi,
      tol = .Machine$double.eps^0.5 * max(1, abs(m)))$root
  }, numeric(1))
}

# Internal: TRUE when the object gets its quantile function from a class-specific
# method rather than from the numerical fallback registered just above. S7 records
# the class a method was registered on in the method's `signature` attribute, so
# inherited fallbacks are recognised by that class being `continuous_distrib`
# itself. Used to decide whether inverse transform sampling is cheap.
has_analytic_quantile <- function(distrib) {
  m <- tryCatch(
    S7::method(distrib_quantile, S7::S7_class(distrib)),
    error = function(e) NULL
  )
  !is.null(m) && !identical(attr(m, "signature")[[1]], continuous_distrib)
}

#' @title Generalized Ratio-of-Uniforms Sampling
#'
#' @description
#' Draws from a continuous distribution using the Generalized Ratio-of-Uniforms
#' (GRoU) method. The sampler needs nothing but the (log) density: no CDF, no
#' quantile function and no inversion. It is the engine behind the default
#' \code{\link{distrib_rng}} method for continuous distributions that provide
#' neither a native RNG nor an analytical quantile function.
#'
#' @param distrib An object inheriting from class \code{"continuous_distrib"}.
#' @param n Number of observations to generate.
#' @param theta A named list of parameters, each of length one.
#' @param r Numeric tuning parameter of the transformation power, default \code{2}.
#'   \code{r = 1} is the classical Ratio-of-Uniforms; larger values enclose
#'   heavier tails (\code{r = 2} keeps the acceptance region bounded for tails as
#'   heavy as the Cauchy's, which is why it is the default).
#'
#' @details
#' Write \eqn{K} for the density kernel. If \eqn{(U, V)} is uniform over
#' \deqn{A_r = \left\{ (u, v) : 0 < u \leq K\!\left(v / u^{r}\right)^{1/(r+1)} \right\}}
#' then \eqn{Y = V / U^{r}} has density proportional to \eqn{K}. Sampling therefore
#' reduces to drawing uniformly in a bounding rectangle
#' \eqn{[0, u_{\max}] \times [v_{\min}, v_{\max}]} and keeping the pairs that fall
#' in \eqn{A_r}, i.e. those with
#' \deqn{(r+1)\log U \leq \log K(Y).}
#'
#' Two devices make this numerically safe for an arbitrary user-supplied density:
#' \itemize{
#'   \item \strong{Recentring.} The kernel is shifted to its mode,
#'     \eqn{K(z) \propto f(m + z)}, and the mode is added back to the draws. Without
#'     this a distribution located far from the origin (say \eqn{\mu = 1000}) gives a
#'     wildly elongated bounding rectangle and an acceptance rate close to zero.
#'     Recentring also splits the box exactly: \eqn{h(z) = z\,K(z)^{r/(r+1)}} is
#'     non-positive for \eqn{z \le 0} and non-negative for \eqn{z \ge 0}, so
#'     \eqn{v_{\min}} and \eqn{v_{\max}} are each found on one side of the mode.
#'   \item \strong{Normalisation.} The kernel is rescaled by its value at the mode, so
#'     that \eqn{\max K = 1} and \eqn{u_{\max} = 1}. Every quantity stays in a safe
#'     numerical range whatever the height of the density, and all computations are
#'     carried out on the log scale.
#' }
#'
#' The bounds \eqn{v_{\min}, v_{\max}} are obtained by expanding geometrically away
#' from the mode until \eqn{h} has clearly turned back towards zero (finite support
#' bounds are used directly), then combining \code{\link[stats]{optimize}} with a grid
#' search over the resulting bracket; the box is enclosing by construction for a
#' unimodal kernel. Candidates are generated and filtered in vectorised batches whose
#' size adapts to the observed acceptance rate.
#'
#' @section Unbounded densities:
#' The method needs the supremum of the kernel to be attained, so a density that
#' diverges at a finite edge of its support would have its acceptance region
#' clipped. Rather than give up, such an edge is transformed away.
#'
#' If \eqn{f(y) \sim |y-a|^{\alpha-1}} near the edge \eqn{a}, with \eqn{\alpha < 1},
#' then \eqn{X = |Y-a|^{1/\lambda}} has density proportional to
#' \eqn{x^{\lambda\alpha-1}} there, which is bounded as soon as
#' \eqn{\lambda\alpha > 1}. Sampling \eqn{X} and mapping back is exact: no
#' quadrature and no inversion are involved. The exponent \eqn{\alpha} is read off
#' the same probe that detects the divergence -- walking towards the edge in
#' decades lifts the log-density by \eqn{(1-\alpha)\log 10} per step -- and
#' \eqn{\lambda} is chosen from it. A Gamma of shape \eqn{0.4}, for instance, is
#' sampled through \eqn{X = Y^{1/4}}.
#'
#' A density diverging at \emph{both} edges, such as a Beta with both shapes below
#' one, is beyond any single power: flattening one edge steepens the other. What
#' does work is a map that behaves like a different power at each end,
#' \deqn{T(u) = \frac{u^{p}}{u^{p} + (1-u)^{q}}, \qquad Y = a + (b-a)\,T(U),}
#' monotone from 0 to 1, with a closed-form derivative, and carrying the exponents
#' \eqn{p\alpha - 1} and \eqn{q\beta - 1} at the two ends. Choosing
#' \eqn{p\alpha > 1} and \eqn{q\beta > 1} makes the density of \eqn{U} vanish at
#' both ends, turning a U-shaped density into a single-peaked one.
#'
#' Near an edge the variable itself stops being resolvable: at a non-zero edge the
#' spacing of doubles is absolute, so a slab of \eqn{u} collapses onto the edge and
#' the density evaluated there returns infinity. The distance from the edge is
#' still exact on the \eqn{u} scale, so the power law is used there instead,
#' calibrated once where the density can be evaluated. This keeps the mass rather
#' than rejecting it. It cannot place it any more finely than the arithmetic
#' allows -- for a Beta(0.9, 0.1) some 2.5\% of the distribution lies within one
#' unit in the last place of 1, and those draws come back equal to 1 exactly, which
#' is also what \code{\link[stats]{rbeta}} does.
#'
#' @return A numeric vector of length \code{n}.
#'
#' @section Requirements:
#' The kernel must be unimodal and the parameters in \code{theta} must be scalars.
#' Densities that diverge at one or at both edges of their support are handled by
#' the reparameterisations described below.
#'
#' Heavy tails are not an obstacle: with the default \code{r = 2} the sampler
#' handles a Student's t with half a degree of freedom and a Pareto with infinite
#' mean, and multimodal densities are usually accepted as well.
#'
#' @seealso \code{\link{distrib_rng}}, \code{\link{check_distrib}}
#'
#' @examples
#' # A distribution defined by its density alone still gets a fast sampler
#' MyDist <- S7::new_class("MyDist", parent = continuous_distrib)
#' S7::method(distrib_pdf, MyDist) <- function(distrib, y, theta, log = FALSE) {
#'   ld <- -log(2 * theta$b) - abs(y - theta$mu) / theta$b
#'   if (log) ld else exp(ld)
#' }
#' d <- MyDist(
#'   distrib_name = "my laplace", dimension = "univariate", bounds = c(-Inf, Inf),
#'   params = c("mu", "b"), params_interpretation = c(mu = "location", b = "scale"),
#'   n_params = 2, params_bounds = list(mu = c(-Inf, Inf), b = c(0, Inf)),
#'   link_params = list(
#'     mu = linkfunctions7::identity_link(),
#'     b = linkfunctions7::log_link()
#'   )
#' )
#' y <- rng_grou(d, 1000, list(mu = 2, b = 1.5))
#' c(mean = mean(y), var = var(y)) # ~ 2 and ~ 2 * 1.5^2
#'
#' @export
rng_grou <- function(distrib, n, theta, r = 2) {
  if (any(lengths(theta) != 1L)) {
    stop("rng_grou() requires scalar parameters; got vector-valued theta.", call. = FALSE)
  }
  if (n <= 0) return(numeric(0))

  lp <- function(y) distrib_pdf(distrib, y, theta, log = TRUE)
  b <- distrib@bounds
  div <- lp_edge_divergence(lp, b)

  if (all(is.na(div))) {
    return(grou_core(lp, b, n, r))
  }
  if (sum(!is.na(div)) > 1L) {
    return(grou_two_sided(lp, b, div, n, r))
  }

  # One divergent edge can be transformed away. If the density behaves like
  # |y - a|^(alpha - 1) near the edge a, then X = |Y - a|^(1/lambda) has density
  # proportional to x^(lambda*alpha - 1) there, which is bounded as soon as
  # lambda*alpha > 1. Sampling X by GRoU and mapping back is exact -- no
  # quadrature, no inversion -- and lambda follows from the exponent, which the
  # same probe that detected the divergence already estimates.
  at_lower <- !is.na(div[["lower"]])
  edge <- if (at_lower) b[1] else b[2]
  alpha <- if (at_lower) div[["lower"]] else div[["upper"]]
  lambda <- ceiling(1 / alpha) + 1

  to_y <- if (at_lower) function(x) edge + x^lambda else function(x) edge - x^lambda
  width <- abs((if (at_lower) b[2] else b[1]) - edge)

  lp_x <- function(x) {
    out <- rep(-Inf, length(x))
    # Once x is small enough that edge +- x^lambda rounds back onto the edge, the
    # original density is evaluated exactly at its singularity and returns
    # infinity, undoing the reparameterisation. The transformed density tends to
    # zero there (lambda*alpha > 1 by construction), so -Inf is both the correct
    # limit and the safe value.
    ok <- is.finite(x) & x > 0 & to_y(x) != edge
    if (any(ok)) {
      v <- lp(to_y(x[ok])) + log(lambda) + (lambda - 1) * log(x[ok])
      v[is.na(v)] <- -Inf
      out[ok] <- v
    }
    out
  }
  b_x <- c(0, if (is.finite(width)) width^(1 / lambda) else Inf)

  to_y(grou_core(lp_x, b_x, n, r))
}

# Internal: both edges of a (necessarily bounded) support diverge. A single
# power cannot repair that, since raising to a power flattens one edge while
# steepening the other. A map that behaves like a *different* power at each end
# does:
#
#   T(u) = u^p / (u^p + (1-u)^q),      Y = a + (b - a) T(U)
#
# It increases monotonically from 0 to 1, behaves like u^p on the left and like
# 1 - (1-u)^q on the right, and its derivative is available in closed form,
#
#   T'(u) = [p u^(p-1) (1-u)^q + q u^p (1-u)^(q-1)] / (u^p + (1-u)^q)^2 .
#
# The density of U therefore carries the exponents p*alpha - 1 and q*beta - 1 at
# the two ends and is bounded as soon as p*alpha > 1 and q*beta > 1 -- indeed it
# vanishes there, which turns the original U-shaped density into a single-peaked
# one, exactly what the sampler wants. Both exponents come from the same probe
# that detected the divergence, so nothing has to be searched for.
grou_two_sided <- function(lp, b, div, n, r) {
  a <- b[1]
  bb <- b[2]
  width <- bb - a
  p <- ceiling(1 / div[["lower"]]) + 1
  q <- ceiling(1 / div[["upper"]]) + 1

  log_sum <- function(x, y) {
    m <- pmax(x, y)
    m + log(exp(x - m) + exp(y - m))
  }
  t_of_u <- function(u) {
    la <- p * log(u)
    lb <- q * log1p(-u)
    exp(la - log_sum(la, lb))
  }
  log_dt <- function(u) {
    lu <- log(u)
    l1u <- log1p(-u)
    log_sum(log(p) + (p - 1) * lu + q * l1u,
            log(q) + p * lu + (q - 1) * l1u) - 2 * log_sum(p * lu, q * l1u)
  }

  # Distances from either edge, computed in u-space and therefore free of the
  # cancellation that reconstructing y would suffer.
  log_gap_lo <- function(u) log(width) + p * log(u) - log_sum(p * log(u), q * log1p(-u))
  log_gap_hi <- function(u) log(width) + q * log1p(-u) - log_sum(p * log(u), q * log1p(-u))

  al <- div[["lower"]]
  be <- div[["upper"]]
  tol_lo <- max(abs(a), 1) * 1e-12
  tol_hi <- max(abs(bb), 1) * 1e-12
  log_c_lo <- lp(a + tol_lo) - (al - 1) * log(tol_lo)
  log_c_hi <- lp(bb - tol_hi) - (be - 1) * log(tol_hi)

  lp_u <- function(u) {
    out <- rep(-Inf, length(u))
    ok <- is.finite(u) & u > 0 & u < 1
    if (!any(ok)) return(out)
    uu <- u[ok]
    lgl <- log_gap_lo(uu)
    lgh <- log_gap_hi(uu)
    v <- numeric(length(uu))

    # Near an edge, y itself is no longer resolvable -- at a non-zero edge the
    # spacing of doubles is absolute, so a whole slab of u-space collapses onto
    # the edge and evaluating the density there would return infinity. The gap
    # from the edge, however, is exact in u-space, and the density is known there
    # to behave like gap^(alpha - 1). Using that law, calibrated once at a point
    # that *is* resolvable, keeps this mass instead of rejecting it.
    near_lo <- lgl < log(tol_lo)
    near_hi <- lgh < log(tol_hi)
    mid <- !near_lo & !near_hi
    if (any(mid)) v[mid] <- lp(a + width * t_of_u(uu[mid]))
    if (any(near_lo)) v[near_lo] <- log_c_lo + (al - 1) * lgl[near_lo]
    if (any(near_hi)) v[near_hi] <- log_c_hi + (be - 1) * lgh[near_hi]

    v <- v + log(width) + log_dt(uu)
    v[is.na(v)] <- -Inf
    out[ok] <- v
    out
  }

  a + width * t_of_u(grou_core(lp_u, c(0, 1), n, r))
}

# Internal: for each finite edge of the support, the exponent alpha of a
# divergence f(y) ~ |y - edge|^(alpha - 1), or NA when the density stays bounded
# there. Walking towards the edge in decades lifts the log-density by
# (1 - alpha) * log(10) per step when it diverges, and by an amount that dies
# away when it does not.
lp_edge_divergence <- function(lp, b) {
  probe <- function(edge, inward) {
    step <- inward * 1e-3 * max(abs(edge), 1) * 10^(-(0:11))
    v <- lp(edge + step)
    v <- v[is.finite(v)]
    if (length(v) < 5L) return(NA_real_)
    dv <- diff(v)
    if (!all(dv > 0) || min(dv[(length(dv) - 2L):length(dv)]) <= 1e-6) return(NA_real_)
    alpha <- 1 - stats::median(dv) / log(10)
    min(max(alpha, 1e-6), 1 - 1e-6)
  }
  c(lower = if (is.finite(b[1])) probe(b[1], 1) else NA_real_,
    upper = if (is.finite(b[2])) probe(b[2], -1) else NA_real_)
}

# Internal: the Generalized Ratio-of-Uniforms itself, on a bare log-density over
# an interval. Kept separate from the distrib object so that it can also be run
# on a reparameterised density.
grou_core <- function(lp, b, n, r) {
  m <- find_lp_anchor(lp, b)
  lp_max <- lp(m)
  if (!is.finite(lp_max)) {
    stop("GRoU needs a bounded density: the log-density at the mode is ", lp_max, ".", call. = FALSE)
  }

  # Log-kernel recentred at the mode and normalised to a maximum of 1
  lk <- function(z) {
    out <- lp(m + z) - lp_max
    out[is.na(out)] <- -Inf
    out
  }
  zl <- b[1] - m
  zu <- b[2] - m

  rr <- r / (r + 1)
  h <- function(z) {
    out <- numeric(length(z))
    ok <- is.finite(z) & z >= zl & z <= zu
    if (any(ok)) out[ok] <- z[ok] * exp(rr * lk(z[ok]))
    out[!is.finite(out)] <- 0
    out
  }

  # Typical width of the bulk of a unimodal density, as in the quantile fallback
  scale_m <- 1 / max(exp(lp_max), 1e-12)

  # Expand away from the mode until h has turned back towards zero
  far_end <- function(bound, sgn) {
    if (is.finite(bound)) return(bound)
    z <- sgn * scale_m
    best <- 0
    for (i in seq_len(100L)) {
      v <- abs(h(z))
      if (v >= best) {
        best <- v
      } else if (v <= 0.5 * best) {
        break
      }
      if (!is.finite(2 * z)) break
      z <- 2 * z
    }
    z
  }

  z_right <- far_end(zu, 1)
  z_left <- far_end(zl, -1)

  # optimize() plus a grid, so that a missed interior peak cannot shrink the box
  extreme <- function(lo, hi, maximum) {
    if (!(hi > lo)) return(0)
    g <- h(seq(lo, hi, length.out = 64L))
    o <- stats::optimize(h, lower = lo, upper = hi, maximum = maximum)$objective
    if (maximum) max(c(g, o, 0)) else min(c(g, o, 0))
  }

  v_max <- extreme(0, z_right, TRUE)
  v_min <- extreme(z_left, 0, FALSE)
  if (!is.finite(v_max) || !is.finite(v_min) || v_max <= v_min) {
    stop("GRoU could not build a bounding rectangle for this density.", call. = FALSE)
  }

  # find_pdf_anchor is only approximate, so allow for a kernel slightly above 1
  lk_grid <- lk(c(seq(z_left, 0, length.out = 32L), seq(0, z_right, length.out = 32L)))
  u_max <- exp(max(0, max(lk_grid[is.finite(lk_grid)])) / (r + 1))

  out <- numeric(n)
  got <- 0L
  rate <- NA_real_
  rounds <- 0L

  while (got < n) {
    rounds <- rounds + 1L
    if (rounds > 1000L) {
      stop("GRoU failed to accept any draw in 1000 rounds; the density may be multimodal.", call. = FALSE)
    }
    need <- n - got
    size <- if (is.na(rate) || rate <= 0) max(64L, 2L * need) else max(64L, ceiling(1.1 * need / rate))
    size <- as.integer(min(size, 1e7))

    u <- stats::runif(size, 0, u_max)
    v <- stats::runif(size, v_min, v_max)
    z <- v / u^r

    ok <- is.finite(z) & z >= zl & z <= zu
    accepted <- numeric(0)
    if (any(ok)) {
      z <- z[ok]
      keep <- (r + 1) * log(u[ok]) <= lk(z)
      accepted <- z[!is.na(keep) & keep]
    }

    rate <- if (is.na(rate)) length(accepted) / size else 0.5 * rate + 0.5 * length(accepted) / size

    if (length(accepted)) {
      take <- min(length(accepted), need)
      out[(got + 1L):(got + take)] <- accepted[seq_len(take)]
      got <- got + take
      rounds <- 0L
    }
  }

  m + out
}

#' @title Default Numerical RNG for Continuous Distributions
#' @name distrib_rng.continuous_distrib
#' @description
#' Fallback method for continuous distributions that do not implement a native RNG.
#' Two strategies are available and the method picks between them automatically:
#' \itemize{
#'   \item \strong{Inverse transform}, \code{distrib_quantile(distrib, runif(n), theta)},
#'     when the distribution provides its own quantile function. This is exact and
#'     costs one quantile evaluation per draw.
#'   \item \strong{Generalized Ratio-of-Uniforms} (\code{\link{rng_grou}}) otherwise.
#'     Inverting a purely numerical CDF costs one root-finding step per draw, which
#'     makes simulation-based tools (\code{approx = "mc"}, \code{\link{check_distrib}})
#'     impractical; GRoU only evaluates the density, so it is orders of magnitude faster.
#' }
#' GRoU requires a bounded, unimodal density; if it cannot build its bounding
#' rectangle the method warns and reverts to inverse transform sampling. Vector-valued
#' \code{theta} is handled by grouping the draws by distinct parameter values.
#' @param distrib An object inheriting from class \code{"continuous_distrib"}.
#' @param n Number of observations to generate.
#' @param theta A named list of parameters.
#' @keywords internal
S7::method(distrib_rng, continuous_distrib) <- function(distrib, n, theta) {
  if (has_analytic_quantile(distrib)) {
    return(distrib_quantile(distrib, stats::runif(n), theta))
  }

  by_inversion <- function() distrib_quantile(distrib, stats::runif(n), theta)
  grou <- function(th, k) {
    tryCatch(rng_grou(distrib, k, th), error = function(e) {
      warning("GRoU sampling failed (", conditionMessage(e),
        "); falling back to inverse transform sampling.", call. = FALSE)
      NULL
    })
  }

  if (all(lengths(theta) == 1L)) {
    res <- grou(theta, n)
    return(if (is.null(res)) by_inversion() else res)
  }

  # Per-observation parameters: draw group by group
  th <- lapply(expand_params(theta), function(x) rep_len(x, n))
  key <- do.call(paste, c(th, sep = "\r"))
  out <- numeric(n)
  for (k in unique(key)) {
    idx <- which(key == k)
    res <- grou(lapply(th, function(x) x[idx[1]]), length(idx))
    if (is.null(res)) return(by_inversion())
    out[idx] <- res
  }
  out
}

# =========================== DISCRETE FALLBACKS ===========================
#
# Default methods registered on `discrete_distrib`: any discrete subclass that
# implements only the pmf (distrib_pdf) gets CDF, quantile and RNG built from a
# cumulative-probability table over the support. The support's lower bound must
# be finite (which covers all standard count distributions).

# Internal: cumulative pmf table from the lower bound, grown geometrically
# until it covers probability `need_p` and support point `need_k`.
disc_cum_table <- function(distrib, theta, need_p = -Inf, need_k = -Inf, kmax = 1e6) {
  b <- distrib@bounds
  if (!is.finite(b[1])) {
    stop("The numerical fallback for discrete distributions requires a finite lower bound.", call. = FALSE)
  }

  need_p <- min(need_p, 1 - 1e-12)
  size <- 256L

  repeat {
    upper_k <- min(b[1] + size - 1, b[2])
    ks <- seq(b[1], upper_k)
    cs <- cumsum(distrib_pdf(distrib, ks, theta))
    covered_k <- upper_k >= min(need_k, b[2])
    covered_p <- cs[length(cs)] >= need_p
    if ((covered_k && covered_p) || upper_k >= b[2] || size >= kmax) break
    size <- size * 2L
  }

  if (!(upper_k >= b[2]) && !(cs[length(cs)] >= need_p) && size >= kmax) {
    warning("Discrete fallback table reached its size limit (", kmax, " points) before covering the requested probability.")
  }

  list(ks = ks, cs = cs)
}

#' @title Default Numerical CDF for Discrete Distributions
#' @name distrib_cdf.discrete_distrib
#' @description
#' Fallback method: discrete distributions that do not implement an analytical CDF
#' get one by summing the pmf from the (finite) lower bound of the support up to
#' \eqn{\lfloor q \rfloor}.
#' @param distrib An object inheriting from class \code{"discrete_distrib"}.
#' @param q A numeric vector of quantiles.
#' @param theta A named list of parameters.
#' @param lower.tail Logical; if \code{TRUE} (default), probabilities are \eqn{P(Y \le q)}.
#' @param log.p Logical; if \code{TRUE}, probabilities are returned as logs.
#' @keywords internal
S7::method(distrib_cdf, discrete_distrib) <- function(distrib, q, theta, lower.tail = TRUE, log.p = FALSE) {
  b <- distrib@bounds
  all_params <- expand_params(c(list(.q = q), theta))
  rows <- transpose_params(all_params)

  prev_th <- NULL
  tab <- NULL

  res <- vapply(rows, function(r) {
    r <- as.list(r)
    th <- r[distrib@params]
    qi <- floor(r$.q)

    if (is.na(qi)) return(NA_real_)
    if (qi < b[1]) return(0)
    if (qi >= b[2]) return(1)

    if (!identical(th, prev_th) || tab$ks[length(tab$ks)] < qi) {
      tab <<- disc_cum_table(distrib, th, need_k = qi)
      prev_th <<- th
    }

    tab$cs[qi - b[1] + 1]
  }, numeric(1))

  res <- pmin(pmax(res, 0), 1)
  if (!lower.tail) res <- 1 - res
  if (log.p) log(res) else res
}

#' @title Default Numerical Quantile Function for Discrete Distributions
#' @name distrib_quantile.discrete_distrib
#' @description
#' Fallback method: discrete distributions that do not implement an analytical
#' quantile function get one by inverting the cumulative pmf table: the quantile is
#' the smallest support point \eqn{k} with \eqn{F(k) \ge p}.
#' @param distrib An object inheriting from class \code{"discrete_distrib"}.
#' @param p A numeric vector of probabilities.
#' @param theta A named list of parameters.
#' @param lower.tail Logical; if \code{TRUE} (default), probabilities are \eqn{P(Y \le p)}.
#' @param log.p Logical; if \code{TRUE}, probabilities are given as logs.
#' @keywords internal
S7::method(distrib_quantile, discrete_distrib) <- function(distrib, p, theta, lower.tail = TRUE, log.p = FALSE) {
  if (log.p) p <- exp(p)
  if (!lower.tail) p <- 1 - p

  b <- distrib@bounds

  # Scalar parameters, which is the case that matters for simulation: build the
  # cumulative table once and invert it for the whole vector in a single call.
  # The step function is inverted exactly, so the loop below buys nothing here --
  # it used to cost one R-level call per draw, which dominated everything else.
  if (all(lengths(theta) <= 1L)) {
    out <- rep(NaN, length(p))
    ok <- !is.na(p) & p >= 0 & p <= 1
    if (any(ok)) {
      tab <- disc_cum_table(distrib, theta, need_p = max(p[ok]))
      idx <- findInterval(p[ok], tab$cs, left.open = TRUE) + 1L
      idx <- pmin(idx, length(tab$ks))
      out[ok] <- tab$ks[idx]
      out[ok & p == 0] <- b[1]
      out[ok & p == 1] <- b[2]
    }
    return(out)
  }

  all_params <- expand_params(c(list(.p = p), theta))
  rows <- transpose_params(all_params)

  prev_th <- NULL
  tab <- NULL

  vapply(rows, function(r) {
    r <- as.list(r)
    th <- r[distrib@params]
    pi <- r$.p

    if (is.na(pi) || pi < 0 || pi > 1) return(NaN)
    if (pi == 0) return(b[1])
    if (pi == 1) return(b[2])

    if (!identical(th, prev_th) || tab$cs[length(tab$cs)] < min(pi, 1 - 1e-12)) {
      tab <<- disc_cum_table(distrib, th, need_p = pi)
      prev_th <<- th
    }

    # Smallest k with F(k) >= p
    idx <- findInterval(pi, tab$cs, left.open = TRUE) + 1L
    if (idx > length(tab$ks)) idx <- length(tab$ks)
    tab$ks[idx]
  }, numeric(1))
}

#' @title Default Numerical RNG for Discrete Distributions
#' @name distrib_rng.discrete_distrib
#' @description
#' Fallback method: discrete distributions that do not implement a native RNG
#' generate draws by inverse transform sampling on the cumulative pmf table.
#'
#' No rejection scheme is needed here, and none would help. Inverting the
#' cumulative mass function is exact, because the distribution function of a
#' lattice-valued variable is a step function: there is nothing to solve. The
#' table is built once per distinct parameter value and the whole sample is
#' located in it with a single binary search, which costs a fraction of a
#' microsecond per draw.
#' @param distrib An object inheriting from class \code{"discrete_distrib"}.
#' @param n Number of observations to generate.
#' @param theta A named list of parameters.
#' @keywords internal
S7::method(distrib_rng, discrete_distrib) <- function(distrib, n, theta) {
  distrib_quantile(distrib, stats::runif(n), theta)
}
