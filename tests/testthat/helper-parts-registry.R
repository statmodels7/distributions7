# The twin comparison of the parts route of the C registry (pt_parts.h),
# shared by test-parts-registry*.R, which split the loop over families so
# that testthat's parallel workers run them at once. The families are those of
# helper-joint-registry.R, with their constants.

# the full P^k array, observations first, from a generic's named components
full_array <- function(comp, params, k, n) {
  P <- length(params)
  out <- array(NA_real_, dim = c(n, rep(P, k)))
  idx <- as.matrix(expand.grid(rep(list(seq_len(P)), k)))
  perms_of <- function(t) {
    if (length(t) <= 1) return(list(t))
    out <- list()
    for (i in seq_along(t)) for (r in perms_of(t[-i])) out[[length(out) + 1L]] <- c(t[i], r)
    unique(out)
  }
  for (r in seq_len(nrow(idx))) {
    t <- idx[r, ]
    key <- NULL
    for (pp in perms_of(t)) {
      nm <- paste(params[pp], collapse = "_")
      if (!is.null(comp[[nm]])) { key <- nm; break }
    }
    if (is.null(key)) stop("no component for ", paste(params[t], collapse = "_"))
    out[cbind(seq_len(n), matrix(rep(t, each = n), nrow = n))] <- rep_len(comp[[key]], n)
  }
  out
}

# points drawn on the link scale, with the two extremes of the response that
# the support admits, as joint_twin() takes them (a far count of 1e6)
parts_points <- function(d, seed, eta_range) {
  set.seed(seed)
  n <- 60
  th <- lapply(d@params, function(p)
    linkfunctions7::linkinv(d@link_params[[p]],
                            stats::runif(n, eta_range[1], eta_range[2])))
  names(th) <- d@params
  y <- vapply(seq_len(n), function(i) distrib_rng(d, 1, lapply(th, `[`, i)), 0)
  hi <- d@bounds[2]
  lo <- d@bounds[1]
  discrete <- S7::S7_inherits(d, discrete_distrib)
  far_y <- c(if (is.finite(hi)) hi * (1 - 1e-12) else if (discrete) 1e6 else 1e40,
             if (is.finite(lo)) lo + max(abs(lo), 1) * 1e-12 else -1e40)
  if (discrete) far_y <- round(far_y)
  y[1:2] <- far_y
  list(th = th, y = y, n = n)
}

# the ranges each family is read on: the skew t to nu = exp(1.5), as in
# test-joint-registry.R; the negative binomials also on a wide one, which
# reaches their sums over the count (theta > 100 max(y, 1))
parts_ranges <- function(cls) {
  if (cls == "SkewTDistrib") return(list(c(-1.5, 1.5)))
  r <- list(c(-1.5, 1.5), c(-4, 4))
  if (cls %in% c("NegBin1Distrib", "NegBin2Distrib")) r <- c(r, list(c(-6, 8)))
  r
}

# the gradient, the Hessian and the expected Hessian against the joint
# entries, the orders three and four against distrib_deriv3()/deriv4()
parts_twin <- function(cls, seed) {
  d <- joint_families[[cls]]()
  route <- distrib_scalar_route(d)
  for (rg in parts_ranges(cls)) {
    pt <- parts_points(d, seed, rg)
    M <- do.call(cbind, c(pt$th, lapply(route$constants, rep, length.out = pt$n)))
    pr <- d7_parts_probe(route$name, pt$y, M)
    jr <- d7_joint_probe(route$name, pt$y, M)
    expect_identical(pr$covered, 4L, label = paste(cls, "coverage"))
    lab <- paste(cls, paste(rg, collapse = ".."))
    expect_agrees_on_scale(pr$grad, jr$grad, label = paste(lab, "gradient"))
    expect_agrees_on_scale(pr$hess, jr$hess, label = paste(lab, "Hessian"))
    expect_agrees_on_scale(pr$info, jr$info, label = paste(lab, "expected"))
    d3 <- full_array(distrib_deriv3(d, pt$y, pt$th), d@params, 3, pt$n)
    d4 <- full_array(distrib_deriv4(d, pt$y, pt$th), d@params, 4, pt$n)
    expect_agrees_on_scale(as.vector(pr$d3), as.vector(d3),
                           label = paste(lab, "third derivative"))
    expect_agrees_on_scale(as.vector(pr$d4), as.vector(d4),
                           label = paste(lab, "fourth derivative"))
  }
}
