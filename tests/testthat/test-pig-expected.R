# The exact expected information of the two Poisson-inverse Gaussians and its
# derivatives, from one pass over the support per observation.

pig_scale_gap <- function(a, b) max(abs(a - b)) / max(abs(b))

# -E[l_a l_b] summed over a support long enough to hold the whole mass, with
# the observed kernels: the second Bartlett identity, sharing with the kernel
# only the row it evaluates at each count
pig_opg_sum <- function(d, th, Y) {
  y <- 0:Y
  f <- distrib_pdf(d, y, th)
  g <- distrib_gradient(d, y, th)
  out <- c(-sum(f * g[[1]]^2), -sum(f * g[[2]]^2), -sum(f * g[[1]] * g[[2]]))
  stats::setNames(out, hess_names(d@params))
}

test_that("the expected information is the summed outer product of the scores", {
  for (cs in list(c(0.5, 0.3), c(3, 0.8), c(3, 3), c(30, 0.3), c(2, 1e-4))) {
    al <- sqrt(1 + 2 * cs[1] * cs[2]) / cs[2]
    for (d in list(pig1_distrib(), pig2_distrib())) {
      th <- stats::setNames(list(cs[1], if (d@params[2] == "sigma") cs[2] else al),
                            d@params)
      got <- unlist(distrib_expected_hessian(d, 1, th))
      ref <- pig_opg_sum(d, th, 3000)
      expect_lt(pig_scale_gap(got[names(ref)], ref), 1e-12)
    }
  }
})

test_that("the derivatives agree with one central difference of the order below", {
  # compared on the scale of the whole vector: a component that happens to
  # be small is not a reason for the difference to be relatively large, and
  # the central difference's own floor at this step is about 1e-8
  fd <- function(fun, d, th, j, h = 1e-4) {
    p <- th; q <- th; hh <- h * th[[j]]
    p[[j]] <- p[[j]] + hh; q[[j]] <- q[[j]] - hh
    (unlist(fun(d, 1, p)) - unlist(fun(d, 1, q))) / (2 * hh)
  }
  for (d in list(pig1_distrib(), pig2_distrib())) {
    P <- d@params
    th <- stats::setNames(list(3, if (P[2] == "sigma") 0.8 else 2), P)
    hn <- hess_names(P)
    d1 <- unlist(distrib_dexpected_hessian(d, 1, th))
    d2 <- unlist(distrib_d2expected_hessian(d, 1, th))
    F1 <- lapply(1:2, function(j) fd(distrib_expected_hessian, d, th, j))
    G1 <- lapply(1:2, function(j) fd(distrib_dexpected_hessian, d, th, j))
    a1 <- b1 <- a2 <- b2 <- numeric(0)
    for (ab in hn) for (k in 1:2) {
      a1 <- c(a1, d1[[paste(ab, P[k], sep = "_")]]); b1 <- c(b1, F1[[k]][[ab]])
    }
    for (ab in hn) {
      a2 <- c(a2, d2[[paste(ab, P[1], P[1], sep = "_")]],
              d2[[paste(ab, P[2], P[2], sep = "_")]],
              d2[[paste(ab, P[1], P[2], sep = "_")]])
      b2 <- c(b2, G1[[1]][[paste(ab, P[1], sep = "_")]],
              G1[[2]][[paste(ab, P[2], sep = "_")]],
              G1[[2]][[paste(ab, P[1], sep = "_")]])
    }
    expect_lt(pig_scale_gap(a1, b1), 1e-7, label = d@distrib_name)
    expect_lt(pig_scale_gap(a2, b2), 1e-7, label = d@distrib_name)
    # and a 1e-4 error in one component is seen
    a1[2] <- a1[2] * (1 + 1e-4) + 1e-4 * max(abs(a1))
    expect_gt(pig_scale_gap(a1, b1), 1e-6)
  }
})

test_that("pig1 reaches the Poisson limit, and the direct form does not", {
  # With mu = 2 the three quantities settle on -2 + 10 sigma, 10 and -52 as
  # the dispersion falls; the first is the Poisson information -mu^2/2 in
  # sigma, and all three are reached monotonically from sigma = 1e-2 down,
  # measured, which a cancelling form cannot do.
  d <- pig1_distrib()
  for (s in c(1e-6, 1e-9, 1e-12)) {
    th <- list(mu = 2, sigma = s)
    expect_equal(distrib_expected_hessian(d, 1, th)$sigma_sigma, -2 + 10 * s,
                 tolerance = 1e-9)
    expect_equal(distrib_dexpected_hessian(d, 1, th)$sigma_sigma_sigma, 10,
                 tolerance = 1e-5)
    expect_equal(distrib_d2expected_hessian(d, 1, th)$sigma_sigma_sigma_sigma,
                 -52, tolerance = 1e-5)
  }
  # the direct form E[l_ab] of the same sum, in pig2's coordinates, has lost
  # every digit by alpha = 1e6, where -E[l_a^2] still reads -mu^2/2 alpha^-4
  d2 <- pig2_distrib(); a <- 1e6
  y <- 0:60; th <- list(mu = 2, alpha = a)
  f <- distrib_pdf(d2, y, th)
  direct <- sum(f * distrib_hessian(d2, y, th)$alpha_alpha) * a^4
  kernel <- distrib_expected_hessian(d2, 1, th)$alpha_alpha * a^4
  expect_equal(kernel, -2, tolerance = 1e-5)
  expect_gt(abs(direct - (-2)), 1)
})

test_that("pig1's observed derivatives in sigma do not cancel at a small dispersion", {
  # Against the per-count chain from pig2, l1_sigma = l2_alpha alpha_sigma,
  # a product and so free of cancellation. The composition through
  # psi(alpha) = -alpha + log S that this replaced lost digits as sigma^-2:
  # 5e-10 at 1e-4, 4.5e-2 at 1e-8, nothing at 1e-9.
  m <- 2; y <- 0:12
  ae <- quote(sqrt(1 + 2 * sg * m) / sg)
  for (s in c(1e-4, 1e-8, 1e-11)) {
    e <- list(m = m, sg = s)
    al <- eval(ae, e); as <- eval(stats::D(ae, "sg"), e)
    g1 <- distrib_gradient(pig1_distrib(), y, list(mu = m, sigma = s))$sigma
    g2 <- distrib_gradient(pig2_distrib(), y, list(mu = m, alpha = al))$alpha
    expect_lt(pig_scale_gap(g1, g2 * as), 1e-13)
  }
  # and the limit itself: l_sigma -> ((y - mu)^2 - y) / 2
  g <- distrib_gradient(pig1_distrib(), y, list(mu = m, sigma = 1e-12))$sigma
  expect_equal(g, ((y - m)^2 - y) / 2, tolerance = 1e-10)
})

test_that("the kernel does not depend on the count of threads", {
  d <- pig2_distrib()
  th <- list(mu = c(3, 30, 0.5), alpha = c(2, 0.4, 8))
  for (fn in list(distrib_expected_hessian, distrib_d2expected_hessian)) {
    expect_identical(fn(d, rep(1, 3), th, threads = 1L),
                     fn(d, rep(1, 3), th, threads = 2L))
  }
  # a parameter outside its domain reads NA rather than a partial sum
  expect_true(all(is.na(unlist(pig2_expected_cpp(1, -1, 2, 0L)))))
})

test_that("the two parametrizations are one law: pig1 is pig2 through the map", {
  # E1 = J' E2 J with J the Jacobian of (mu, alpha(mu, sigma)); the two
  # kernels run different recurrences (in w and in alpha) and different rows.
  for (cs in list(c(3, 0.8), c(0.5, 3), c(20, 0.05))) {
    m <- cs[1]; s <- cs[2]
    ae <- quote(sqrt(1 + 2 * sg * m) / sg); e <- list(m = m, sg = s)
    al <- eval(ae, e)
    J <- rbind(c(1, 0), c(eval(stats::D(ae, "m"), e), eval(stats::D(ae, "sg"), e)))
    E2 <- distrib_expected_hessian(pig2_distrib(), 1, list(mu = m, alpha = al))
    M <- matrix(c(E2$mu_mu, E2$mu_alpha, E2$mu_alpha, E2$alpha_alpha), 2)
    N <- t(J) %*% M %*% J
    E1 <- distrib_expected_hessian(pig1_distrib(), 1, list(mu = m, sigma = s))
    expect_lt(pig_scale_gap(c(E1$mu_mu, E1$sigma_sigma, E1$mu_sigma),
                            c(N[1, 1], N[2, 2], N[1, 2])), 1e-12)
  }
})
