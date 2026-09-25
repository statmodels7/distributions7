# The expected information of the location-scale families that have none in
# elementary form -- skewnormal1, skewnormal2, skewt, pseudohuber -- is one
# quadrature over the standardized response per distinct shape.

azzalini_info <- function(a) {
  # Azzalini (1985), with a_k = E[Z^k {phi(aZ)/Phi(aZ)}^2]; an independent
  # route, sharing nothing with the family's own derivatives
  ak <- sapply(0:2, function(k) stats::integrate(function(z)
    2 * z^k * exp(stats::dnorm(z, log = TRUE) + 2 * stats::dnorm(a * z, log = TRUE) -
                  stats::pnorm(a * z, log.p = TRUE)),
    -Inf, Inf, rel.tol = 1e-12)$value)
  b <- sqrt(2 / pi)
  c(mu_mu = 1 + a^2 * ak[1], sigma_sigma = 2 + a^2 * ak[3], alpha_alpha = ak[3],
    mu_sigma = b * a * (1 + 2 * a^2) / (1 + a^2)^1.5 + a^2 * ak[2],
    mu_alpha = b / (1 + a^2)^1.5 - a * ak[2],
    sigma_alpha = -a * ak[3])
}

test_that("the skew normal reproduces Azzalini's information", {
  d <- skewnormal1_distrib()
  for (a in c(-5, 0.7, 3)) {
    e <- -unlist(distrib_expected_hessian(d, 0, list(mu = 0, sigma = 1, alpha = a)))
    ref <- azzalini_info(a)
    expect_lt(max(abs(e[names(ref)] - ref)) / max(abs(ref)), 1e-10)
  }
})

test_that("the rule is what carries the accuracy", {
  # a rule whose weights are 1e-6 out fails the check above
  good <- loc_scale_rule()
  bad <- good; bad$w <- bad$w * (1 + 1e-6)
  local_mocked_bindings(loc_scale_rule = function() bad)
  e <- -unlist(distrib_expected_hessian(skewnormal1_distrib(), 0,
                                        list(mu = 0, sigma = 1, alpha = 3)))
  ref <- azzalini_info(3)
  expect_gt(max(abs(e[names(ref)] - ref)) / max(abs(ref)), 1e-8)
})

test_that("the pseudo-Huber and the skew t agree with an adaptive quadrature", {
  skip_on_cran()
  cases <- list(
    list(pseudohuber_distrib(), list(mu = 0.3, sigma = 2, nu = 0.2)),
    list(pseudohuber_distrib(), list(mu = -1, sigma = 0.5, nu = 50)),
    list(skewt_distrib(), list(mu = 1, sigma = 2, alpha = 2, nu = 5)),
    list(skewt_distrib(), list(mu = 0, sigma = 0.3, alpha = -1, nu = 40)))
  for (cs in cases) {
    d <- cs[[1]]; th <- cs[[2]]
    e <- unlist(distrib_expected_hessian(d, 0, th))
    ref <- sapply(names(e), function(nm) {
      g <- function(x, i) {
        dm <- dim(x); xv <- as.vector(x)
        matrix(distrib_hessian(d, xv, th)[[nm]] * distrib_pdf(d, xv, th), dm[1], dm[2])
      }
      suppressWarnings(numericals7::quad_vec(g, -Inf, th$mu, rtol = 1e-10) +
                         numericals7::quad_vec(g, th$mu, Inf, rtol = 1e-10))
    })
    expect_lt(max(abs(e - ref)) / max(abs(ref)), 1e-9)
  }
})

test_that("the information ignores the observations and scales with sigma", {
  d <- skewt_distrib()
  th <- list(mu = 0.4, sigma = 1.3, alpha = 1.5, nu = 6)
  e1 <- distrib_expected_hessian(d, c(-3, 0, 7), th)
  e2 <- distrib_expected_hessian(d, c(100, 2, -1), th)
  expect_identical(e1, e2)
  e3 <- distrib_expected_hessian(d, 0, modifyList(th, list(mu = -8, sigma = 2.6)))
  k <- c(mu_mu = 2, sigma_sigma = 2, alpha_alpha = 0, nu_nu = 0, mu_sigma = 2,
         mu_alpha = 1, mu_nu = 1, sigma_alpha = 1, sigma_nu = 1, alpha_nu = 0)
  for (nm in names(k)) {
    expect_equal(e3[[nm]], e1[[nm]][1] / 2^k[[nm]], tolerance = 1e-12, info = nm)
  }
})

test_that("a shape shared by several observations is integrated once and read back", {
  d <- skewnormal1_distrib()
  a <- c(2, -1, 2, 0.5, -1)
  e <- distrib_expected_hessian(d, rep(0, 5), list(mu = 0, sigma = 1, alpha = a))
  for (i in seq_along(a)) {
    ei <- distrib_expected_hessian(d, 0, list(mu = 0, sigma = 1, alpha = a[i]))
    expect_equal(vapply(e, `[`, 0, i), unlist(ei), tolerance = 0)
  }
})

test_that("the pseudo-Huber's components odd in the location vanish", {
  d <- pseudohuber_distrib()
  th <- list(mu = 0.3, sigma = 1.7, nu = 2)
  e <- distrib_expected_hessian(d, 0, th)
  expect_lt(abs(e$mu_sigma), 1e-15)
  expect_lt(abs(e$mu_nu), 1e-15)
})

# the first derivative against a Richardson difference of the analytic
# expected information, the second against one of the analytic first
check_dexpected <- function(d, th, tol1, tol2) {
  P <- d@params; hn <- hess_names(P); n1 <- dexpected_names(P)
  x0 <- unlist(th)
  at <- function(x) stats::setNames(as.list(x), P)
  J1 <- numDeriv::jacobian(function(x) unlist(distrib_expected_hessian(d, 0, at(x))[hn]), x0)
  J2 <- numDeriv::jacobian(function(x) unlist(distrib_dexpected_hessian(d, 0, at(x))[n1]), x0)
  d1 <- distrib_dexpected_hessian(d, 0, th)
  d2 <- distrib_d2expected_hessian(d, 0, th)
  g1 <- 0
  for (i in seq_along(hn)) for (k in seq_along(P)) {
    g1 <- max(g1, abs(J1[i, k] - d1[[paste0(hn[i], "_", P[k])]]))
  }
  g2 <- 0
  pairs <- which(upper.tri(diag(length(P)), diag = TRUE), arr.ind = TRUE)
  for (r in seq_len(nrow(pairs))) for (c in seq_along(P)) for (k in seq_along(P)) {
    a <- pairs[r, 1L]; b <- pairs[r, 2L]
    i <- match(dexpected_key(P, a, b, c), n1)
    g2 <- max(g2, abs(J2[i, k] - d2[[d2expected_key(P, a, b, c, k)]]))
  }
  expect_lt(g1 / max(abs(J1)), tol1)
  expect_lt(g2 / max(abs(J2)), tol2)
}

test_that("the derivatives of the information agree with a difference of the order below", {
  skip_if_not_installed("numDeriv")
  skip_on_cran()
  check_dexpected(skewnormal1_distrib(), list(mu = 0.5, sigma = 1.4, alpha = 2.5), 1e-7, 1e-7)
  check_dexpected(pseudohuber_distrib(), list(mu = 0.5, sigma = 1.4, nu = 0.8), 1e-7, 1e-7)
  check_dexpected(skewnormal2_distrib(), list(mu = 0.5, sigma = 1.4, gamma1 = 0.4), 1e-7, 1e-7)
  # the skew t's components in nu come from stencils, so the second order
  # carries their accuracy
  check_dexpected(skewt_distrib(), list(mu = 0.5, sigma = 1.4, alpha = 1.5, nu = 8), 1e-5, 1e-3)
})

test_that("the link scale is carried by the shared rule", {
  d <- skewnormal1_distrib()
  th <- list(mu = 0.5, sigma = 1.4, alpha = 2.5)
  el <- distrib_expected_hessian(d, 0, th, scale = "link")
  d1 <- distrib_dexpected_hessian(d, 0, th, scale = "link")
  eta <- log(1.4); h <- 1e-5
  up <- distrib_expected_hessian(d, 0, modifyList(th, list(sigma = exp(eta + h))), scale = "link")
  dn <- distrib_expected_hessian(d, 0, modifyList(th, list(sigma = exp(eta - h))), scale = "link")
  for (nm in names(el)) {
    expect_equal(d1[[paste0(nm, "_sigma")]], (up[[nm]] - dn[[nm]]) / (2 * h),
                 tolerance = 1e-7, info = nm)
  }
})

test_that("the four families declare the quadrature, and a wrapper answers for its parent", {
  for (d in list(skewnormal1_distrib(), skewnormal2_distrib(), skewt_distrib(),
                 pseudohuber_distrib())) {
    expect_true(expected_hessian_by_quadrature(d), label = d@distrib_name)
    expect_true(expected_hessian_exact(d), label = d@distrib_name)
  }
  for (d in list(gaussian1_distrib(), gamma1_distrib(), pig1_distrib())) {
    expect_false(expected_hessian_by_quadrature(d), label = d@distrib_name)
  }
  expect_true(expected_hessian_by_quadrature(fixed(skewt_distrib(), nu = 6)))
  expect_false(expected_hessian_by_quadrature(fixed(gaussian1_distrib(), mu = 0)))
})
