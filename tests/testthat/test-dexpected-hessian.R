# The derivative of the expected information.

test_that("the gaussian's components are the hand-written ones", {
  # E[l_mumu] = -1/s^2, E[l_musigma] = 0, E[l_sigmasigma] = -2/s^2, so every
  # component is known outright. This is the one reference that shares no
  # arithmetic at all with the route under test.
  s <- 1.7
  g <- distrib_dexpected_hessian(gaussian1_distrib(), 0,
                                 list(mu = 0.4, sigma = s),
                                 scale = "parameter")
  expect_equal(as.numeric(g[["mu_mu_sigma"]]), 2 / s^3, tolerance = 1e-8)
  expect_equal(as.numeric(g[["sigma_sigma_sigma"]]), 4 / s^3, tolerance = 1e-8)
  for (nm in c("mu_mu_mu", "sigma_sigma_mu", "mu_sigma_mu", "mu_sigma_sigma")) {
    expect_lt(abs(as.numeric(g[[nm]])), 1e-8)
  }
})

test_that("the components obey E[l_abc] + E[l_ab l_c]", {
  # The identity the closed forms would be derived from, and a route that
  # shares nothing with the stencil: distrib_deriv3(expected = TRUE) for the
  # first term and one quadrature per component for the second. The beta is
  # the family to ask it of, every one of whose six components is non-zero --
  # a family whose components mostly vanish would satisfy this trivially.
  d <- beta1_distrib()
  th <- list(mu = 0.42, phi = 5.5)
  params <- d@params
  a <- distrib_dexpected_hessian(d, 0, th, scale = "parameter")
  e3 <- distrib_deriv3(d, 0, th, expected = TRUE, scale = "parameter")
  for (i in seq_along(params)) {
    for (j in seq_along(params)) {
      for (k in seq_along(params)) {
        nm <- dexpected_key(params, i, j, k)
        mix <- expectation(d, function(y, theta, ab, cc) {
          distrib_hessian(d, y, theta)[[ab[1]]] *
            distrib_gradient(d, y, theta)[[cc[1]]]
        }, theta = th, ab = sub(paste0("_", params[k], "$"), "", nm),
        cc = params[k])
        want <- e3[[paste(params[sort(c(i, j, k))], collapse = "_")]] + mix
        expect_equal(as.numeric(a[[nm]]), as.numeric(want), tolerance = 1e-6)
      }
    }
  }
})

test_that("the link scale is the expected information's own", {
  # It is differenced along the free scale of the parameter differentiated in
  # and read off distrib_expected_hessian(scale = "link") at each point, so
  # nothing here can disagree with the chain rule that generic applies. At an
  # IDENTITY link the two scales must coincide, which is what pins that.
  d <- gaussian1_distrib(link_mu = linkfunctions7::identity_link())
  th <- list(mu = 0.3, sigma = 1.4)
  p <- distrib_dexpected_hessian(d, 0, th, scale = "parameter")
  l <- distrib_dexpected_hessian(d, 0, th, scale = "link")
  expect_equal(as.numeric(p[["mu_mu_mu"]]), as.numeric(l[["mu_mu_mu"]]),
               tolerance = 1e-8)
})

test_that("a family that approximates its expected information is refused", {
  # Not an accuracy judgement: measured at 100 observations these six cost
  # 1880 to 147300 ms against a median of 0.183 ms for the thirty-four that
  # write the expectation out, so 2p of those calls per evaluation is
  # unusable rather than merely slow.
  for (d in list(skewnormal1_distrib(), skewnormal2_distrib(),
                 pseudohuber_distrib(), skewt_distrib())) {
    th <- generate_random_theta(d)
    expect_error(distrib_dexpected_hessian(d, distrib_rng(d, 1L, th), th),
                 "approximates its expected information")
  }
})

test_that("the exactness predicate follows the arithmetic, not the owner", {
  # skewnormal2 registers its own method and pseudohuber registers one that
  # patches two components of the fallback, so reading the owning class said
  # "written out" about a quadrature in both. The consequences were live:
  # fit_distrib() rejected a legitimate fisher_scoring(approx = ) on them with
  # a message that was untrue.
  expect_false(has_exact_expected_hessian(skewnormal2_distrib()))
  expect_false(has_exact_expected_hessian(pseudohuber_distrib()))
  expect_false(has_exact_expected_hessian(skewnormal1_distrib()))
  # and the ones that do write it out still say so, including a family
  # reached through a reparametrization whose parent is exact
  expect_true(has_exact_expected_hessian(gaussian1_distrib()))
  expect_true(has_exact_expected_hessian(negbin2_distrib()))
  expect_true(has_exact_expected_hessian(weibull3_distrib()))
  # the refusal that had been turned into an error is available again
  expect_silent(fit_distrib(pseudohuber_distrib(),
                            distrib_rng(pseudohuber_distrib(), 5L,
                                        list(mu = 0, sigma = 1, nu = 1)),
                            method = fisher_scoring(approx = "mc", nsim = 20L),
                            n_start = 1L))
})

test_that("the keys are built and enumerated by the same rule", {
  params <- c("mu", "sigma")
  nm <- dexpected_names(params)
  expect_length(nm, length(hess_names(params)) * length(params))
  built <- unlist(lapply(seq_along(params), function(a)
    lapply(seq_along(params), function(b)
      lapply(seq_along(params), function(k) dexpected_key(params, a, b, k)))))
  expect_true(all(built %in% nm))
  # the component is symmetric in (a, b) and NOT in the parameter it is
  # differentiated in, which is what makes its keying differ from deriv_names'
  expect_identical(dexpected_key(params, 1, 2, 2),
                   dexpected_key(params, 2, 1, 2))
  expect_false(identical(dexpected_key(params, 1, 1, 2),
                         dexpected_key(params, 1, 2, 1)))
  d <- gamma1_distrib()
  expect_named(distrib_dexpected_hessian(d, 0, list(mu = 2, phi = 1.5)),
               dexpected_names(d@params), ignore.order = TRUE)
})

# The five families whose derivatives of the expected information are written
# out in a compiled kernel, at a point where every component that is not zero
# by construction is of order one.
analytic_dexpected_cases <- function() list(
  list(d = gaussian1_distrib(), th = list(mu = 0.4, sigma = 1.7)),
  list(d = poisson_distrib(),   th = list(mu = 3.2)),
  list(d = gamma1_distrib(),    th = list(mu = 2.1, phi = 0.37)),
  list(d = negbin2_distrib(),   th = list(mu = 4.3, theta = 2.6)),
  list(d = beta1_distrib(),     th = list(mu = 0.42, phi = 5.5))
)

# |a - b| against the size of the whole list, so that a component that is zero
# by construction is compared absolutely: the quadrature leaves 1e-10 there.
expect_close_list <- function(a, b, tol) {
  sc <- max(1, vapply(b, function(v) max(abs(v)), numeric(1)))
  gap <- vapply(names(b), function(k) max(abs(a[[k]] - b[[k]])), numeric(1))
  expect_lt(max(gap) / sc, tol)
}

test_that("the analytic first derivatives obey E[l_abc] + E[l_ab l_c]", {
  # The right-hand side by expectation() -- a quadrature or an exact sum over
  # the support -- which shares no arithmetic with the kernels.
  for (cs in analytic_dexpected_cases()) {
    d <- cs$d; th <- cs$th; P <- d@params; np <- length(P)
    a <- distrib_dexpected_hessian(d, 0, th)
    expect_identical(names(a), dexpected_names(P))
    want <- list()
    for (i in 1:np) for (j in i:np) for (k in 1:np) {
      hab <- hess_pair_name(P, i, j)
      want[[dexpected_key(P, i, j, k)]] <- expectation(d, function(y, theta)
        distrib_deriv3(d, y, theta)[[paste(P[sort(c(i, j, k))], collapse = "_")]] +
          distrib_hessian(d, y, theta)[[hab]] *
          distrib_gradient(d, y, theta)[[P[k]]], theta = th)
    }
    expect_close_list(a, want, 1e-7)
  }
})

test_that("the analytic second derivatives obey the five-moment identity", {
  for (cs in analytic_dexpected_cases()) {
    d <- cs$d; th <- cs$th; P <- d@params; np <- length(P)
    a <- distrib_d2expected_hessian(d, 0, th)
    expect_identical(names(a), d2expected_names(P))
    srt <- function(ix) paste(P[sort(ix)], collapse = "_")
    want <- list()
    for (i in 1:np) for (j in i:np) for (k in 1:np) for (l in k:np) {
      hab <- hess_pair_name(P, i, j)
      want[[d2expected_key(P, i, j, k, l)]] <- expectation(d, function(y, theta) {
        g <- distrib_gradient(d, y, theta)
        h <- distrib_hessian(d, y, theta)
        d3 <- distrib_deriv3(d, y, theta)
        distrib_deriv4(d, y, theta)[[srt(c(i, j, k, l))]] +
          d3[[srt(c(i, j, l))]] * g[[P[k]]] + d3[[srt(c(i, j, k))]] * g[[P[l]]] +
          h[[hab]] * h[[hess_pair_name(P, k, l)]] +
          h[[hab]] * g[[P[k]]] * g[[P[l]]]
      }, theta = th)
    }
    expect_close_list(a, want, 1e-7)
  }
})

test_that("the analytic first derivatives agree with the stencil they replace", {
  for (cs in analytic_dexpected_cases()) {
    for (sc in c("parameter", "link")) {
      a <- distrib_dexpected_hessian(cs$d, 0, cs$th, scale = sc)
      n <- numerical_dexpected_hessian(cs$d, 0, cs$th, sc)
      expect_close_list(a, n, 1e-6)
    }
  }
})

test_that("the link-scale second derivatives are the derivative of the first", {
  # One central difference of the ANALYTIC link-scale first derivative along
  # the free coordinate, which is a reference and not a route: it shares the
  # kernels and nothing of dexpected_link()'s second-order Leibniz terms.
  for (cs in analytic_dexpected_cases()) {
    d <- cs$d; th <- cs$th; P <- d@params; np <- length(P)
    a <- distrib_d2expected_hessian(d, 0, th, scale = "link")
    want <- list()
    for (k in 1:np) {
      lk <- d@link_params[[P[k]]]
      e0 <- linkfunctions7::linkfun(lk, th[[k]])
      h <- 1e-5 * max(1, abs(e0))
      tp <- tm <- th
      tp[[k]] <- linkfunctions7::linkinv(lk, e0 + h)
      tm[[k]] <- linkfunctions7::linkinv(lk, e0 - h)
      up <- distrib_dexpected_hessian(d, 0, tp, scale = "link")
      dn <- distrib_dexpected_hessian(d, 0, tm, scale = "link")
      for (i in 1:np) for (j in i:np) for (l in k:np) {
        k1 <- dexpected_key(P, i, j, l)
        want[[d2expected_key(P, i, j, k, l)]] <- (up[[k1]] - dn[[k1]]) / (2 * h)
      }
    }
    expect_close_list(a, want, 1e-6)
  }
})

test_that("the analytic kernels are identical at any thread count", {
  n <- 600L
  for (cs in analytic_dexpected_cases()) {
    th <- lapply(cs$th, function(v) v * seq(0.9, 1.1, length.out = n))
    if (!is.null(th$mu) && S7::S7_inherits(cs$d, Beta1Distrib)) {
      th$mu <- seq(0.3, 0.6, length.out = n)
    }
    y <- rep(0, n)
    for (fn in list(distrib_dexpected_hessian, distrib_d2expected_hessian)) {
      expect_identical(fn(cs$d, y, th, threads = 1L), fn(cs$d, y, th, threads = 2L))
    }
  }
})

test_that("a family without a second derivative of its expected information is refused", {
  d <- weibull1_distrib()
  expect_error(distrib_d2expected_hessian(d, 1, list(mu = 1, sigma = 2)),
               "no analytic second derivative")
})

test_that("the second-derivative keys are symmetric in each pair", {
  P <- c("mu", "sigma")
  expect_length(d2expected_names(P), length(hess_names(P))^2)
  expect_identical(d2expected_key(P, 1, 2, 2, 1), d2expected_key(P, 2, 1, 1, 2))
  expect_false(identical(d2expected_key(P, 1, 1, 2, 2), d2expected_key(P, 2, 2, 1, 1)))
  expect_true(all(d2expected_key(P, 2, 1, 1, 2) %in% d2expected_names(P)))
})
