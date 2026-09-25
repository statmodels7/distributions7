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
  # Not an accuracy judgement: an expected information approximated per
  # observation costs seconds at 100 observations against a median of 0.183 ms
  # for a family that writes it out, so 2p of those calls per evaluation is
  # unusable rather than merely slow.
  for (d in list(pig1_distrib(), pig2_distrib())) {
    th <- generate_random_theta(d)
    expect_error(distrib_dexpected_hessian(d, distrib_rng(d, 1L, th), th),
                 "approximates its expected information")
  }
})

test_that("the exactness predicate follows the arithmetic, not the owner", {
  # skewnormal2 chains onto skewnormal1 and answers for it
  expect_identical(has_exact_expected_hessian(skewnormal2_distrib()),
                   has_exact_expected_hessian(skewnormal1_distrib()))
  expect_true(has_exact_expected_hessian(skewnormal2_distrib()))
  expect_true(has_exact_expected_hessian(pseudohuber_distrib()))
  expect_false(has_exact_expected_hessian(pig1_distrib()))
  # and the ones that do write it out still say so, including a family
  # reached through a reparametrization whose parent is exact
  expect_true(has_exact_expected_hessian(gaussian1_distrib()))
  expect_true(has_exact_expected_hessian(negbin2_distrib()))
  expect_true(has_exact_expected_hessian(weibull3_distrib()))
  # a strategy for the approximation is accepted where there is one to make,
  # and refused where the family computes its information exactly
  set.seed(1)
  y <- distrib_rng(pig1_distrib(), 30L, list(mu = 3, sigma = 0.5))
  expect_no_error(fit_distrib(pig1_distrib(), y,
                              method = fisher_scoring(approx = "mc", nsim = 20L),
                              n_start = 1L))
  expect_error(fit_distrib(pseudohuber_distrib(),
                           distrib_rng(pseudohuber_distrib(), 5L,
                                       list(mu = 0, sigma = 1, nu = 1)),
                           method = fisher_scoring(approx = "mc", nsim = 20L),
                           n_start = 1L),
               "computes its expected information exactly")
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
  list(d = beta1_distrib(),     th = list(mu = 0.42, phi = 5.5)),
  list(d = bernoulli_distrib(), th = list(mu = 0.37)),
  list(d = binomial_distrib(size = 10), th = list(mu = 0.37)),
  list(d = exponential_distrib(), th = list(mu = 1.7)),
  list(d = geometric_distrib(), th = list(mu = 2.3)),
  list(d = chisq_distrib(),     th = list(mu = 3.4)),
  list(d = cauchy_distrib(),    th = list(mu = 0.3, sigma = 1.6)),
  list(d = logistic_distrib(),  th = list(mu = 0.3, sigma = 1.6)),
  # the quadrature of the moment identity fails in the Gumbel's left tail,
  # where the observed derivatives carry exp(-z) against a vanishing density;
  # Richardson covers it
  list(d = gumbel_distrib(),    th = list(mu = 0.3, sigma = 1.6), quad = FALSE),
  list(d = gaussian2_distrib(), th = list(mu = 0.3, sigma2 = 1.8)),
  list(d = gaussian3_distrib(), th = list(mu = 0.3, tau = 0.7)),
  list(d = lognormal1_distrib(), th = list(mu = 0.3, sigma2 = 0.6)),
  list(d = invgauss1_distrib(), th = list(mu = 1.4, phi = 0.8)),
  list(d = invgauss2_distrib(), th = list(mu = 1.4, lambda = 2.2)),
  list(d = gamma2_distrib(),    th = list(mu = 2.1, sigma2 = 1.3)),
  list(d = beta2_distrib(),     th = list(alpha = 2.3, beta = 3.1)),
  list(d = weibull1_distrib(),  th = list(mu = 1.7, sigma = 2.2)),
  list(d = student_t1_distrib(), th = list(mu = 0.3, sigma = 1.4, nu = 6.5)),
  # above nu = 30 the derivatives of E[l_nu_nu] come from the series
  list(d = student_t1_distrib(), th = list(mu = 0.3, sigma = 1.4, nu = 45)),
  list(d = gengamma1_distrib(), th = list(a = 1.7, d = 2.3, p = 1.6)),
  # the reparametrized families, through reparam_dexpected()
  list(d = student_t2_distrib(), th = list(mu = -1.99, sigma = 4.86, nu = 4.37)),
  list(d = gengamma2_distrib(), th = list(mean = 1.09, d = 4.86, p = 2.37)),
  list(d = lognormal2_distrib(), th = list(mean = 1.09, var = 4.86)),
  list(d = weibull3_distrib(),  th = list(mean = 1.09, sigma = 4.86)),
  list(d = vonmises1_distrib(), th = list(mu = 0.4, kappa = 2.3)),
  # past kappa = 20 the Bessel ratio's derivatives come from its series
  list(d = vonmises1_distrib(), th = list(mu = 0.4, kappa = 40)),
  list(d = vonmises2_distrib(), th = list(mu = 0.4, rho = 0.6)),
  list(d = vonmises2_distrib(), th = list(mu = 0.4, rho = 0.97)),
  list(d = negbin1_distrib(),   th = list(mu = 3.4, theta = 0.7)),
  list(d = betabinom1_distrib(size = 10), th = list(mu = 0.35, sigma = 0.4)),
  list(d = betabinom2_distrib(size = 10), th = list(alpha = 1.3, beta = 2.4)),
  list(d = gpd_distrib(),       th = list(sigma = 1.5, xi = 0.2)),
  # a bounded support: the quadrature of the identities does not converge
  list(d = gpd_distrib(),       th = list(sigma = 1.5, xi = -0.2), quad = FALSE)
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
    if (isFALSE(cs$quad)) next
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
    if (isFALSE(cs$quad)) next
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
    th <- lapply(cs$th, function(v) v * seq(0.97, 1.03, length.out = n))
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
  d <- pig1_distrib()
  expect_error(distrib_d2expected_hessian(d, 1, list(mu = 2, sigma = 0.5)),
               "no analytic second derivative")
})

test_that("the second-derivative keys are symmetric in each pair", {
  P <- c("mu", "sigma")
  expect_length(d2expected_names(P), length(hess_names(P))^2)
  expect_identical(d2expected_key(P, 1, 2, 2, 1), d2expected_key(P, 2, 1, 1, 2))
  expect_false(identical(d2expected_key(P, 1, 1, 2, 2), d2expected_key(P, 2, 2, 1, 1)))
  expect_true(all(d2expected_key(P, 2, 1, 1, 2) %in% d2expected_names(P)))
})

# Richardson on the parameter scale: the first order against one difference of
# the family's own analytic expected information, the second against one
# difference of the analytic first. Neither nests a difference.
dexpected_richardson <- function(d, th) {
  P <- d@params; np <- length(P); hn <- hess_names(P); n1 <- dexpected_names(P)
  x0 <- unlist(th)
  at <- function(x) stats::setNames(as.list(x), P)
  J1 <- numDeriv::jacobian(function(x)
    unlist(distrib_expected_hessian(d, 0, at(x))[hn]), x0)
  J2 <- numDeriv::jacobian(function(x)
    unlist(distrib_dexpected_hessian(d, 0, at(x))[n1]), x0)
  r1 <- list(); r2 <- list()
  for (i in seq_along(hn)) for (k in seq_len(np)) r1[[paste0(hn[i], "_", P[k])]] <- J1[i, k]
  pairs <- lapply(seq_along(hn), function(i) {
    w <- match(strsplit(hn[i], "_")[[1]], P)
    if (anyNA(w)) stop("parameter names with underscores are not covered here")
    w
  })
  for (ab in pairs) for (cd in pairs) {
    k1 <- dexpected_key(P, ab[1], ab[2], cd[2])
    r2[[d2expected_key(P, ab[1], ab[2], cd[1], cd[2])]] <- J2[match(k1, n1), cd[1]]
  }
  list(r1 = r1, r2 = r2)
}

test_that("the analytic derivatives agree with Richardson on the analytic order below", {
  skip_if_not_installed("numDeriv")
  for (cs in analytic_dexpected_cases()) {
    r <- dexpected_richardson(cs$d, cs$th)
    expect_close_list(distrib_dexpected_hessian(cs$d, 0, cs$th), r$r1, 1e-8)
    expect_close_list(distrib_d2expected_hessian(cs$d, 0, cs$th), r$r2, 1e-8)
  }
})

test_that("a component 5 per cent out fails the Richardson comparison", {
  # The negative control of the test above: the gamma2 kernel with its
  # largest second-order component inflated by 5 per cent.
  skip_if_not_installed("numDeriv")
  d <- gamma2_distrib(); th <- list(mu = 2.1, sigma2 = 1.3)
  r <- dexpected_richardson(d, th)
  good <- distrib_d2expected_hessian(d, 0, th)
  real <- gamma2_dexpected_cpp
  local_mocked_bindings(gamma2_dexpected_cpp = function(y, mu, sigma2, order, threads = 1) {
    out <- real(y, mu, sigma2, order, threads)
    if (order == 2L) {
      k <- names(out)[which.max(vapply(out, function(v) max(abs(v)), 0))]
      out[[k]] <- 1.05 * out[[k]]
    }
    out
  })
  bad <- distrib_d2expected_hessian(d, 0, th)
  gap <- function(a) max(vapply(names(r$r2), function(k) max(abs(a[[k]] - r$r2[[k]])), 0))
  expect_lt(gap(good), 1e-8)
  expect_gt(gap(bad), 1e-3)
})

test_that("the gamma by its variance keeps its digits at a large shape", {
  # With a = mu^2 / sigma2 the function g = a psi'(a) - 1 tends to 1/(2a), and
  # the components reach the derivatives of the limiting information, among
  # them d E_{mu sigma2} / d mu -> -1 / (mu^2 sigma2), at the rate 1/a. The
  # direct form a psi'(a) - 1 would lose every digit of g by a = 1e8.
  v <- 0.5
  gaps <- vapply(c(1e4, 1e6, 1e8), function(a) {
    m <- sqrt(a * v)
    g <- distrib_dexpected_hessian(gamma2_distrib(), 0, list(mu = m, sigma2 = v))
    abs(g[["mu_sigma2_mu"]] * (-m^2 * v) - 1)
  }, 0)
  expect_equal(gaps, c(1e-4, 1e-6, 1e-8), tolerance = 1e-2)
})

test_that("the Student t's nu derivatives keep their digits at a large nu", {
  # E[l_nu_nu] = -7/(2 nu^4) + ..., so its derivatives tend to 14/nu^5 and
  # -70/nu^6; the direct form cancels from order nu^-2 and would read noise
  # long before nu = 1e6.
  d <- student_t1_distrib()
  for (v in c(1e3, 1e6, 1e8)) {
    th <- list(mu = 0, sigma = 1, nu = v)
    a <- distrib_dexpected_hessian(d, 0, th)[["nu_nu_nu"]]
    b <- distrib_d2expected_hessian(d, 0, th)[["nu_nu_nu_nu"]]
    expect_equal(a * v^5 / 14, 1, tolerance = 10 / v)
    expect_equal(b * v^6 / -70, 1, tolerance = 10 / v)
  }
  # and the link scale is finite at a nu far past double precision's square
  # root, which is where a chain dividing by nu^k would have overflowed
  l <- distrib_d2expected_hessian(d, 0, list(mu = 0, sigma = 1, nu = 1e150),
                                  scale = "link")
  expect_true(all(is.finite(unlist(l))))
})

test_that("the von Mises derivatives keep their digits at a large concentration", {
  # E[l_kappa_kappa] = -A'(kappa), whose derivatives tend to 1/kappa^3 and
  # -3/kappa^4; the recursion for the third derivative of A was 0.61 out at
  # kappa = 1e4 before the Bessel ratio's derivatives took the asymptotic
  # series past kappa = 20
  d <- vonmises1_distrib()
  for (k in c(1e3, 1e4, 1e6)) {
    th <- list(mu = 0, kappa = k)
    a <- distrib_dexpected_hessian(d, 0, th)[["kappa_kappa_kappa"]]
    b <- distrib_d2expected_hessian(d, 0, th)[["kappa_kappa_kappa_kappa"]]
    expect_equal(a * k^3, 1, tolerance = 5 / k)
    expect_equal(b * k^4 / -3, 1, tolerance = 5 / k)
  }
})

test_that("the beta-binomial kernel agrees with the sum over the family's own derivatives", {
  # two routes sharing no arithmetic: the compiled power sums, and the
  # support sum of the family's analytic observed derivatives
  set.seed(5)
  for (d in list(betabinom1_distrib(size = 12), betabinom2_distrib(size = 12))) {
    th <- generate_random_theta(d)
    n <- 30
    thv <- lapply(th, function(v) v * exp(stats::rnorm(n, 0, 0.1)))
    y <- rep(0, n)
    for (k in 1:2) {
      a <- if (k == 1L) distrib_dexpected_hessian(d, y, thv) else
        distrib_d2expected_hessian(d, y, thv)
      b <- support_dexpected(d, y, thv, k, d@size)[names(a)]
      for (nm in names(a)) {
        expect_lt(max(abs(a[[nm]] - b[[nm]])) / max(abs(b[[nm]]), 1e-300), 1e-10,
                  label = paste(class(d)[1], nm))
      }
    }
  }
})
