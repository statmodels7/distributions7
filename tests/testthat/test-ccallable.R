# The scalar C entry points of the fast route: the score and the (k, k)
# second derivative of the log-density in one parameter, and the (k, k)
# expected second derivative with its derivative in the same parameter.
# Each family's header carries one function per component, called by the
# vector kernels and by the registry alike. The two routes agree to the bit
# on gcc x86_64, but on arm64 macOS clang contracts multiply-adds into FMAs
# differently in each context the component is inlined into (and a wrapper's
# R method is a second route outright), so the values are compared on the
# scale of the expected vector, at 1e-10 (expect_agrees_on_scale(), in
# helper-distributions.R). Names, ids and flags stay identical().
# The families, the wrappers and ccallable_twin() are in helper-ccallable.R.

test_that("the route names every covered class and answers NULL elsewhere", {
  for (cls in names(ccallable_families)) {
    r <- distrib_scalar_route(ccallable_families[[cls]]())
    expect_identical(r$name, cls)
    expect_true(is.list(r$constants))
    expect_identical(d7_scalar_thread_safe_probe(cls), 1L)
  }
  expect_identical(distrib_scalar_route(binomial_distrib(size = 7))$constants,
                   list(size = 7))
  expect_null(distrib_scalar_route(mvgaussian1_distrib(n_dim = 2)))
  expect_identical(d7_scalar_thread_safe_probe("NoSuchDistrib"), -1L)
})

test_that("every covered class has a constructor in this file", {
  expect_setequal(d7_scalar_classes_covered(), names(ccallable_families))
})


test_that("the twin comparison passes the platform's last bits and fails a wrong value", {
  # the worst gap the arm64 runner produced (2026-10-08), read against the
  # entry itself as the scale, which is the strictest reading of it
  x <- 0.0076170969780378
  expect_lte(scale_gap(x + 1.14e-13, x), 1e-10)
  expect_lte(scale_gap(c(-0.152444567735, 5.5e-16), c(-0.152444567735, 0)), 1e-10)
  # non-finite entries must match, and a curvature off by one per cent fails
  expect_identical(scale_gap(c(NaN, -Inf, 1), c(NaN, -Inf, 1)), 0)
  expect_identical(scale_gap(c(-Inf, 1), c(-1e300, 1)), Inf)
  expect_gt(scale_gap(c(-0.1524 * 1.01, 0.3), c(-0.1524, 0.3)), 1e-3)
  # a component that is zero analytically: rounding noise on both routes
  # passes (the arm64 runner's worst, 4.2e-15 against 3.7e-15), and a value
  # of 1e-12 against zeros fails
  expect_lte(scale_gap(c(4.163e-15, -1.67e-16, 4.72e-16),
                       c(3.719e-15, -1.39e-16, -1.11e-16)), 1e-10)
  expect_gt(scale_gap(c(1e-12, 0), c(0, 0)), 1e-10)
})

test_that("the scalar entries are the vector kernels", {
  # the Poisson-inverse Gaussian's support sums grow as sigma mu, which the
  # common wide range takes past 1e4; its ranges below reach the series
  # branch near the Poisson limit and the heavy tail at a moderate cost
  pig <- c("Pig1Distrib", "Pig2Distrib")
  for (cls in setdiff(names(ccallable_families), pig)) {
    ccallable_twin(cls, c(-0.5, 0.5), 1)
    # wider, so that the series branches of the remainders are reached
    ccallable_twin(cls, c(-4.5, 5), 2)
  }
  for (cls in pig) ccallable_twin(cls, c(-0.5, 0.5), 1)
  ccallable_twin("Pig1Distrib", NULL, 2,
                 ranges = list(mu = c(-4.5, 1.5), sigma = c(-4.5, 1.5)))
  ccallable_twin("Pig2Distrib", NULL, 2,
                 ranges = list(mu = c(-4.5, 1.5), alpha = c(-1.5, 5)))
  for (cls in loc_scale_compiled()) ccallable_twin(cls, c(-1, 1), 3, TRUE)
})

test_that("the centered skew normal's series region is the R one", {
  expect_identical(sn2_ge_cpp(), sn2_ge())
  d <- skewnormal2_distrib()
  g <- c(-2.9e-3, -1e-3, -2e-5, 3e-6, 4e-4, 2.99e-3, 3.01e-3, -0.2, 0.6)
  th <- list(mu = rep(0.3, 9), sigma = seq(0.5, 2, length.out = 9), gamma1 = g)
  y <- seq(-1, 2, length.out = 9)
  tm <- do.call(cbind, th)
  E <- distrib_expected_hessian(d, y, th)
  dE <- distrib_dexpected_hessian(d, y, th)
  for (k in 1:3) {
    p <- d@params[k]
    e <- d7_info_probe("SkewNormal2Distrib", k, y, tm)
    expect_agrees_on_scale(e$expected, E[[paste(p, p, sep = "_")]])
    expect_agrees_on_scale(e$dexpected, dE[[paste(p, p, p, sep = "_")]])
  }
})

test_that("the compiled quadrature rule is the R rule", {
  expect_identical(loc_scale_rule_cpp(), loc_scale_rule())
})

test_that("an unknown family answers -1", {
  pr <- d7_scalar_probe("NoSuchDistrib", 1L, c(1, 2),
                        cbind(c(1, 1), c(1, 1)))
  expect_identical(pr$id, -1L)
  pr <- d7_info_probe("NoSuchDistrib", 1L, c(1, 2),
                      cbind(c(1, 1), c(1, 1)))
  expect_identical(pr$id, -1L)
})

test_that("the log-density entry is distrib_pdf()'s", {
  for (cls in names(ccallable_families)) {
    d <- ccallable_families[[cls]]()
    set.seed(4)
    n <- 60
    th <- lapply(d@params, function(p)
      linkfunctions7::linkinv(d@link_params[[p]], runif(n, -1, 1)))
    names(th) <- d@params
    y <- vapply(seq_len(n), function(i) distrib_rng(d, 1, lapply(th, `[`, i)), 0)
    # off the support, and between counts for a discrete family
    y[1:4] <- c(-3, -0.5, 1.5, 1e6)
    tm <- do.call(cbind, c(th, ccallable_constants(d, n)))
    got <- d7_logpdf_probe(cls, y, tm)$logpdf
    ref <- suppressWarnings(distrib_pdf(d, y, th, log = TRUE))
    expect_agrees_on_scale(got, ref, label = cls)
  }
  # a near-tie of the square root in D = sqrt(nu + z^2), which MinGW's
  # library sqrt misrounded by one ulp in a build at -O0; this one stays
  # identical(), since the last bit is what it checks
  y <- 1.5092648866770235
  mu <- 0.33886842332719114
  expect_identical(
    d7_logpdf_probe("PseudoHuberDistrib", y, cbind(mu, 1.2, 6))$logpdf,
    distrib_pdf(pseudohuber_distrib(), y, list(mu = mu, sigma = 1.2, nu = 6),
                log = TRUE))
})

# wrapped families: fixed() over families with and without constants

test_that("a wrapped family's entries are the wrapper's methods", {
  for (mk in ccallable_wrappers) {
    d <- mk()
    expect_false(is.null(distrib_scalar_route(d)))
    ccallable_twin(NULL, c(-0.5, 0.5), 5, d = d)
  }
})

test_that("a route the registry cannot read is rejected", {
  expect_identical(d7_scalar_thread_safe_probe("FixedContinuousDistrib:0|Gaussian1Distrib"), -1L)
  expect_identical(d7_scalar_thread_safe_probe("FixedContinuousDistrib:4|Gaussian1Distrib"), -1L)
  expect_identical(d7_scalar_thread_safe_probe("NoWrapper|Gaussian1Distrib"), -1L)
  expect_identical(d7_scalar_thread_safe_probe("ZeroInflatedDistrib:1|PoissonDistrib"), -1L)
  expect_null(distrib_scalar_route(folded(vonmises1_distrib())))
  # a transformer built by hand, or a ready-made one with a function
  # replaced, has no code; a Box-Cox at lambda = 0 is the log transformer
  tr <- log_transform()
  tr@trans_abs_jac <- function(y, log = TRUE) if (log) 2 * y else exp(2 * y)
  expect_null(distrib_scalar_route(transformation(gamma1_distrib(), tr)))
  tr <- transformer(name = "cube", trans_fun = function(x) x^3,
                    trans_inv = function(y) sign(y) * abs(y)^(1 / 3),
                    trans_abs_jac = function(y, log = TRUE) -log(3) - 2 / 3 * log(abs(y)),
                    trans_inv_hessian = function(y) 0, grad_log_jac = function(y) 0,
                    hess_log_jac = function(y) 0, bounds_fun = function(b) b^3,
                    valid_support = function(b) TRUE, decreasing = FALSE)
  expect_null(distrib_scalar_route(transformation(gaussian1_distrib(), tr)))
  expect_identical(
    distrib_scalar_route(transformation(gamma1_distrib(), bc_transform(0)))$name,
    "TransformedDistrib:1|Gamma1Distrib")
  expect_identical(
    distrib_scalar_route(transformation(gaussian1_distrib(), affine_transform(1, -2)))$constants,
    list(loc = 1, scale = -2))
  # a chain is entered once: the same name returns the same id
  nm <- distrib_scalar_route(fixed(zero_inflated(poisson_distrib()), mu = 2))$name
  expect_identical(nm, "FixedDiscreteDistrib:1|ZeroInflatedDistrib|PoissonDistrib")
  id1 <- d7_scalar_probe(nm, 1L, 0, cbind(0.3, 2))$id
  expect_identical(d7_scalar_probe(nm, 1L, 0, cbind(0.3, 2))$id, id1)
  expect_identical(d7_scalar_thread_safe_probe(nm), 1L)
  expect_identical(d7_scalar_thread_safe_probe(
    "FixedContinuousDistrib:2|FixedContinuousDistrib:1|Gaussian1Distrib"), -1L)
  expect_identical(d7_scalar_thread_safe_probe(
    "FoldedDistrib|NoWrapper|Gaussian1Distrib"), -1L)
  # folded() over a parent whose center or scale the registry cannot read
  expect_null(distrib_scalar_route(
    folded(transformation(gaussian1_distrib(), affine_transform(1, 2)))))
  expect_identical(d7_scalar_thread_safe_probe("TransformedDistrib:0|Gaussian1Distrib"), -1L)
  expect_identical(d7_scalar_thread_safe_probe("TransformedDistrib:13|Gaussian1Distrib"), -1L)
})

test_that("folded()'s expected information meets 30-digit values", {
  # mpmath's quadrature of -E[l_k^2] at 30 digits, split at |mu|
  d <- folded(gaussian1_distrib())
  e <- distrib_expected_hessian(d, 1, list(mu = 3, sigma = 0.5))
  expect_equal(e$mu_mu, -3.9999999992522573708, tolerance = 1e-14)
  expect_equal(e$sigma_sigma, -7.9999998923250613929, tolerance = 1e-14)
  e <- distrib_expected_hessian(d, 1, list(mu = 0.5, sigma = 1.2))
  expect_equal(e$mu_mu, -0.18228458070607484782, tolerance = 1e-14)
  expect_equal(e$sigma_sigma, -1.0332223168483544876, tolerance = 1e-14)
  # the Laplace parent, whose kink at its center the rule is split at
  l <- folded(laplace_distrib())
  m <- vapply(list(c(0.5, 1.2), c(3, 0.5), c(-1, 2)), function(p)
    distrib_expected_hessian(l, 1, list(mu = p[1], sigma = p[2]))$mu_mu, 0)
  expect_equal(m, c(-0.273693449951955, -3.99995084660318, -0.115529289315002),
               tolerance = 1e-11)
})

test_that("a truncated discrete family meets 50-digit values", {
  # mpmath, y = 5 under a Poisson truncated below 4 at mu = 0.05, where
  # 1 - F(3) has lost nine digits to the difference of the cdf
  d <- truncated(poisson_distrib(), lower = 4)
  th <- list(mu = 0.05)
  expect_equal(distrib_pdf(d, 5, th, log = TRUE), -4.6152036146879817761,
               tolerance = 1e-15)
  expect_equal(distrib_gradient(d, 5, th)$mu, 19.798660942114209393,
               tolerance = 1e-14)
  # 1.25e-14 relative on arm64 macOS, where 1e-14 was the bound
  expect_equal(distrib_hessian(d, 5, th)$mu_mu, -400.02689585307639218,
               tolerance = 1e-13)
  # the information of a negative binomial truncated below 1, by exact sums
  d <- truncated(negbin2_distrib(), lower = 1)
  th <- list(mu = 2, theta = 0.5)
  expect_equal(distrib_expected_hessian(d, 1, th)$theta_theta,
               -0.58159342305398457498, tolerance = 1e-14)
  expect_equal(distrib_dexpected_hessian(d, 1, th)$theta_theta_mu,
               -0.0063572261128758559268, tolerance = 1e-13)
})

test_that("a truncated continuous family's entries are its methods", {
  for (mk in list(
    function() truncated(gaussian1_distrib(), lower = 0),
    function() truncated(gaussian1_distrib(), upper = 0.5),
    function() truncated(gaussian1_distrib(), lower = -1, upper = 2),
    function() truncated(gaussian1_distrib(), lower = 1.5)
  )) {
    d <- mk()
    expect_false(is.null(distrib_scalar_route(d)))
    ccallable_twin(NULL, c(-0.5, 0.5), 6, d = d)
  }
})

test_that("a continuous family truncated over fixed() has entries equal to its methods", {
  for (mk in list(
    function() truncated(fixed(student_t1_distrib(), nu = 5), lower = 0),
    function() truncated(fixed(gamma1_distrib(), phi = 0.5), lower = 0.5),
    function() truncated(fixed(beta2_distrib(), alpha = 2), upper = 0.8),
    function() truncated(fixed(laplace_distrib(), sigma = 1), lower = -1, upper = 1.5),
    # distrib is named: d alone would partially match it
    function() truncated(fixed(distrib = gengamma1_distrib(), d = 2), lower = 0.5),
    function() truncated(fixed(skewt_distrib(), nu = 6), lower = -0.2),
    function() truncated(fixed(fixed(student_t1_distrib(), nu = 5), sigma = 1.5),
                         upper = 1.2)
  )) {
    d <- mk()
    expect_false(is.null(distrib_scalar_route(d)))
    ccallable_twin(NULL, c(-0.5, 0.5), 8, d = d)
  }
})

test_that("truncating over fixed() is truncating the family at the fixed value", {
  # the information of the free parameters equals the corresponding block of
  # the family truncated at the same point with the parameter held there
  th <- list(mu = c(0.3, -0.4), sigma = c(1.2, 0.7))
  a <- truncated(fixed(student_t1_distrib(), nu = 5), lower = 0)
  b <- truncated(student_t1_distrib(), lower = 0)
  ea <- distrib_expected_hessian(a, c(0.5, 1), th)
  eb <- distrib_expected_hessian(b, c(0.5, 1), c(th, list(nu = 5)))
  for (nm in names(ea)) expect_equal(ea[[nm]], eb[[nm]], tolerance = 1e-14)
  expect_equal(distrib_pdf(a, c(0.5, 1), th, log = TRUE),
               distrib_pdf(b, c(0.5, 1), c(th, list(nu = 5)), log = TRUE),
               tolerance = 1e-15)
})

test_that("a beta truncated up to 1 with a second shape below one keeps its mass at 1", {
  # mpmath at 40 digits in v = -log(1 - y), where the (1 - y)^(b - 1)
  # singularity is a tail; the rule's nodes at y = 1 in double precision were
  # dropped and the information was off by 1.5e-2 at b = 0.24, 0.48 at 0.1
  d <- truncated(beta2_distrib(), lower = 0.3)
  e <- distrib_expected_hessian(d, 0.9, list(alpha = 2, beta = 0.24))
  expect_equal(-e$alpha_alpha, 0.04363137083140498541, tolerance = 1e-14)
  expect_equal(-e$beta_beta, 17.93315887487820139007, tolerance = 1e-14)
  expect_equal(-e$alpha_beta, -0.4494034163784670764, tolerance = 1e-14)
  d <- truncated(beta2_distrib(), lower = 0.5)
  e <- distrib_expected_hessian(d, 0.9, list(alpha = 0.7, beta = 0.1))
  expect_equal(-e$beta_beta, 99.72952477407535529800, tolerance = 1e-14)
  d <- truncated(beta1_distrib(), lower = 0.5)
  e <- distrib_expected_hessian(d, 0.9, list(mu = 0.875, phi = 0.8))
  expect_equal(-e$mu_mu, 64.51204880532808960035, tolerance = 1e-14)
  expect_equal(-e$phi_phi, 1.456160382975719921717, tolerance = 1e-14)
  expect_equal(-e$mu_phi, -9.644939888374605743253, tolerance = 1e-14)
  # and the twins agree, over fixed() too
  for (tr in list(truncated(beta2_distrib(), lower = 0.5),
                  truncated(fixed(beta1_distrib(), phi = 0.8), lower = 0.5))) {
    ccallable_twin(NULL, c(-2, -0.5), 9, d = tr)
  }
})

test_that("a continuous family truncated over folded() or transformation() has entries equal to its methods", {
  for (mk in list(
    function() truncated(folded(gaussian1_distrib()), upper = 2),
    function() truncated(folded(laplace_distrib()), lower = 0.5),
    function() truncated(folded(fixed(student_t1_distrib(), nu = 4)), lower = 0.2, upper = 3),
    function() truncated(transformation(gaussian1_distrib(), exp_transform()), lower = 0.5),
    function() truncated(transformation(gamma1_distrib(), log_transform()), upper = 1),
    function() truncated(transformation(gaussian1_distrib(), affine_transform(1, -2)), lower = 0),
    function() truncated(transformation(gaussian1_distrib(), expit_transform()), upper = 0.8),
    function() truncated(transformation(beta2_distrib(), logit_transform()), lower = -0.5),
    function() truncated(transformation(gamma1_distrib(), inverse_transform()), lower = 0.4)
  )) {
    d <- mk()
    expect_false(is.null(distrib_scalar_route(d)), label = d@distrib_name)
    ccallable_twin(NULL, c(-0.5, 0.5), 10, d = d)
  }
})

test_that("a truncated transformation is the parent truncated at the mapped points", {
  # exp of a gaussian truncated below 0.5 is the gaussian truncated below
  # log(0.5); a decreasing affine map swaps the ends
  th <- list(mu = c(0.2, -0.3), sigma = c(0.8, 1.3))
  a <- truncated(transformation(gaussian1_distrib(), exp_transform()), lower = 0.5)
  b <- truncated(gaussian1_distrib(), lower = log(0.5))
  ea <- distrib_expected_hessian(a, c(1, 2), th)
  eb <- distrib_expected_hessian(b, log(c(1, 2)), th)
  for (nm in names(ea)) expect_equal(ea[[nm]], eb[[nm]], tolerance = 1e-12)
  a <- truncated(transformation(gaussian1_distrib(), affine_transform(1, -2)), lower = 0)
  b <- truncated(gaussian1_distrib(), upper = 0.5)
  ea <- distrib_expected_hessian(a, c(-1, 0.5), th)
  eb <- distrib_expected_hessian(b, c(1, 0.25), th)
  for (nm in names(ea)) expect_equal(ea[[nm]], eb[[nm]], tolerance = 1e-12)
})
