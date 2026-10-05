# The scalar C entry points of the fast route: the score and the (k, k)
# second derivative of the log-density in one parameter, and the (k, k)
# expected second derivative with its derivative in the same parameter.
# Each family's header carries one function per component, called by the
# vector kernels and by the registry alike, so the two routes must agree
# with identical(), not within a tolerance: the consumer's twin test rests
# on these being the same numbers to the bit.

# one constructor per class the registry covers
ccallable_families <- list(
  Gaussian1Distrib = function() gaussian1_distrib(),
  Gamma1Distrib = function() gamma1_distrib(),
  PoissonDistrib = function() poisson_distrib(),
  NegBin2Distrib = function() negbin2_distrib(),
  Beta1Distrib = function() beta1_distrib(),
  BernoulliDistrib = function() bernoulli_distrib(),
  BinomialDistrib = function() binomial_distrib(size = 7),
  ExponentialDistrib = function() exponential_distrib(),
  GeometricDistrib = function() geometric_distrib(),
  ChisqDistrib = function() chisq_distrib(),
  CauchyDistrib = function() cauchy_distrib(),
  LogisticDistrib = function() logistic_distrib(),
  Gaussian2Distrib = function() gaussian2_distrib(),
  Gaussian3Distrib = function() gaussian3_distrib(),
  Lognormal1Distrib = function() lognormal1_distrib(),
  InvGauss1Distrib = function() invgauss1_distrib(),
  InvGauss2Distrib = function() invgauss2_distrib(),
  Gamma2Distrib = function() gamma2_distrib(),
  GPDDistrib = function() gpd_distrib(),
  VonMises1Distrib = function() vonmises1_distrib(),
  VonMises2Distrib = function() vonmises2_distrib(),
  Weibull3Distrib = function() weibull3_distrib(),
  Lognormal2Distrib = function() lognormal2_distrib(),
  StudentT1Distrib = function() student_t1_distrib(),
  StudentT2Distrib = function() student_t2_distrib(),
  GenGamma1Distrib = function() gengamma1_distrib(),
  GenGamma2Distrib = function() gengamma2_distrib(),
  GumbelDistrib = function() gumbel_distrib(),
  LaplaceDistrib = function() laplace_distrib(),
  Laplace2Distrib = function() laplace2_distrib(),
  Weibull1Distrib = function() weibull1_distrib(),
  Beta2Distrib = function() beta2_distrib(),
  EnetDistrib = function() enet_distrib(),
  NegBin1Distrib = function() negbin1_distrib(),
  SkewNormal1Distrib = function() skewnormal1_distrib(),
  SkewTDistrib = function() skewt_distrib(),
  PseudoHuberDistrib = function() pseudohuber_distrib(),
  PseudoHuber2Distrib = function() pseudohuber2_distrib(),
  SkewNormal2Distrib = function() skewnormal2_distrib(),
  Pig1Distrib = function() pig1_distrib(),
  Pig2Distrib = function() pig2_distrib(),
  BetaBinom1Distrib = function() betabinom1_distrib(size = 9),
  BetaBinom2Distrib = function() betabinom2_distrib(size = 9)
)

# the constants a family carries besides its parameters follow the
# parameters in the vector the registry reads, as distrib_scalar_route()
# lists them
ccallable_constants <- function(d, n) {
  lapply(distrib_scalar_route(d)$constants, rep, length.out = n)
}

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
  expect_null(distrib_scalar_route(truncated(poisson_distrib(), lower = 1)))
  expect_identical(d7_scalar_thread_safe_probe("NoSuchDistrib"), -1L)
})

test_that("every covered class has a constructor in this file", {
  expect_setequal(d7_scalar_classes_covered(), names(ccallable_families))
})

ccallable_twin <- function(cls, eta_range, seed, const_shape = FALSE,
                           ranges = NULL, d = NULL) {
  # a wrapped family is passed as `d` and addressed by its route's name, and
  # its log-density is compared as well
  wrapped <- !is.null(d)
  if (wrapped) cls <- distrib_scalar_route(d)$name
  else d <- ccallable_families[[cls]]()
  set.seed(seed)
  n <- 80
  # `ranges`, where given, holds one eta range per parameter
  th <- lapply(d@params, function(p) {
    rg <- if (is.null(ranges)) eta_range else ranges[[p]]
    linkfunctions7::linkinv(d@link_params[[p]], runif(n, rg[1], rg[2]))
  })
  names(th) <- d@params
  # one shape for every observation, which the quadrature families' cache
  # serves after the first call
  if (const_shape) th[-(1:2)] <- lapply(th[-(1:2)], function(v) rep(v[1], n))
  # one draw per observation at its own parameters: the von Mises rng
  # recycles a per-observation kappa against its proposals
  y <- vapply(seq_len(n), function(i) distrib_rng(d, 1, lapply(th, `[`, i)), 0)
  tm <- do.call(cbind, c(th, ccallable_constants(d, n)))
  g <- distrib_gradient(d, y, th)
  h <- distrib_hessian(d, y, th)
  E <- distrib_expected_hessian(d, y, th)
  dE <- distrib_dexpected_hessian(d, y, th)
  for (k in seq_along(d@params)) {
    p <- d@params[k]
    s <- d7_scalar_probe(cls, k, y, tm)
    e <- d7_info_probe(cls, k, y, tm)
    expect_identical(s$id, e$id)
    expect_identical(s$score, g[[p]], label = paste(cls, p, "score"))
    expect_identical(s$curvature, h[[paste(p, p, sep = "_")]],
                     label = paste(cls, p, "curvature"))
    expect_identical(e$expected, E[[paste(p, p, sep = "_")]],
                     label = paste(cls, p, "expected"))
    expect_identical(e$dexpected, dE[[paste(p, p, p, sep = "_")]],
                     label = paste(cls, p, "dexpected"))
  }
  if (wrapped) {
    expect_identical(d7_logpdf_probe(cls, y, tm)$logpdf,
                     distrib_pdf(d, y, th, log = TRUE),
                     label = paste(cls, "log-density"))
    expect_identical(d7_scalar_thread_safe_probe(cls), 1L)
  }
}

test_that("the scalar entries are the vector kernels, bit for bit", {
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
    expect_identical(e$expected, E[[paste(p, p, sep = "_")]])
    expect_identical(e$dexpected, dE[[paste(p, p, p, sep = "_")]])
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

test_that("the log-density entry is distrib_pdf()'s, bit for bit", {
  loose <- character()
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
    if (cls %in% loose) {
      expect_equal(got, ref, tolerance = 1e-14, label = cls)
    } else {
      expect_identical(got, ref, label = cls)
    }
  }
  # a near-tie of the square root in D = sqrt(nu + z^2), which MinGW's
  # library sqrt misrounded by one ulp in a build at -O0
  y <- 1.5092648866770235
  mu <- 0.33886842332719114
  expect_identical(
    d7_logpdf_probe("PseudoHuberDistrib", y, cbind(mu, 1.2, 6))$logpdf,
    distrib_pdf(pseudohuber_distrib(), y, list(mu = mu, sigma = 1.2, nu = 6),
                log = TRUE))
})

# wrapped families: fixed() over families with and without constants
fix_last <- function(d, value) {
  args <- list(d)
  args[[d@params[d@n_params]]] <- value
  do.call(fixed, args)
}
ccallable_wrappers <- list(
  function() fixed(student_t1_distrib(), nu = 5),
  function() fixed(gaussian1_distrib(), mu = 0.3),
  function() fix_last(gamma1_distrib(), 0.4),
  function() fixed(skewt_distrib(), alpha = 1.5, nu = 6),
  function() fix_last(betabinom1_distrib(size = 9), 0.2),
  function() fix_last(negbin2_distrib(), 3),
  function() zero_inflated(poisson_distrib()),
  function() zero_inflated(negbin2_distrib()),
  function() zero_inflated(negbin1_distrib()),
  function() zero_inflated(betabinom1_distrib(size = 9)),
  function() zero_inflated(pig1_distrib()),
  function() zero_adjusted(poisson_distrib()),
  function() zero_adjusted(negbin2_distrib()),
  function() zero_adjusted(binomial_distrib(size = 7)),
  function() zero_adjusted(gamma1_distrib()),
  function() zero_adjusted(lognormal1_distrib()),
  function() folded(gaussian1_distrib()),
  function() folded(student_t1_distrib()),
  function() folded(laplace_distrib()),
  function() folded(cauchy_distrib()),
  function() folded(skewnormal1_distrib()),
  function() folded(skewnormal2_distrib()),
  # transformation() with each of the twelve ready-made transformers
  function() transformation(gamma1_distrib(), log_transform()),
  function() transformation(gaussian1_distrib(), exp_transform()),
  function() transformation(gamma2_distrib(), inverse_transform()),
  function() transformation(gamma1_distrib(), sqrt_transform()),
  function() transformation(weibull1_distrib(), power_transform(3)),
  function() transformation(gamma1_distrib(), power_transform(0.5)),
  function() transformation(student_t1_distrib(), asinh_transform()),
  function() transformation(gamma1_distrib(), bc_transform(0.5)),
  function() transformation(gamma1_distrib(), bc_transform(-0.7)),
  function() transformation(gaussian1_distrib(), yj_transform(0.3)),
  function() transformation(gaussian1_distrib(), yj_transform(0)),
  function() transformation(gaussian1_distrib(), yj_transform(2)),
  function() transformation(logistic_distrib(), affine_transform(1, -2)),
  function() transformation(beta1_distrib(), logit_transform()),
  function() transformation(gaussian1_distrib(), expit_transform()),
  function() transformation(gamma1_distrib(), softplus_transform(2)),
  # wrappers of wrappers, as far as the constructors allow them
  function() fixed(zero_inflated(poisson_distrib()), mu = 2),
  function() fixed(zero_inflated(negbin2_distrib()), theta = 3),
  function() zero_inflated(fixed(negbin2_distrib(), theta = 3)),
  function() zero_adjusted(fixed(negbin2_distrib(), theta = 3)),
  function() zero_adjusted(fixed(gamma1_distrib(), phi = 0.5)),
  function() zero_adjusted(folded(gaussian1_distrib())),
  function() zero_adjusted(transformation(gaussian1_distrib(), exp_transform())),
  function() folded(fixed(gaussian1_distrib(), sigma = 1.2)),
  function() folded(fixed(laplace_distrib(), sigma = 1.2)),
  function() folded(fixed(student_t1_distrib(), nu = 5)),
  function() transformation(fixed(gaussian1_distrib(), mu = 0), exp_transform()),
  function() transformation(folded(gaussian1_distrib()), log_transform()),
  function() transformation(transformation(gaussian1_distrib(), exp_transform()),
                            log_transform()),
  function() fixed(transformation(gaussian1_distrib(), exp_transform()), sigma = 0.8),
  function() fixed(folded(gaussian1_distrib()), mu = 0.5),
  function() fixed(fixed(skewt_distrib(), nu = 6), alpha = 1.5),
  function() fixed(zero_adjusted(folded(gaussian1_distrib())), sigma = 0.9)
)

test_that("a wrapped family's entries are the wrapper's methods, bit for bit", {
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
  expect_null(distrib_scalar_route(zero_inflated(truncated(poisson_distrib(), upper = 50))))
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
