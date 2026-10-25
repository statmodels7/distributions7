# The families, wrappers and twin comparison of the scalar registry's tests,
# shared by test-ccallable*.R, which split the loop over truncated families
# so that testthat's parallel workers run it at once.

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
    expect_agrees_on_scale(s$score, g[[p]], label = paste(cls, p, "score"))
    expect_agrees_on_scale(s$curvature, h[[paste(p, p, sep = "_")]],
                           label = paste(cls, p, "curvature"))
    expect_agrees_on_scale(e$expected, E[[paste(p, p, sep = "_")]],
                           label = paste(cls, p, "expected"))
    expect_agrees_on_scale(e$dexpected, dE[[paste(p, p, p, sep = "_")]],
                           label = paste(cls, p, "dexpected"))
  }
  if (wrapped) {
    expect_agrees_on_scale(d7_logpdf_probe(cls, y, tm)$logpdf,
                           distrib_pdf(d, y, th, log = TRUE),
                           label = paste(cls, "log-density"))
    expect_identical(d7_scalar_thread_safe_probe(cls), 1L)
  }
}


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
  function() fixed(zero_adjusted(folded(gaussian1_distrib())), sigma = 0.9),
  # truncated() over the discrete families: below only (the removed points,
  # or the retained series when the removed mass is most of the total),
  # above only, and on both sides
  function() truncated(poisson_distrib(), lower = 1),
  function() truncated(poisson_distrib(), lower = 4),
  function() truncated(poisson_distrib(), upper = 6),
  function() truncated(poisson_distrib(), lower = 1, upper = 5),
  function() truncated(negbin2_distrib(), lower = 1),
  function() truncated(negbin2_distrib(), lower = 3),
  function() truncated(negbin2_distrib(), upper = 8),
  function() truncated(negbin1_distrib(), lower = 1),
  function() truncated(negbin1_distrib(), lower = 2, upper = 12),
  function() truncated(geometric_distrib(), lower = 1),
  function() truncated(geometric_distrib(), upper = 9),
  function() truncated(binomial_distrib(size = 7), lower = 1),
  function() truncated(binomial_distrib(size = 7), upper = 5),
  function() truncated(betabinom1_distrib(size = 9), lower = 1),
  function() truncated(betabinom2_distrib(size = 9), lower = 1, upper = 7),
  function() truncated(pig1_distrib(), lower = 1),
  function() truncated(pig2_distrib(), upper = 10),
  function() truncated(fixed(negbin2_distrib(), theta = 3), lower = 1),
  function() zero_inflated(truncated(poisson_distrib(), upper = 50))
)


# the truncation points of a family in the tests of truncated continuous families
ccallable_cut <- function(d) {
  b <- d@bounds
  if (b[1] == 0 && b[2] == 1) c(0.3, 0.8)
  else if (b[1] == 0) c(0.5, 3)
  else if (b[1] < -3 && b[2] < 4) c(-1, 1.5)
  else c(-0.2, 1.2)
}
