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
  EnetDistrib = function() enet_distrib()
)

# the constants a family carries besides its parameters, which follow the
# parameters in the vector the registry reads
ccallable_constants <- function(d, n) {
  if (S7::S7_inherits(d, BinomialDistrib)) return(list(size = rep(d@size, length.out = n)))
  list()
}

test_that("every covered class has a constructor in this file", {
  expect_setequal(d7_scalar_classes_covered(), names(ccallable_families))
})

ccallable_twin <- function(cls, eta_range, seed) {
  d <- ccallable_families[[cls]]()
  set.seed(seed)
  n <- 80
  th <- lapply(d@params, function(p)
    linkfunctions7::linkinv(d@link_params[[p]],
                            runif(n, eta_range[1], eta_range[2])))
  names(th) <- d@params
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
}

test_that("the scalar entries are the vector kernels, bit for bit", {
  for (cls in names(ccallable_families)) {
    ccallable_twin(cls, c(-0.5, 0.5), 1)
    # wider, so that the series branches of the remainders are reached
    ccallable_twin(cls, c(-4.5, 5), 2)
  }
})

test_that("an unknown family answers -1", {
  pr <- d7_scalar_probe("NoSuchDistrib", 1L, c(1, 2),
                        cbind(c(1, 1), c(1, 1)))
  expect_identical(pr$id, -1L)
  pr <- d7_info_probe("NoSuchDistrib", 1L, c(1, 2),
                      cbind(c(1, 1), c(1, 1)))
  expect_identical(pr$id, -1L)
})
