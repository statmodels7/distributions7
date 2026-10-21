# The joint entries of the C registry: every component of the gradient, the
# Hessian, the expected information and its derivative at one observation.
# The per-family functions are written from the expressions of the vector
# kernels, so they are compared with the R generics on the scale of the
# expected vector, at 1e-10, as test-ccallable.R compares the scalar entries,
# at interior points and in the tails of the response.

# one constructor per class the joint registry covers
joint_families <- list(
  Gaussian1Distrib = function() gaussian1_distrib(),
  Gamma1Distrib = function() gamma1_distrib(),
  NegBin2Distrib = function() negbin2_distrib(),
  Beta1Distrib = function() beta1_distrib(),
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

# the classes whose expected information is approximated, covered for the
# gradient and the Hessian only
joint_grad_hess_only <- c("SkewNormal1Distrib", "SkewTDistrib",
                          "PseudoHuberDistrib", "PseudoHuber2Distrib",
                          "SkewNormal2Distrib")

test_that("the covered classes are the ones listed here", {
  expect_setequal(d7_joint_classes_covered(), names(joint_families))
})

joint_twin <- function(cls, seed, eta_range = c(-1.5, 1.5)) {
  d <- joint_families[[cls]]()
  P <- length(d@params)
  set.seed(seed)
  n <- 60
  th <- lapply(d@params, function(p)
    linkfunctions7::linkinv(d@link_params[[p]],
                            stats::runif(n, eta_range[1], eta_range[2])))
  names(th) <- d@params
  y <- vapply(seq_len(n), function(i) distrib_rng(d, 1, lapply(th, `[`, i)), 0)
  # the extremes of the response that the family's support admits
  lo <- d@bounds[1]
  hi <- d@bounds[2]
  discrete <- S7::S7_inherits(d, discrete_distrib)
  far <- c(if (is.finite(hi)) hi * (1 - 1e-12) else if (discrete) 1e6 else 1e40,
           if (is.finite(lo)) lo + max(abs(lo), 1) * 1e-12 else -1e40)
  # a count of 1e12 or more is beyond the Poisson-inverse Gaussian's kernels
  # (they allocate one term per count), so the far count is 1e6
  if (discrete) far <- round(far)
  y[1:2] <- far
  consts <- lapply(distrib_scalar_route(d)$constants, rep, length.out = n)
  pr <- d7_joint_probe(cls, y, do.call(cbind, c(th, consts)))
  full <- !cls %in% joint_grad_hess_only
  expect_identical(pr$covered, if (full) 4L else 2L)
  g <- distrib_gradient(d, y, th)
  h <- distrib_hessian(d, y, th)
  if (full) {
    E <- distrib_expected_hessian(d, y, th)
    dE <- distrib_dexpected_hessian(d, y, th)
  } else {
    # the orders beyond the coverage are NaN
    expect_true(all(is.nan(pr$info)) && all(is.nan(pr$dinfo)))
  }
  hn <- hess_names(d@params)
  for (a in seq_len(P)) {
    expect_agrees_on_scale(pr$grad[, a], g[[d@params[a]]],
                           label = paste(cls, "gradient", a))
    for (b in seq_len(P)) {
      key <- hn[hn %in% c(paste(d@params[a], d@params[b], sep = "_"),
                          paste(d@params[b], d@params[a], sep = "_"))][1]
      expect_agrees_on_scale(pr$hess[, a, b], h[[key]],
                             label = paste(cls, "Hessian", a, b))
      if (!full) next
      expect_agrees_on_scale(pr$info[, a, b],
                             rep_len(E[[key]], n),
                             label = paste(cls, "expected", a, b))
      for (c in seq_len(P)) {
        expect_agrees_on_scale(pr$dinfo[, a, b, c],
                               rep_len(dE[[dexpected_key(d@params, a, b, c)]], n),
                               label = paste(cls, "dexpected", a, b, c))
      }
    }
  }
}

test_that("every joint entry is the generic's component", {
  for (cls in names(joint_families)) {
    joint_twin(cls, 11)
    joint_twin(cls, 12, c(-4, 4))
  }
})

test_that("a family the joint registry does not cover answers 0 and NaN", {
  pr <- d7_joint_probe("PoissonDistrib", c(1, 2), cbind(c(1, 1)))
  expect_identical(pr$covered, 0L)
  expect_identical(dim(pr$grad), c(2L, 0L))
  pr <- d7_joint_probe("NoSuchDistrib", c(1, 2), cbind(c(1, 1), c(1, 1)))
  expect_identical(pr$covered, -1L)
})
