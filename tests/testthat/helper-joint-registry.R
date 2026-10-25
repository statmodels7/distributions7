# The families, wrappers and twin comparison of the joint registry's tests,
# shared by test-joint-registry*.R, which split the loops over families so
# that testthat's parallel workers run them at once.

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
  BetaBinom2Distrib = function() betabinom2_distrib(size = 9),
  PoissonDistrib = function() poisson_distrib(),
  BernoulliDistrib = function() bernoulli_distrib(),
  BinomialDistrib = function() binomial_distrib(size = 9),
  ExponentialDistrib = function() exponential_distrib(),
  GeometricDistrib = function() geometric_distrib(),
  ChisqDistrib = function() chisq_distrib()
)

# the classes covered for the gradient and the Hessian only: none since the
# location-scale quadrature of the skew normals, the skew t and the
# pseudo-Huber families is compiled in full (pt_loc_scale.h)
joint_grad_hess_only <- character(0)

joint_twin <- function(cls, seed, eta_range = c(-1.5, 1.5), d = NULL,
                       full = !cls %in% joint_grad_hess_only, far = TRUE) {
  if (is.null(d)) d <- joint_families[[cls]]()
  route <- distrib_scalar_route(d)
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
  far_y <- c(if (is.finite(hi)) hi * (1 - 1e-12) else if (discrete) 1e6 else 1e40,
           if (is.finite(lo)) lo + max(abs(lo), 1) * 1e-12 else -1e40)
  # a count of 1e12 or more is beyond the Poisson-inverse Gaussian's kernels
  # (they allocate one term per count), so the far count is 1e6
  if (discrete) far_y <- round(far_y)
  if (far) y[1:2] <- far_y
  consts <- lapply(route$constants, rep, length.out = n)
  pr <- d7_joint_probe(route$name, y, do.call(cbind, c(th, consts)))
  expect_identical(pr$covered, if (full) 4L else 2L)
  g <- distrib_gradient(d, y, th)
  h <- distrib_hessian(d, y, th)
  # the families whose expected information is the quadrature of
  # loc_scale_expected(), alone or inside a wrapper
  loc_scale <- grepl("skew|pseudohuber", tolower(paste(cls, d@distrib_name)))
  if (full) {
    E <- distrib_expected_hessian(d, y, th)
    dE <- distrib_dexpected_hessian(d, y, th)
    dscale <- max(1e-4, abs(unlist(lapply(dE, function(v) v[is.finite(v)]))))
    escale <- max(1e-4, abs(unlist(lapply(E, function(v) v[is.finite(v)]))))
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
      if (loc_scale) {
        # E_mu,sigma of a symmetric location-scale family is zero
        # analytically: it is read on the scale of the whole information,
        # where its two roundings differ between platforms (clang's fused
        # multiply-adds on macOS: 4e-14 against a floor of 1e-4)
        want <- rep_len(E[[key]], n)
        ok <- is.finite(want)
        expect_lte(max(abs(pr$info[ok, a, b] - want[ok])) / escale, 1e-10,
                   label = paste(cls, "expected", a, b,
                                 "(gap on the scale of the information)"))
      } else {
        expect_agrees_on_scale(pr$info[, a, b],
                               rep_len(E[[key]], n),
                               label = paste(cls, "expected", a, b))
      }
      for (c in seq_len(P)) {
        want <- rep_len(dE[[dexpected_key(d@params, a, b, c)]], n)
        if (loc_scale) {
          # a location-scale family's d_mu E_ab is zero analytically, so it is
          # read on the scale of the whole derivative, not on its own
          ok <- is.finite(want)
          expect_lte(max(abs(pr$dinfo[ok, a, b, c] - want[ok])) / dscale, 1e-10,
                     label = paste(cls, "dexpected", a, b, c,
                                   "(gap on the scale of the derivative)"))
        } else {
          expect_agrees_on_scale(pr$dinfo[, a, b, c], want,
                                 label = paste(cls, "dexpected", a, b, c))
        }
      }
    }
  }
}

# the wrappers, whose joint entries are composed from the inner family's, with
# their coverage: TRUE for all four entries, FALSE for the gradient and the
# Hessian only; and, where given, the two ranges of the predictors. The
# discrete zero-adjusted and truncated families are read at predictors in
# +-0.5 and +-1, as test-ccallable.R reads them: where the retained mass sits
# at the first retained point (f(L)/Z > 0.9) their R generics take
# head_override(), and the registry, scalar and joint, does not. The
# continuous truncations are read there too, so that the interval keeps its
# mass.
narrow <- list(c(-0.5, 0.5), c(-1, 1))
joint_wrappers <- list(
  fixed_t1_nu = list(function() fixed(student_t1_distrib(), nu = 5), TRUE),
  fixed_t1_sigma = list(function() fixed(student_t1_distrib(), sigma = 1.3), TRUE),
  fixed_gengamma1_a = list(function() fixed(gengamma1_distrib(), a = 1.2), TRUE),
  fixed_skewt_nu = list(function() fixed(skewt_distrib(), nu = 6), TRUE),
  zi_poisson = list(function() zero_inflated(poisson_distrib()), TRUE),
  zi_binomial = list(function() zero_inflated(binomial_distrib(size = 9)), TRUE),
  zi_negbin2 = list(function() zero_inflated(negbin2_distrib()), TRUE),
  zi_negbin1 = list(function() zero_inflated(negbin1_distrib()), TRUE),
  zi_betabinom1 = list(function() zero_inflated(betabinom1_distrib(size = 9)), TRUE),
  zi_pig2 = list(function() zero_inflated(pig2_distrib()), TRUE),
  zi_fixed_negbin2 = list(function()
    zero_inflated(fixed(negbin2_distrib(), theta = 2)), TRUE),
  fixed_zi_negbin2 = list(function()
    fixed(zero_inflated(negbin2_distrib()), theta = 2), TRUE),
  fixed_zi_negbin2_zi = list(function()
    fixed(zero_inflated(negbin2_distrib()), zi = 0.3), TRUE),
  za_poisson = list(function() zero_adjusted(poisson_distrib()), TRUE, narrow),
  za_negbin2 = list(function() zero_adjusted(negbin2_distrib()), TRUE, narrow),
  za_binomial = list(function() zero_adjusted(binomial_distrib(size = 7)), TRUE,
                     narrow),
  za_gamma1 = list(function() zero_adjusted(gamma1_distrib()), TRUE),
  za_lognormal1 = list(function() zero_adjusted(lognormal1_distrib()), TRUE),
  za_skewnormal1 = list(function() zero_adjusted(skewnormal1_distrib()), TRUE),
  za_fixed_negbin2 = list(function()
    zero_adjusted(fixed(negbin2_distrib(), theta = 3)), TRUE, narrow),
  za_folded_gaussian1 = list(function()
    zero_adjusted(folded(gaussian1_distrib())), TRUE),
  tf_gamma1_log = list(function()
    transformation(gamma1_distrib(), log_transform()), TRUE),
  tf_gaussian1_exp = list(function()
    transformation(gaussian1_distrib(), exp_transform()), TRUE),
  tf_t1_asinh = list(function()
    transformation(student_t1_distrib(), asinh_transform()), TRUE),
  tf_gaussian1_yj = list(function()
    transformation(gaussian1_distrib(), yj_transform(0.3)), TRUE),
  tf_beta1_logit = list(function()
    transformation(beta1_distrib(), logit_transform()), TRUE),
  tf_skewnormal1_exp = list(function()
    transformation(skewnormal1_distrib(), exp_transform()), TRUE),
  tf_fixed_gaussian1 = list(function()
    transformation(fixed(gaussian1_distrib(), mu = 0), exp_transform()), TRUE),
  fo_gaussian1 = list(function() folded(gaussian1_distrib()), TRUE),
  fo_t1 = list(function() folded(student_t1_distrib()), TRUE),
  fo_laplace = list(function() folded(laplace_distrib()), TRUE),
  fo_cauchy = list(function() folded(cauchy_distrib()), TRUE),
  fo_skewnormal1 = list(function() folded(skewnormal1_distrib()), TRUE),
  fo_fixed_laplace = list(function()
    folded(fixed(laplace_distrib(), sigma = 1.2)), TRUE),
  tr_poisson_lo1 = list(function() truncated(poisson_distrib(), lower = 1),
                        TRUE, narrow),
  tr_poisson_lo4 = list(function() truncated(poisson_distrib(), lower = 4),
                        TRUE, narrow),
  tr_poisson_up6 = list(function() truncated(poisson_distrib(), upper = 6),
                        TRUE),
  tr_negbin2_lo1 = list(function() truncated(negbin2_distrib(), lower = 1),
                        TRUE, narrow),
  tr_negbin2_up8 = list(function() truncated(negbin2_distrib(), upper = 8),
                        TRUE),
  tr_negbin1_both = list(function()
    truncated(negbin1_distrib(), lower = 2, upper = 12), TRUE, narrow),
  tr_betabinom2_both = list(function()
    truncated(betabinom2_distrib(size = 9), lower = 1, upper = 7), TRUE, narrow),
  tr_pig1_lo1 = list(function() truncated(pig1_distrib(), lower = 1), TRUE,
                     narrow),
  tr_fixed_negbin2 = list(function()
    truncated(fixed(negbin2_distrib(), theta = 3), lower = 1), TRUE, narrow),
  zi_tr_poisson = list(function()
    zero_inflated(truncated(poisson_distrib(), upper = 50)), TRUE),
  tr_gaussian1 = list(function()
    truncated(gaussian1_distrib(), lower = -1, upper = 3), TRUE, narrow),
  tr_gamma1_lo = list(function() truncated(gamma1_distrib(), lower = 0.5), TRUE, narrow),
  tr_weibull1_up = list(function() truncated(weibull1_distrib(), upper = 2), TRUE, narrow),
  tr_beta1_up = list(function() truncated(beta1_distrib(), upper = 0.7), TRUE, narrow),
  tr_beta2_lo = list(function() truncated(beta2_distrib(), lower = 0.3), TRUE, narrow),
  tr_t1 = list(function() truncated(student_t1_distrib(), lower = 0), TRUE, narrow),
  tr_fixed_t1 = list(function()
    truncated(fixed(student_t1_distrib(), nu = 5), lower = -0.5, upper = 4), TRUE,
    narrow),
  tr_laplace = list(function()
    truncated(laplace_distrib(), lower = -2, upper = 1), TRUE, narrow)
)
