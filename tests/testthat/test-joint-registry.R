# The joint entries of the C registry: every component of the gradient, the
# Hessian, the expected information and its derivative at one observation.
# The per-family functions are written from the expressions of the vector
# kernels, so they are compared with the R generics on the scale of the
# expected vector, at 1e-10, as test-ccallable.R compares the scalar entries,
# at interior points and in the tails of the response.
# The families, the wrappers and joint_twin() are in helper-joint-registry.R.

test_that("the covered classes are the ones listed here", {
  expect_setequal(d7_joint_classes_covered(), names(joint_families))
})

test_that("every joint entry is the generic's component (part 1 of 2)", {
  for (cls in test_part(names(joint_families), 1, 2)) {
    joint_twin(cls, 11)
    # the skew t below nu = 0.2 is beyond its quadrature on both routes, the R
    # one more (the identity d_sigma E_sigma,sigma = -2 E_sigma,sigma/sigma
    # fails by 4e-6 in C and 1e-4 in R at nu = 0.083, by 2.3e-2 and 1.0e-1 at
    # nu = 0.019): it is read to nu = exp(-1.5)
    joint_twin(cls, 12, if (cls == "SkewTDistrib") c(-1.5, 1.5) else c(-4, 4))
  }
})

test_that("a zero-adjusted family keeps its digits when the mass at zero nears one", {
  # the parent's mass away from zero is -expm1(log f0): formed as 1 - f0 it
  # lost its digits, and the information at mu = 1e-8 was 22 per cent off
  for (mk in list(function() zero_adjusted(poisson_distrib()),
                  function() zero_adjusted(negbin2_distrib()))) {
    d <- mk()
    for (mu in c(1e-4, 1e-6, 1e-8)) {
      th <- list(mu = mu, theta = 5, za = 0.3)[d@params]
      route <- distrib_scalar_route(d)
      M <- matrix(c(unlist(th), unlist(route$constants)), nrow = 1)
      pr <- d7_joint_probe(route$name, 1, M)
      E <- distrib_expected_hessian(d, 1, th)
      g <- distrib_gradient(d, 1, th)
      expect_equal(pr$info[1, 1, 1], E$mu_mu, tolerance = 1e-6,
                   label = paste(d@distrib_name, mu, "information"))
      expect_equal(pr$grad[1, 1], g$mu, tolerance = 1e-6,
                   label = paste(d@distrib_name, mu, "gradient"))
      ip <- d7_info_probe(route$name, 1, 1, M)
      expect_equal(ip$expected, E$mu_mu, tolerance = 1e-6,
                   label = paste(d@distrib_name, mu, "scalar information"))
    }
  }
})

test_that("a zero-inflated family keeps its digits when the mass at zero nears one", {
  # the score in zi at zero and its information, against their closed forms
  # with 1 - f0 = -expm1(-mu): formed as 1 - f0 they were eps/mu off
  d <- zero_inflated(poisson_distrib())
  route <- distrib_scalar_route(d)
  for (mu in c(1e-4, 1e-6, 1e-8)) {
    zi <- 0.3
    q <- -expm1(-mu)
    l0 <- zi + (1 - zi) * exp(-mu)
    score <- q / l0
    info <- -(q^2 / l0 + q / (1 - zi))
    M <- matrix(c(mu, zi, unlist(route$constants)), nrow = 1)
    pr <- d7_joint_probe(route$name, 0, M)
    th <- list(mu = mu, zi = zi)
    expect_equal(pr$grad[1, 2], score, tolerance = 1e-13)
    expect_equal(distrib_gradient(d, 0, th)$zi, score, tolerance = 1e-13)
    expect_equal(pr$info[1, 2, 2], info, tolerance = 1e-13)
    expect_equal(distrib_expected_hessian(d, 0, th)$zi_zi, info,
                 tolerance = 1e-13)
    expect_equal(d7_info_probe(route$name, 2, 0, M)$expected, info,
                 tolerance = 1e-13)
  }
})

test_that("a folded family far in the tail has the parent's Hessian", {
  # far in the tail of a family whose tails separate (a Gaussian, a skew
  # normal), with a positive center, the weight of the positive preimage is
  # exactly one, so the folded score and Hessian are the parent's there
  for (mk in list(function() gaussian1_distrib(),
                  function() skewnormal1_distrib())) {
    par <- mk()
    th <- list(mu = 0.5, sigma = 1.2, alpha = 2)[par@params]
    # (at 1e40, y - mu and -y - mu round to the same magnitude and the
    # weight is one half on any route; at 1e10 the subtraction of R2 - R1^2
    # would still leave a rounding of order 1e4 against a Hessian of order 1)
    y <- c(1e10, 1e6, 1e3)
    M <- matrix(unlist(th), nrow = 3, ncol = length(th), byrow = TRUE)
    pr <- d7_joint_probe(distrib_scalar_route(folded(par))$name, y, M)
    pp <- d7_joint_probe(distrib_scalar_route(par)$name, y, M)
    expect_identical(pr$grad, pp$grad)
    expect_identical(pr$hess, pp$hess)
  }
})

test_that("a folded family without a center and a scale has no information", {
  # folded() has no route over a transformed family; the registry's node is
  # named directly
  tr <- distrib_scalar_route(transformation(gaussian1_distrib(),
                                            asinh_transform()))
  th <- do.call(cbind, c(list(c(1, 1), c(0.5, 0.5)),
                         lapply(tr$constants, rep, length.out = 2)))
  pr <- d7_joint_probe(paste0("FoldedDistrib|", tr$name), c(1, 2), th)
  expect_identical(pr$covered, 2L)
  expect_true(all(is.nan(pr$info)) && all(is.nan(pr$dinfo)))
})

test_that("an unknown family answers -1", {
  pr <- d7_joint_probe("NoSuchDistrib", c(1, 2), cbind(c(1, 1), c(1, 1)))
  expect_identical(pr$covered, -1L)
})
