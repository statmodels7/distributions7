# A moment answers with one value per parameter setting.
#
# Every other vectorized quantity in the package does, and the moment methods
# mostly did: 147 of the 168 force the shape with `moment_const(theta, k, 0)`,
# which recycles to the length the family's k parameters imply. Eighteen were
# written without it and answered with the length of whichever components
# entered the formula, so a quantity that does not read the location came back
# length 1 when only the location varied. Nothing warned; a caller binding
# moments to the rows of a data frame got one number recycled down the column.
#
# The sweep below is the check that cannot be satisfied by repairing one
# family: it varies EVERY parameter of EVERY univariate family in turn. Q36's
# own survey varied the first parameter alone and so could not see the three
# `mean` methods that read the first parameter and not a later one.

uni_families <- function() {
  ctors <- sort(grep("_distrib$", getNamespaceExports("distributions7"),
                     value = TRUE))
  out <- list()
  for (nm in ctors) {
    d <- tryCatch(
      if (nm %in% c("betabinom1_distrib", "betabinom2_distrib")) {
        get(nm)(size = 12)
      } else {
        get(nm)()
      },
      error = function(e) NULL)
    if (is.null(d)) next
    if (S7::S7_inherits(d, multivariate_distrib)) next
    out[[nm]] <- d
  }
  out
}


test_that("every moment recycles to the length its parameters imply", {
  set.seed(20260826L)
  fams <- uni_families()
  expect_gte(length(fams), 42L)          # 42 univariate families ship

  gens <- list(mean = mean, variance = variance,
               skewness = skewness, kurtosis = kurtosis)
  for (nm in names(fams)) {
    d <- fams[[nm]]
    th <- generate_random_theta(d)
    for (j in seq_along(th)) {
      t3 <- th
      t3[[j]] <- rep(th[[j]], 3L)         # three identical, valid settings
      for (g in names(gens)) {
        expect_length(gens[[g]](d, t3), 3L)
      }
    }
  }
})


test_that("the eighteen that did not recycle now do, at a written-out theta", {
  # Named one at a time so a failure says which family and which parameter,
  # and with the thetas written out so the case does not move with the seed.
  three <- function(v) rep(v, 3L)

  expect_length(mean(enet_distrib(),
                     list(mu = 0.3, lambda = three(2), alpha = 0.7)), 3L)
  expect_length(variance(enet_distrib(),
                         list(mu = three(0.3), lambda = 2, alpha = 0.7)), 3L)
  expect_length(skewness(enet_distrib(),
                         list(mu = three(0.3), lambda = 2, alpha = 0.7)), 3L)

  expect_length(variance(gumbel_distrib(),
                         list(mu = three(1), sigma = 2)), 3L)
  expect_length(variance(laplace_distrib(),
                         list(mu = three(1), sigma = 2)), 3L)
  expect_length(variance(laplace2_distrib(),
                         list(mu = three(1), lambda = 2)), 3L)

  expect_length(mean(pig1_distrib(), list(mu = 2, sigma = three(1))), 3L)
  expect_length(mean(pig2_distrib(), list(mu = 2, alpha = three(1))), 3L)

  ph <- list(mu = three(0), sigma = 1, nu = 2)
  expect_length(variance(pseudohuber_distrib(), ph), 3L)
  expect_length(kurtosis(pseudohuber_distrib(), ph), 3L)

  sn <- list(mu = three(0), sigma = 1, alpha = 2)
  expect_length(variance(skewnormal1_distrib(), sn), 3L)
  expect_length(skewness(skewnormal1_distrib(), sn), 3L)
  expect_length(kurtosis(skewnormal1_distrib(), sn), 3L)
  expect_length(kurtosis(skewnormal2_distrib(),
                         list(mu = three(0), sigma = 1, gamma1 = 0.3)), 3L)

  st <- list(mu = three(0), sigma = 1, alpha = 2, nu = 9)
  expect_length(variance(skewt_distrib(), st), 3L)
  expect_length(skewness(skewt_distrib(), st), 3L)
  expect_length(kurtosis(skewt_distrib(), st), 3L)

  wb <- list(mu = three(2), sigma = 1.5)
  expect_length(skewness(weibull1_distrib(), wb), 3L)
  expect_length(kurtosis(weibull1_distrib(), wb), 3L)
  expect_length(skewness(weibull3_distrib(), list(mean = three(2), sigma = 1.5)), 3L)
  expect_length(kurtosis(weibull3_distrib(), list(mean = three(2), sigma = 1.5)), 3L)
})


test_that("recycling did not move any value", {
  # `moment_const(theta, k, 0)` is a vector of zeros, so the answer at a scalar
  # theta must be what it was. These are the closed forms, transcribed here
  # rather than read from the methods, so agreement is evidence.
  expect_equal(variance(laplace_distrib(), list(mu = 1, sigma = 2)), 8)
  expect_equal(variance(laplace2_distrib(), list(mu = 1, lambda = 2)), 0.5)
  expect_equal(variance(gumbel_distrib(), list(mu = 1, sigma = 2)),
               pi^2 * 4 / 6)
  expect_equal(mean(pig1_distrib(), list(mu = 2, sigma = 1)), 2)
  expect_equal(mean(pig2_distrib(), list(mu = 2, alpha = 1)), 2)
  expect_equal(mean(enet_distrib(), list(mu = 0.3, lambda = 2, alpha = 0.7)),
               0.3)
  expect_equal(skewness(enet_distrib(), list(mu = 0.3, lambda = 2, alpha = 0.7)),
               0)

  # A Weibull at shape 1 is an exponential, whose skewness is 2 and whose
  # excess kurtosis is 6.
  expect_equal(skewness(weibull1_distrib(), list(mu = 2, sigma = 1)), 2,
               tolerance = 1e-8)
  expect_equal(kurtosis(weibull1_distrib(), list(mu = 2, sigma = 1)), 6,
               tolerance = 1e-8)

  # A skew normal at zero shape is a Gaussian.
  sn0 <- list(mu = 0, sigma = 3, alpha = 0)
  expect_equal(variance(skewnormal1_distrib(), sn0), 9)
  expect_equal(skewness(skewnormal1_distrib(), sn0), 0)
  expect_equal(kurtosis(skewnormal1_distrib(), sn0), 0)

  # A skew t at zero shape and large nu is one too.
  st0 <- list(mu = 0, sigma = 1, alpha = 0, nu = 1e6)
  expect_equal(skewness(skewt_distrib(), st0), 0, tolerance = 1e-6)
  expect_equal(kurtosis(skewt_distrib(), st0), 0, tolerance = 1e-4)
})


test_that("a recycled moment is the scalar one repeated", {
  # The property the fix is for: three identical settings give three identical
  # answers, and each is the answer at one setting.
  d <- skewnormal1_distrib()
  one <- variance(d, list(mu = 0, sigma = 2, alpha = 3))
  many <- variance(d, list(mu = c(-1, 0, 1), sigma = 2, alpha = 3))
  expect_equal(many, rep(one, 3L))

  # And a component that DOES enter the value still varies it.
  vary <- variance(d, list(mu = 0, sigma = c(1, 2, 3), alpha = 3))
  expect_equal(vary, one * c(1, 4, 9) / 4)
})
