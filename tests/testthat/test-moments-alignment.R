# A moment reads theta by NAME.
#
# Every generic here takes theta as a named list and `align_theta()` reorders
# it, strips stray names off the values and validates it against the open
# `params_bounds`. No moment generic aligns before it dispatches, so each
# method has to, and 165 of 168 did. Four did not: they read theta[[1]] and
# theta[[2]] straight, or handed theta to a helper that does, so a caller who
# named the components in another order got another parameter's value back
# with no warning.
#
#   mean(enet_distrib(), list(mu = 0.3, lambda = 2, alpha = 0.7))   #> 0.3
#   mean(enet_distrib(), list(alpha = 0.7, lambda = 2, mu = 0.3))   #> 0.7
#
# The sweep below is the check that cannot be satisfied by repairing one
# family, and the named block beside it says which methods carried the defect.

test_that("every moment is the same however the components are ordered", {
  ctors <- sort(grep("_distrib$", getNamespaceExports("distributions7"),
                     value = TRUE))
  gens <- list(mean = mean, variance = variance, std_dev = std_dev,
               skewness = skewness, kurtosis = kurtosis)
  seen <- 0L

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

    set.seed(20260826L)
    th <- tryCatch(generate_random_theta(d), error = function(e) NULL)
    if (is.null(th) || length(th) < 2L) next     # one component cannot shuffle
    rev_th <- th[rev(seq_along(th))]
    seen <- seen + 1L

    for (g in names(gens)) {
      a <- tryCatch(gens[[g]](d, th), error = function(e) "error")
      b <- tryCatch(gens[[g]](d, rev_th), error = function(e) "error")
      if (identical(a, "error") && identical(b, "error")) next
      expect_equal(b, a, label = sprintf("%s(%s) reversed", g, nm))
    }
  }
  expect_gte(seen, 35L)          # families with at least two components
})


test_that("the four that read theta by position now read it by name", {
  d <- enet_distrib()
  ord <- list(mu = 0.3, lambda = 2, alpha = 0.7)
  rev <- list(alpha = 0.7, lambda = 2, mu = 0.3)

  # The mean of the elastic net is its location, whichever way it is named.
  expect_equal(mean(d, ord), 0.3)
  expect_equal(mean(d, rev), 0.3)

  expect_equal(variance(d, rev), variance(d, ord))
  expect_equal(std_dev(d, rev), std_dev(d, ord))   # follows from variance
  expect_equal(skewness(d, rev), skewness(d, ord))

  v <- vonmises2_distrib()
  vo <- list(mu = 0.3, rho = 0.6)
  vr <- list(rho = 0.6, mu = 0.3)
  expect_equal(mean(v, vr), mean(v, vo))
})


test_that("aligning also validates, which is what the other 165 do", {
  # `align_theta()` checks against `params_bounds` treated as OPEN intervals,
  # so a component outside its domain is an error rather than a number. Before
  # the fix these four returned a value for anything.
  d <- enet_distrib()
  expect_error(mean(d, list(mu = 0.3, lambda = -1, alpha = 0.7)))
  expect_error(variance(d, list(mu = 0.3, lambda = 2, alpha = 1.5)))

  v <- vonmises2_distrib()
  expect_error(mean(v, list(mu = 0.3, rho = 2)))

  # and a missing component is named rather than silently recycled
  expect_error(mean(d, list(mu = 0.3, lambda = 2)))
})


test_that("no value moved at a theta in the family's own order", {
  # Aligning an already-ordered list is the identity, so the numbers a caller
  # who wrote them in order was getting are the numbers they get now. These are
  # transcribed rather than read from the methods.
  d <- enet_distrib()
  th <- list(mu = 0.3, lambda = 2, alpha = 0.7)
  expect_equal(mean(d, th), 0.3)          # the location
  expect_equal(skewness(d, th), 0)        # symmetric about it
  expect_equal(std_dev(d, th), sqrt(variance(d, th)))

  # The elastic net at alpha -> 1 is the Laplace in its rate chart.
  lap <- fixed(laplace2_distrib(), mu = 0)
  expect_equal(variance(d, list(mu = 0, lambda = 2, alpha = 1 - 1e-9)),
               variance(lap, list(lambda = 2)), tolerance = 1e-6)
})
