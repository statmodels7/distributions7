# Drawing a family at several parameter settings on one panel. A component of
# theta given as a vector asks for one curve per element; the univariate
# methods overlay them, and the multivariate one, whose picture is already a
# matrix of panels, rejects the request.

test_that("a vector component becomes one setting per element", {
  d <- gaussian1_distrib()
  ps <- distributions7:::plot_settings(d, list(mu = 0, sigma = c(1, 2, 4)))

  expect_identical(ps$k, 3L)
  expect_identical(ps$varying, "sigma")
  # the held parameter travels into every setting, the varying one splits
  expect_identical(vapply(ps$settings, function(s) s$mu, numeric(1)),
                   c(0, 0, 0))
  expect_identical(vapply(ps$settings, function(s) s$sigma, numeric(1)),
                   c(1, 2, 4))
  # every setting is scalar, which is what distrib_pdf() must be handed: a
  # vector component there means one value per observation instead
  expect_true(all(vapply(ps$settings,
                         function(s) all(lengths(s) == 1L), logical(1))))
})

test_that("one setting is the case k = 1 and nothing varies", {
  d <- gaussian1_distrib()
  ps <- distributions7:::plot_settings(d, list(mu = 0, sigma = 1))
  expect_identical(ps$k, 1L)
  expect_identical(ps$varying, character(0))
  expect_null(distributions7:::plot_labels(d, ps)$legend)

  # a numeric theta is accepted as the methods have always accepted it
  ps2 <- distributions7:::plot_settings(d, c(mu = 0, sigma = 1))
  expect_identical(ps2$settings, ps$settings)
})

test_that("a length that is neither 1 nor k is rejected, not recycled", {
  d <- gaussian1_distrib()
  # 2 divides no part of 3, but the point is stronger than that: R would
  # recycle a length that DOES divide without a word, and a partial setting
  # is far more likely to be a mistake than a request
  expect_error(
    distributions7:::plot_settings(d, list(mu = c(0, 1), sigma = c(1, 2, 3))),
    "length 1 or 3")
  expect_error(
    distributions7:::plot_settings(d, list(mu = numeric(0), sigma = 1)),
    "Empty parameters")
  expect_error(
    distributions7:::plot_settings(d, list(mu = 0)), "Missing parameters")
})

test_that("several settings are told apart without relying on colour", {
  # a continuous family by line type, a discrete one by symbol, so both
  # survive a printed copy that has no colour
  k3 <- distributions7:::plot_keys(3)
  expect_length(k3$col, 3L)
  expect_length(k3$lty, 3L)
  expect_length(k3$pch, 3L)
  expect_identical(anyDuplicated(k3$col), 0L)
  expect_identical(anyDuplicated(k3$lty), 0L)
  expect_identical(anyDuplicated(k3$pch), 0L)

  # one curve keeps the plain appearance the method always had
  k1 <- distributions7:::plot_keys(1)
  expect_identical(k1$col, "black")
  expect_identical(k1$lty, 1L)
  expect_identical(k1$pch, 16L)

  # and the caller's choice wins, recycled over the settings
  own <- distributions7:::plot_keys(3, list(col = "red", lty = 1, pch = 3))
  expect_identical(own$col, rep("red", 3))
  expect_identical(own$lty, rep(1, 3))
  expect_identical(own$pch, rep(3, 3))
})

test_that("the varying parameters go in the legend and the fixed in the title", {
  d <- gaussian1_distrib()
  ps <- distributions7:::plot_settings(d, list(mu = 0, sigma = c(1, 2)))
  labs <- distributions7:::plot_labels(d, ps)
  expect_length(labs$legend, 2L)
  expect_match(labs$legend[1L], "sigma = 1")
  expect_match(labs$main, "mu = 0")
  expect_false(grepl("sigma", labs$main))

  # two parameters moving together are named together, once per curve
  ps2 <- distributions7:::plot_settings(d, list(mu = c(0, 5), sigma = c(1, 2)))
  labs2 <- distributions7:::plot_labels(d, ps2)
  expect_match(labs2$legend[2L], "mu = 5, sigma = 2")
})

test_that("the univariate methods draw every setting without error", {
  f <- tempfile(fileext = ".png")
  grDevices::png(f, width = 400, height = 300)
  on.exit({grDevices::dev.off(); unlink(f)}, add = TRUE)

  expect_silent(plot(gaussian1_distrib(), list(mu = 0, sigma = c(0.7, 1.5, 3))))
  expect_silent(plot(poisson_distrib(), list(mu = c(1, 4, 10))))
  # and one setting still works, which is what every existing caller passes
  expect_silent(plot(gaussian1_distrib(), list(mu = 0, sigma = 1)))
  expect_silent(plot(poisson_distrib(), list(mu = 4)))
})

test_that("the window covers every setting, not the first", {
  # the widest setting has to fit, or the curves the caller asked to compare
  # are drawn cut off at the panel's edge
  d <- gaussian1_distrib()
  f <- tempfile(fileext = ".png")
  grDevices::png(f, width = 400, height = 300)
  plot(d, list(mu = 0, sigma = c(1, 5)))
  got <- graphics::par("usr")[1:2]
  grDevices::dev.off()
  unlink(f)

  # sigma = 5 reaches about -12.9 at its 0.5% quantile
  expect_lt(got[1L], distrib_quantile(d, 0.005, list(mu = 0, sigma = 5)) + 1)
  expect_gt(got[2L], distrib_quantile(d, 0.995, list(mu = 0, sigma = 5)) - 1)
})

test_that("a multivariate family rejects several settings", {
  d <- mvgaussian_distrib(2)
  th <- generate_random_theta(d)
  th[[1L]] <- c(0, 1)
  expect_error(plot(d, th), "one setting")
})
