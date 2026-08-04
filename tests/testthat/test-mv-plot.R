# The panel matrix of a multivariate density, and the two arguments that decide
# how the expected information behind a Fisher scoring fit is obtained.

test_that("a marginal is refused for a family that has no closed form", {
  # The base class refuses rather than integrating out the other coordinates
  # numerically, and the plot goes with it, since a panel shows a marginal.
  Odd <- S7::new_class("OddMv", parent = multivariate_distrib, package = NULL)
  d <- Odd(
    distrib_name = "odd", dimension = "multivariate", n_dim = 2L,
    bounds = c(-Inf, Inf), params = c("a", "b"),
    params_interpretation = c(a = "location", b = "location"), n_params = 2L,
    params_bounds = list(a = c(-Inf, Inf), b = c(-Inf, Inf)),
    link_params = list(
      a = linkfunctions7::identity_link(), b = linkfunctions7::identity_link()
    )
  )
  expect_error(mv_marginal(d, list(a = 0, b = 0), 1L), "no closed form")
})

test_that("the panel matrix refuses to be unreadable", {
  d <- mvgaussian_distrib(5)
  th <- generate_random_theta(d)
  expect_error(plot(d, th), "not readable")
  # and says what to do about it
  expect_error(plot(d, th), "which = ")
  # three is the most that is drawn, and choosing three of five is allowed
  expect_silent({
    grDevices::pdf(NULL)
    on.exit(grDevices::dev.off(), add = TRUE)
    plot(d, th, which = c(1, 3, 5))
  })
})

test_that("the panel matrix draws a density and a fit", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)

  d <- mvgaussian_distrib(3)
  th <- as.list(stats::setNames(
    c(0, 1, -1, 0, 0, 0, 0.6, -0.3, 0.2), d@params
  ))
  expect_silent(plot(d, th))
  expect_silent(plot(d, th, which = c(1, 2)))
  # theta may be left out, as for a univariate distribution
  expect_message(plot(d), "random parameters")

  set.seed(51)
  y <- distrib_rng(d, 300, th)
  fit <- fit_distrib(d, y)
  # the fitted picture adds the observations and their kernel density, so the
  # model and the data are read in the same frame
  expect_silent(plot(fit))
  expect_silent(plot(fit, mv_which = c(2, 3)))

  # a heavy-tailed family is drawn from its scale matrix, which exists where
  # the covariance does not
  dt <- mvstudent_t_distrib(2)
  expect_silent(plot(dt, list(mu1 = 0, mu2 = 0, sigma_log_L1 = 0, sigma_log_L2 = 0,
                              sigma_L2.1 = 0.3, nu = 1.5)))
})

test_that("the strategy for the expected information lives on fisher_scoring()", {
  set.seed(52)
  d <- mvstudent_t_distrib(2)
  th <- list(mu1 = 0, mu2 = 0, sigma_log_L1 = 0, sigma_log_L2 = 0, sigma_L2.1 = 0.3, nu = 6)
  y <- distrib_rng(d, 200, th)

  # One argument says how to optimise, and it takes either an optimiser or a
  # fisher_scoring() specification. How the expectation is approximated is a
  # property of Fisher scoring, so it cannot be handed to an optimiser that
  # would never look at it.
  fs <- fisher_scoring(approx = "mc", nsim = 2000)
  expect_s3_class(fs, "distributions7::FisherScoring")
  expect_identical(fs@approx, "mc")
  expect_identical(fs@nsim, 2000)
  expect_false("approx" %in% names(formals(fit_distrib)))
  expect_false("nsim" %in% names(formals(fit_distrib)))

  # The gaussian's expected information is a closed form, so a strategy for
  # approximating it would be ignored, and an argument that is silently
  # ignored is how a caller comes to believe a fit used a method it did not.
  g <- mvgaussian_distrib(2)
  yg <- distrib_rng(g, 200, th[1:5])
  expect_error(
    fit_distrib(g, yg, method = fisher_scoring(approx = "mc")), "closed form"
  )
  expect_error(
    fit_distrib(gaussian_distrib(), rnorm(50),
                method = fisher_scoring(approx = "mc")),
    "closed form"
  )
  # the default carries no claim and is accepted everywhere
  expect_silent(fit_distrib(g, yg, method = fisher_scoring()))

  # and where the strategy does change something it is taken
  set.seed(53)
  expect_silent(fit_distrib(d, y, method = fisher_scoring(approx = "opg")))

  # the object validates its own arguments
  expect_error(fisher_scoring(approx = "nonsense"), "should be one of")
  expect_error(fisher_scoring(nsim = -1), "positive")
  expect_error(fisher_scoring(criterion = "grad"), "stopping rule")

  # and it can carry its own stopping rule and budget
  fs2 <- fisher_scoring(criterion = optimizers7::crit_grad(1e-4), maxit = 3)
  f <- fit_distrib(g, yg, method = fs2)
  expect_lte(f@iterations, 3)
})

test_that("a family without a closed form still fits by Fisher scoring", {
  # Nothing is refused when the argument is left alone: the default strategy
  # applies, whatever the family can do.
  set.seed(54)
  d <- mvstudent_t_distrib(2)
  y <- distrib_rng(d, 400, list(mu1 = 0, mu2 = 0, sigma_log_L1 = 0, sigma_log_L2 = 0,
                                sigma_L2.1 = 0.3, nu = 5))
  fit <- fit_distrib(d, y, method = "fisher")
  expect_true(fit@converged)
})


test_that("a family is asked correctly whether its expected information is exact", {
  # has_exact_expected_hessian() decides whether fisher_scoring()'s 'approx'
  # means anything, so getting it backwards makes the refusal fire on exactly
  # the families the argument exists for. Its argument used to be called
  # 'distrib', which SHADOWS the base class of the same name, so the comparison
  # against that class compared the owning class with the distribution object
  # and every base-class fallback was reported as a closed form.
  exact <- list(gaussian_distrib(), pseudohuber_distrib(), laplace_distrib(),
                weibull_distrib(), gumbel_distrib(), mvgaussian_distrib(2))
  for (d in exact) {
    expect_true(distributions7:::has_exact_expected_hessian(d),
                label = d@distrib_name)
  }

  # These inherit the approximating method -- the skew families from the base
  # class, the multivariate t from the multivariate one -- so the strategy is
  # theirs to choose.
  approximated <- list(skewnormal_distrib(), skewt_distrib(),
                       mvstudent_t_distrib(2))
  for (d in approximated) {
    expect_false(distributions7:::has_exact_expected_hessian(d),
                 label = d@distrib_name)
  }

  # and the consequence, at the level a caller sees it
  set.seed(71)
  ys <- distrib_rng(skewnormal_distrib(), 200,
                    list(mu = 0, sigma = 1, alpha = 3))
  expect_silent(
    fit_distrib(skewnormal_distrib(), ys, method = fisher_scoring(approx = "opg"))
  )
  expect_error(
    fit_distrib(gaussian_distrib(), stats::rnorm(100),
                method = fisher_scoring(approx = "opg")),
    "closed form"
  )
})
