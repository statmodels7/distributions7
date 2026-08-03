# Maximum-likelihood estimation on the link scale, checked against closed-form
# MLEs where they exist.

test_that("gaussian fit reproduces the closed-form MLE and its standard errors", {
  set.seed(42)
  d <- gaussian_distrib()
  y <- distrib_rng(d, 2000, list(mu = 2, sigma = 3))
  f <- fit_distrib(d, y, start = list(mu = 0, sigma = 1))

  expect_true(f@converged)
  expect_equal(unname(coef(f)["mu"]), mean(y), tolerance = 1e-6)
  expect_equal(unname(coef(f)["sigma"]), sqrt(mean((y - mean(y))^2)), tolerance = 1e-6)

  # se(mu_hat) = sigma / sqrt(n)
  expect_equal(unname(f@se["mu"]), unname(coef(f)["sigma"]) / sqrt(length(y)), tolerance = 1e-4)
  expect_equal(as.numeric(logLik(f)), sum(distrib_pdf(d, y, as.list(coef(f)), log = TRUE)))
})

test_that("poisson fit reproduces the closed-form MLE and information", {
  set.seed(7)
  d <- poisson_distrib()
  y <- distrib_rng(d, 1000, list(mu = 4.5))
  f <- fit_distrib(d, y, start = list(mu = 1))

  expect_true(f@converged)
  expect_equal(unname(coef(f)), mean(y), tolerance = 1e-7)
  expect_equal(unname(f@se), sqrt(mean(y) / length(y)), tolerance = 1e-5)
})

test_that("bernoulli fit is exact and its interval respects the (0,1) domain", {
  set.seed(11)
  d <- bernoulli_distrib()
  y <- rbinom(60, 1, 0.92)
  f <- fit_distrib(d, y, start = list(mu = 0.5))

  expect_true(f@converged)
  expect_equal(unname(coef(f)), mean(y), tolerance = 1e-7)
  expect_true(all(f@ci > 0 & f@ci < 1))
  expect_lt(f@ci[1, "lower"], coef(f)[[1]])
  expect_gt(f@ci[1, "upper"], coef(f)[[1]])
})

test_that("laplace fit works despite the degenerate observed Hessian in mu", {
  set.seed(13)
  d <- laplace_distrib()
  y <- distrib_rng(d, 2000, list(mu = 1, b = 2))
  f <- fit_distrib(d, y, start = list(mu = 0, b = 1))

  expect_true(f@converged)
  # b_hat is the mean absolute deviation about the fitted location
  expect_equal(unname(coef(f)["b"]), mean(abs(y - coef(f)[["mu"]])), tolerance = 1e-4)
  # the location is at (essentially) the sample median
  expect_equal(unname(coef(f)["mu"]), median(y), tolerance = 0.05)
  expect_true(all(f@se > 0))
})

test_that("all optimisation methods agree", {
  set.seed(3)
  d <- gaussian_distrib()
  y <- distrib_rng(d, 800, list(mu = 1, sigma = 2))
  fits <- lapply(c("fisher", "newton", "bfgs"), function(m) {
    fit_distrib(d, y, start = list(mu = 0, sigma = 1), method = m)
  })
  for (f in fits) expect_true(f@converged)
  expect_equal(coef(fits[[1]]), coef(fits[[2]]), tolerance = 1e-5)
  expect_equal(coef(fits[[1]]), coef(fits[[3]]), tolerance = 1e-4)
})

test_that("extractors and the fitted object are coherent", {
  set.seed(5)
  d <- gamma_distrib()
  y <- distrib_rng(d, 800, list(mu = 3, sigma2 = 2))
  f <- fit_distrib(d, y, start = list(mu = 1, sigma2 = 1))

  expect_s3_class(logLik(f), "logLik")
  expect_equal(attr(logLik(f), "nobs"), length(y))
  expect_equal(dim(vcov(f)), c(2L, 2L))
  expect_equal(dim(vcov(f, scale = "link")), c(2L, 2L))
  expect_named(coef(f), d@params)
  expect_named(coef(f, scale = "link"), d@params)

  # delta method: V_theta = diag(h') V_eta diag(h')
  J <- vapply(seq_along(d@params), function(i) {
    linkfunctions7::linkinvderiv(d@link_params[[d@params[i]]], f@eta[i], order = 1)
  }, numeric(1))
  expect_equal(vcov(f), diag(J) %*% vcov(f, scale = "link") %*% diag(J),
               ignore_attr = TRUE, tolerance = 1e-10)

  # AIC / BIC definitions
  expect_equal(f@aic, -2 * f@loglik + 2 * length(coef(f)))
  expect_equal(f@bic, -2 * f@loglik + log(f@n) * length(coef(f)))
})

test_that("fitting works from the default random starting values", {
  set.seed(9)
  d <- gaussian_distrib()
  y <- distrib_rng(d, 500, list(mu = 2, sigma = 3))
  f <- fit_distrib(d, y)
  expect_true(f@converged)
  expect_equal(unname(coef(f)["mu"]), mean(y), tolerance = 1e-4)
})

test_that("the print method reports both scales", {
  set.seed(2)
  d <- gaussian_distrib()
  y <- distrib_rng(d, 200, list(mu = 0, sigma = 1))
  f <- fit_distrib(d, y, start = list(mu = 0, sigma = 1))
  out <- paste(utils::capture.output(print(f)), collapse = "\n")
  expect_match(out, "Maximum-likelihood fit")
  expect_match(out, "Parameter scale")
  expect_match(out, "Link scale")
  expect_match(out, "Fisher scoring")
})


# --- the optimisation is delegated to optimizers7 ---------------------------

test_that("an optimizers7 optimiser can be passed as the method", {
  set.seed(11)
  d <- gaussian_distrib()
  y <- distrib_rng(d, 400, list(mu = 2, sigma = 3))

  ref <- fit_distrib(d, y, start = list(mu = 1, sigma = 1))

  for (o in list(optimizers7::lbfgs(), optimizers7::cg(),
                 optimizers7::bb(), optimizers7::nelder_mead(maxit = 5000))) {
    f <- fit_distrib(d, y, start = list(mu = 1, sigma = 1), method = o)
    expect_true(f@converged, label = o@name)
    expect_equal(f@loglik, ref@loglik, tolerance = 1e-6, label = o@name)
    # the method reported is the one that actually ran
    expect_identical(f@method, o@name)
  }
})


test_that("the optimiser brings its own stopping rule", {
  set.seed(12)
  d <- gaussian_distrib()
  y <- distrib_rng(d, 300, list(mu = 0, sigma = 2))

  # crit_never() can never fire, so the run must end on the iteration budget
  # and report that it did not converge rather than claiming success.
  f <- fit_distrib(d, y, start = list(mu = 0, sigma = 2),
                   method = optimizers7::lbfgs(
                     criterion = optimizers7::crit_never(), maxit = 3))
  expect_false(f@converged)

  # and a rule the method cannot evaluate is refused by optimizers7, not
  # silently ignored
  expect_error(
    fit_distrib(d, y, start = list(mu = 0, sigma = 2),
                method = optimizers7::nelder_mead(
                  criterion = optimizers7::crit_grad())),
    "does not provide"
  )
})


test_that("an explicitly chosen optimiser is never replaced by the fallback", {
  # fisher and newton fall back to BFGS; a supplied optimiser does not, so a
  # run that fails is reported as a failure under its own name.
  set.seed(13)
  d <- gaussian_distrib()
  y <- distrib_rng(d, 200, list(mu = 1, sigma = 1))

  f <- fit_distrib(d, y, start = list(mu = 1, sigma = 1),
                   method = optimizers7::gd(maxit = 2))
  expect_identical(f@method, optimizers7::gd()@name)
  expect_false(f@converged)
})


test_that("the three named strategies agree with each other and the closed form", {
  set.seed(14)
  d <- gaussian_distrib()
  y <- distrib_rng(d, 500, list(mu = 2, sigma = 3))
  st <- list(mu = 0, sigma = 1)

  fits <- lapply(c("fisher", "newton", "bfgs"),
                 function(m) fit_distrib(d, y, start = st, method = m))
  for (f in fits) expect_true(f@converged)

  lls <- vapply(fits, function(f) f@loglik, numeric(1))
  expect_equal(max(lls) - min(lls), 0, tolerance = 1e-8)

  # the Gaussian MLE is available in closed form
  expect_equal(unname(coef(fits[[1]])["mu"]), mean(y), tolerance = 1e-6)
  expect_equal(unname(coef(fits[[1]])["sigma"]),
               sqrt(mean((y - mean(y))^2)), tolerance = 1e-6)
})


test_that("a non-smooth distribution still fits, and quietly", {
  # The Laplace location sits at a kink. Fisher scoring uses the expected
  # information and reaches the median; the gradient check optimizers7 runs on
  # every call must not mistake the kink for an inconsistent gradient.
  set.seed(15)
  d <- laplace_distrib()
  y <- distrib_rng(d, 400, list(mu = 1, b = 2))

  expect_silent(f <- fit_distrib(d, y, start = list(mu = 0, b = 1)))
  expect_true(f@converged)
  expect_equal(unname(coef(f)["mu"]), stats::median(y), tolerance = 1e-3)
  expect_equal(unname(coef(f)["b"]),
               mean(abs(y - unname(coef(f)["mu"]))), tolerance = 1e-6)
})


# --- confidence intervals ---------------------------------------------------

test_that("confint() agrees with the intervals the fit computed", {
  set.seed(31)
  d <- gamma_distrib()
  y <- distrib_rng(d, 400, list(mu = 4, sigma2 = 6))
  f <- fit_distrib(d, y, start = list(mu = 4, sigma2 = 6))

  expect_equal(unname(confint(f)), unname(f@ci))
  expect_equal(unname(confint(f, scale = "link")), unname(f@ci_eta))

  # the parameter-scale interval is the image of the link-scale one
  ce <- confint(f, scale = "link")
  mapped <- t(vapply(seq_along(d@params), function(i) {
    lk <- d@link_params[[d@params[i]]]
    sort(linkfunctions7::linkinv(lk, ce[i, ]))
  }, numeric(2)))
  expect_equal(unname(mapped), unname(confint(f)))
})


test_that("confint() honours level and parm", {
  set.seed(32)
  d <- gaussian_distrib()
  y <- distrib_rng(d, 300, list(mu = 1, sigma = 2))
  f <- fit_distrib(d, y, start = list(mu = 1, sigma = 2))

  wide <- confint(f, level = 0.99)
  narrow <- confint(f, level = 0.95)
  expect_true(all(wide[, 1] < narrow[, 1]))
  expect_true(all(wide[, 2] > narrow[, 2]))

  expect_equal(rownames(confint(f, "sigma")), "sigma")
  expect_equal(rownames(confint(f, 1)), "mu")
  expect_error(confint(f, "nope"), "does not name a parameter")
})


test_that("an interval never leaves the parameter's domain, either direction", {
  # A bounded-above link is DECREASING, so mapping the two ends swaps them; the
  # limits must still come back in order and inside the domain.
  set.seed(33)
  du <- gaussian_distrib(link_mu = linkfunctions7::bounded_link(upr = 10))
  yu <- distrib_rng(du, 300, list(mu = 6, sigma = 2))
  fu <- fit_distrib(du, yu, start = list(mu = 6, sigma = 2))

  ci <- confint(fu)
  expect_true(all(ci[, 1] < ci[, 2]))
  expect_true(ci["mu", 2] < 10)
  expect_true(ci["sigma", 1] > 0)

  # and on the link scale the interval is symmetric about the estimate
  ce <- confint(fu, scale = "link")
  expect_equal(unname(ce[, 2] - fu@eta), unname(fu@eta - ce[, 1]))

  # a probability stays inside (0, 1)
  set.seed(34)
  fb <- fit_distrib(bernoulli_distrib(), rbinom(60, 1, 0.9),
                    start = list(mu = 0.5))
  cb <- confint(fb)
  expect_true(all(cb > 0 & cb < 1))
})


test_that("the print method shows the interval on both scales", {
  set.seed(35)
  d <- gaussian_distrib()
  y <- distrib_rng(d, 200, list(mu = 0, sigma = 1))
  f <- fit_distrib(d, y, start = list(mu = 0, sigma = 1))
  out <- utils::capture.output(print(f))

  # both blocks carry the two limit columns, not just the estimate
  for (block in c("^Parameter scale", "^Link scale")) {
    i <- grep(block, out)
    expect_length(i, 1)
    header <- out[i + 1]
    expect_match(header, "2\\.5%")
    expect_match(header, "97\\.5%")
  }

  # the headings name the scale and the links, and nothing else
  expect_match(out[grep("^Parameter scale", out)], "^Parameter scale:$")
  expect_match(out[grep("^Link scale", out)], "^Link scale \\(identity, log\\):$")

  # a different level renames the columns
  f90 <- fit_distrib(d, y, start = list(mu = 0, sigma = 1), level = 0.90)
  out90 <- utils::capture.output(print(f90))
  expect_match(out90[grep("^Link scale", out90) + 1], "95%")
  expect_match(out90[grep("^Link scale", out90) + 1], "5%")
})
