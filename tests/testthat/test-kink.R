test_that("the three kinked families declare a decomposition and the rest do not", {
  declared <- c("laplace_distrib", "laplace2_distrib", "enet_distrib")
  for (f in declared) {
    d <- do.call(f, list())
    kd <- kink_decomposition(d)
    expect_s3_class(kd, "distributions7::kink_spec")
    expect_identical(kd@phi, "abs")
  }

  # every other family inherits the base method.  The list is read from the
  # namespace rather than typed, so a family added later joins the check.
  ctors <- sort(grep("_distrib$", getNamespaceExports("distributions7"), value = TRUE))
  n_null <- 0L
  for (f in setdiff(ctors, declared)) {
    args <- if (f %in% c("betabinom1_distrib", "betabinom2_distrib", "binomial_distrib"))
              list(size = 10) else if (f == "multinomial_distrib") list(2L, 5) else list()
    d <- tryCatch(do.call(f, args), error = function(e) tryCatch(do.call(f, list(2L)),
                                                                error = function(e) NULL))
    if (is.null(d)) next
    expect_null(kink_decomposition(d), info = f)
    n_null <- n_null + 1L
  }
  expect_gt(n_null, 40L)
})

test_that("params_order deduces the order rather than reading a declaration", {
  expect_identical(params_order(laplace_distrib()),  c(mu = 0, sigma = Inf))
  expect_identical(params_order(laplace2_distrib()), c(mu = 0, lambda = Inf))
  expect_identical(params_order(enet_distrib()),     c(mu = 0, lambda = Inf, alpha = Inf))
  expect_identical(params_order(gaussian1_distrib()), c(mu = Inf, sigma = Inf))

  # the order comes from the composition and the set of parameters v moves
  # with, so a decomposition whose v touches the scale reports the scale too.
  Hub <- S7::new_class("Hub", parent = LaplaceDistrib)
  S7::method(kink_decomposition, Hub) <- function(distrib) {
    kink_spec(phi  = "hinge2",
              v    = function(y, theta) abs(y - theta$mu) - 1.5 * theta$sigma,
              dv   = function(y, theta) list(mu = -sign(y - theta$mu), sigma = -1.5),
              coef = function(theta) -0.5)
  }
  h <- Hub(distrib_name = "hub", dimension = "univariate", bounds = c(-Inf, Inf),
           params = c("mu", "sigma"), n_params = 2L,
           params_bounds = list(mu = c(-Inf, Inf), sigma = c(0, Inf)),
           params_interpretation = c(mu = "location", sigma = "scale"),
           link_params = laplace_distrib()@link_params,
           params_smooth = c(mu = FALSE, sigma = FALSE))
  expect_identical(params_order(h), c(mu = 1, sigma = 1))
})

test_that("params_order agrees with the params_smooth every family declares", {
  ctors <- sort(grep("_distrib$", getNamespaceExports("distributions7"), value = TRUE))
  checked <- 0L
  for (f in ctors) {
    args <- if (f %in% c("betabinom1_distrib", "betabinom2_distrib", "binomial_distrib"))
              list(size = 10) else if (f == "multinomial_distrib") list(2L, 5) else list()
    d <- tryCatch(do.call(f, args), error = function(e) tryCatch(do.call(f, list(2L)),
                                                                error = function(e) NULL))
    if (is.null(d)) next
    ord <- params_order(d)
    sm  <- param_smoothness(d)
    # Inf exactly where the family says smooth, and finite (or NA) where not.
    expect_identical(unname(is.infinite(ord)), unname(sm[names(ord)]), info = f)
    checked <- checked + 1L
  }
  expect_gt(checked, 40L)
})

test_that("a wrapper reports NA rather than Inf, which is the honest state", {
  # No wrapper propagates the declaration yet, so the order is not established.
  # NA says so; Inf would say the parameter is smooth, which it is not.
  w <- truncated(laplace_distrib(), lower = -3)
  expect_null(kink_decomposition(w))
  expect_true(is.na(params_order(w)[["mu"]]))
  expect_identical(params_order(w)[["sigma"]], Inf)

  # fixed() REMOVES the location, and with it the only non-smooth parameter,
  # so what is left is smooth and says so.
  fx <- fixed(laplace_distrib(), mu = 0)
  expect_false("mu" %in% fx@params)
  expect_true(all(is.infinite(params_order(fx))))
})

test_that("check_kink holds a declaration to the family that declares it", {
  for (f in c("laplace_distrib", "laplace2_distrib", "enet_distrib")) {
    r <- check_kink(do.call(f, list()), verbose = FALSE)
    expect_true(r[["dv"]],     info = f)
    expect_true(r[["jump"]],   info = f)
    expect_true(r[["smooth"]], info = f)
  }
  expect_error(check_kink(gaussian1_distrib()), "declares no kink_decomposition")
})

test_that("check_kink REFUSES a wrong declaration", {
  # The injection goes on a SUBCLASS: registering a method mutates the generic
  # in place, so overwriting the real family would break it for the session.
  Wrong <- S7::new_class("Wrong", parent = LaplaceDistrib)
  mk <- function() Wrong(
    distrib_name = "wrong", dimension = "univariate", bounds = c(-Inf, Inf),
    params = c("mu", "sigma"), n_params = 2L,
    params_bounds = list(mu = c(-Inf, Inf), sigma = c(0, Inf)),
    params_interpretation = c(mu = "location", sigma = "scale"),
    link_params = laplace_distrib()@link_params,
    params_smooth = c(mu = FALSE, sigma = TRUE))

  # (a) a v that puts the kink in the wrong place
  S7::method(kink_decomposition, Wrong) <- function(distrib) {
    kink_spec(phi = "abs", v = function(y, theta) y - 2 * theta$mu,
              dv = function(y, theta) list(mu = -2, sigma = 0),
              coef = function(theta) -1 / theta$sigma)
  }
  r <- check_kink(mk(), verbose = FALSE)
  expect_false(r[["jump"]])

  # (b) a coefficient wrong by a factor
  S7::method(kink_decomposition, Wrong) <- function(distrib) {
    kink_spec(phi = "abs", v = function(y, theta) y - theta$mu,
              dv = function(y, theta) list(mu = -1, sigma = 0),
              coef = function(theta) -2 / theta$sigma)
  }
  expect_false(check_kink(mk(), verbose = FALSE)[["jump"]])

  # (c) a dv that does not differentiate v
  S7::method(kink_decomposition, Wrong) <- function(distrib) {
    kink_spec(phi = "abs", v = function(y, theta) y - theta$mu,
              dv = function(y, theta) list(mu = 1, sigma = 0),
              coef = function(theta) -1 / theta$sigma)
  }
  expect_false(check_kink(mk(), verbose = FALSE)[["dv"]])

  # (d) a declaration that names too FEW parameters: the scale is called smooth
  #     while v moves with it.  Without the third check this would pass.
  S7::method(kink_decomposition, Wrong) <- function(distrib) {
    kink_spec(phi = "abs", v = function(y, theta) abs(y - theta$mu) - theta$sigma,
              dv = function(y, theta) list(mu = -sign(y - theta$mu), sigma = 0),
              coef = function(theta) -1 / theta$sigma)
  }
  r <- check_kink(mk(), verbose = FALSE)
  expect_false(isTRUE(r[["dv"]]) && isTRUE(r[["jump"]]) && isTRUE(r[["smooth"]]))
})

test_that("kink_order is the first derivative that jumps", {
  expect_identical(kink_order("abs"),    0)
  expect_identical(kink_order("hinge"),  0)
  expect_identical(kink_order("hinge2"), 1)
  expect_identical(kink_order("hinge3"), 2)
  expect_identical(kink_order("step"),  -1)
  expect_error(kink_order("quadratic"), "must be one of")

  # and the table is what a measurement gives: the jump of the k-th derivative
  # across the origin, kept only where it does not shrink with the step.
  jump <- function(f, k, e) {
    d <- function(x, h) switch(as.character(k),
      "0" = f(x),
      "1" = (f(x + h) - f(x - h)) / (2 * h),
      "2" = (f(x + h) - 2 * f(x) + f(x - h)) / h^2,
      "3" = (f(x + 2*h) - 2*f(x + h) + 2*f(x - h) - f(x - 2*h)) / (2 * h^3))
    abs(d(e, e / 10) - d(-e, e / 10))
  }
  phis <- list(abs = abs, hinge = function(v) pmax(v, 0),
               hinge2 = function(v) pmax(v, 0)^2, hinge3 = function(v) pmax(v, 0)^3,
               step = function(v) as.numeric(v > 0))
  for (nm in names(phis)) {
    f <- phis[[nm]]
    j1 <- vapply(0:3, jump, numeric(1), f = f, e = 1e-4)
    j2 <- vapply(0:3, jump, numeric(1), f = f, e = 1e-5)
    real <- j1 > 1e-6 & j2 > 0.5 * j1          # a real jump does not shrink
    expect_identical(which(real)[1] - 2L, as.integer(kink_order(nm)), info = nm)
  }
})

test_that("the declaration says what a curvature measurement cannot", {
  # The second Bartlett identity fails POINTWISE exactly where the order is
  # finite: the mean of the observed Hessian misses the mass on the kink.
  set.seed(4)
  for (f in c("laplace_distrib", "laplace2_distrib")) {
    d  <- do.call(f, list())
    th <- list(mu = 0.4, sigma = 1.3)
    if (f == "laplace2_distrib") th <- list(mu = 0.4, lambda = 1 / 1.3)
    y  <- distrib_rng(d, 4000L, th)
    H  <- distrib_hessian(d, y, th)
    E  <- distrib_expected_hessian(d, y, th)
    expect_identical(mean(H$mu_mu), 0)             # nothing, pointwise
    expect_lt(mean(E$mu_mu), 0)                    # and the expectation is not
    expect_identical(params_order(d)[["mu"]], 0)
  }
  # the elastic net's kink is real at every alpha and SMALL at a small one,
  # which is why a curvature detector finds it in one sweep and not another.
  d <- enet_distrib()
  for (al in c(0.9, 0.05)) {
    th <- list(mu = 0, lambda = 1, alpha = al)
    expect_identical(params_order(d)[["mu"]], 0)
    expect_equal(kink_decomposition(d)@coef(th), -al, tolerance = 1e-12)
  }
})
