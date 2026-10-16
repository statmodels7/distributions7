# The adapters on numericals7::fd_tensor() against the R implementations of
# distributions7 0.78.1 (helper-tensor-twins.R), and their accuracy against
# the closed forms.

normgap <- function(a, b) {
  a <- unlist(a[names(b)])
  b <- unlist(b)
  sqrt(sum((a - b)^2)) / max(sqrt(sum(b^2)), 1e-300)
}

twin_cases <- function() {
  c(all_distrib_cases(), list(
    truncated = list(d = truncated(gaussian1_distrib(), lower = -1, upper = 2),
                     theta = list(mu = 0.4, sigma = 1.3)),
    zero_inflated = list(d = zero_inflated(poisson_distrib()),
                         theta = list(mu = 3, zi = 0.2)),
    zero_adjusted = list(d = zero_adjusted(gamma2_distrib()),
                         theta = list(mu = 2, sigma2 = 0.6, za = 0.15)),
    folded = list(d = folded(gaussian1_distrib()),
                  theta = list(mu = 0.5, sigma = 1.1)),
    transformation = list(d = transformation(gamma2_distrib(), log_transform()),
                          theta = list(mu = 2, sigma2 = 0.5)),
    fixed = list(d = fixed(student_t1_distrib(), nu = 5),
                 theta = list(mu = 0.3, sigma = 1.2))))
}

test_that("at accuracy 2 the adapters reproduce the R implementations", {
  # The same nodes and steps: every route is bit-identical except the fourth
  # order from the log-density, whose one-dimensional weights fd_weights()
  # solves to 1.8e-15 and Fornberg's recursion gives exactly; the gap is the
  # rounding of the stencil itself (the error against the closed forms is
  # 2e-5 to 2e-4 on both routes).
  set.seed(3)
  cases <- twin_cases()
  for (nm in names(cases)) {
    d <- cases[[nm]]$d
    th <- cases[[nm]]$theta
    y <- distrib_rng(d, 6, th)
    if (nm == "zero_adjusted") y <- y[y > 0]
    for (order in 3:4) for (base in 0:(order - 1L)) {
      if (base > 3L) next
      g <- normgap(tensor_derivatives(d, y, th, order, base, accuracy = 2L),
                   twin_tensor_derivatives(d, y, th, order, base))
      tol <- if (order == 4L && base == 0L) 1e-4 else 1e-13
      expect_lte(g, tol, label = sprintf("%s order %d base %d", nm, order, base))
    }
    if (has_exact_deriv4(d)) {
      for (sc in c("parameter", "link")) {
        g <- normgap(numerical_deriv5(d, y, th, scale = sc, accuracy = 2L),
                     twin_numerical_deriv5(d, y, th, scale = sc))
        expect_lte(g, 1e-13, label = sprintf("%s deriv5 %s", nm, sc))
      }
    }
    if (S7::S7_inherits(d, continuous_distrib)) {
      for (ac in list(c(1, 1), c(2, 1), c(3, 1), c(1, 2), c(2, 2))) {
        g <- normgap(
          mixed_tensor_derivatives(d, y, th, ac[1], ac[2], accuracy = 2L),
          twin_mixed_tensor_derivatives(d, y, th, ac[1], ac[2]))
        expect_lte(g, 1e-13, label = sprintf("%s mixed %d,%d", nm, ac[1], ac[2]))
      }
    }
  }
})

bare_family <- function(name, pdf, template) {
  cls <- S7::new_class(name, parent = continuous_distrib)
  S7::method(distrib_pdf, cls) <- function(distrib, y, theta, log = FALSE, ...) {
    pdf(y, theta, log)
  }
  cls(distrib_name = name, dimension = "univariate", bounds = template@bounds,
      params = template@params, n_params = template@n_params,
      params_bounds = template@params_bounds,
      link_params = template@link_params,
      params_interpretation = template@params_interpretation)
}

test_that("a density alone meets the closed forms at orders 3 to 5", {
  g1 <- gaussian1_distrib()
  t1 <- student_t1_distrib()
  bare_g <- bare_family("BareGaussTwin", function(y, th, log)
    stats::dnorm(y, th$mu, th$sigma, log = log), g1)
  bare_t <- bare_family("BareTTwin", function(y, th, log) {
    l <- stats::dt((y - th$mu) / th$sigma, th$nu, log = TRUE) - log(th$sigma)
    if (log) l else exp(l)
  }, t1)
  set.seed(1)
  thg <- list(mu = 0.3, sigma = 1.2)
  tht <- list(mu = 0.3, sigma = 1.2, nu = 5)
  yg <- distrib_rng(g1, 500, thg)
  yt <- distrib_rng(t1, 500, tht)
  # measured at n = 2000 with numericals7 0.23.0: gaussian 1.1e-7, 1.4e-6,
  # 5.1e-5; t 3.3e-8, 4.8e-7, 1.8e-5
  expect_lt(normgap(distrib_deriv3(bare_g, yg, thg), distrib_deriv3(g1, yg, thg)), 1e-6)
  expect_lt(normgap(distrib_deriv4(bare_g, yg, thg), distrib_deriv4(g1, yg, thg)), 1e-5)
  expect_lt(normgap(distrib_deriv5(bare_g, yg, thg), distrib_deriv5(g1, yg, thg)), 3e-4)
  expect_lt(normgap(distrib_deriv3(bare_t, yt, tht), distrib_deriv3(t1, yt, tht)), 1e-6)
  expect_lt(normgap(distrib_deriv4(bare_t, yt, tht), distrib_deriv4(t1, yt, tht)), 1e-5)
  expect_lt(normgap(distrib_deriv5(bare_t, yt, tht), distrib_deriv5(t1, yt, tht)), 2e-4)
})

test_that("the fifth order of a density-only family is no longer nested", {
  g1 <- gaussian1_distrib()
  bare <- bare_family("BareGaussNest", function(y, th, log)
    stats::dnorm(y, th$mu, th$sigma, log = log), g1)
  expect_identical(deriv5_base(bare), 0L)
  expect_identical(deriv5_base(g1), 4L)
  # a wrapper of a density-only parent differences its own log-density
  expect_identical(deriv5_base(fixed(bare, mu = 0)), 0L)
  set.seed(2)
  th <- list(mu = 0.3, sigma = 1.2)
  y <- distrib_rng(g1, 500, th)
  ref <- distrib_deriv5(g1, y, th)
  # the route of 0.78.1: a first difference of the fourth-order stencil
  old <- normgap(twin_numerical_deriv5(bare, y, th, deriv4 = function(t)
    twin_tensor_derivatives(bare, y, t, 4L, 0L)), ref)
  new <- normgap(numerical_deriv5(bare, y, th), ref)
  expect_gt(old, 1e-2)
  expect_lt(new, 3e-4)
})

test_that("accuracy 4 improves the mixed derivatives of a density alone", {
  g1 <- gaussian1_distrib()
  bare <- bare_family("BareGaussMixed", function(y, th, log)
    stats::dnorm(y, th$mu, th$sigma, log = log), g1)
  set.seed(4)
  th <- list(mu = 0.3, sigma = 1.2)
  y <- distrib_rng(g1, 200, th)
  refs <- list(distrib_cross_y(g1, y, th), distrib_cross2_y(g1, y, th),
               distrib_grad_y_hess(g1, y, th))
  ac <- list(c(1, 1), c(2, 1), c(1, 2))
  for (k in seq_along(ac)) {
    e2 <- normgap(mixed_tensor_derivatives(bare, y, th, ac[[k]][1], ac[[k]][2],
                                           accuracy = 2L), refs[[k]])
    e4 <- normgap(mixed_tensor_derivatives(bare, y, th, ac[[k]][1], ac[[k]][2]),
                  refs[[k]])
    expect_lt(e4, e2)
    expect_lt(e4, 1e-6)
  }
})
