# Validates every analytic derivative in the package against numerical
# references:
#   - gradient and observed Hessian vs central finite differences of the
#     log-density;
#   - expected Hessian vs E[observed Hessian] computed with expectation()
#     (numerical integration / series summation).

test_that("analytic gradients match finite differences for all distributions", {
  set.seed(11)
  for (nm in names(all_distrib_cases())) {
    case <- all_distrib_cases()[[nm]]
    y <- distrib_rng(case$d, 25, case$theta)

    a_grad <- distrib_gradient(case$d, y, case$theta)
    n_grad <- fd_gradient_ref(case$d, y, case$theta)

    expect_named(a_grad, names(case$theta), label = paste0(nm, " gradient names"))
    for (p in names(case$theta)) {
      expect_equal(a_grad[[p]], n_grad[[p]],
        tolerance = 1e-5,
        label = sprintf("%s: analytic gradient '%s'", nm, p),
        expected.label = "finite differences"
      )
    }
  }
})

test_that("analytic observed Hessians match finite differences for all distributions", {
  set.seed(22)
  for (nm in names(all_distrib_cases())) {
    case <- all_distrib_cases()[[nm]]
    y <- distrib_rng(case$d, 25, case$theta)

    a_hess <- distrib_hessian(case$d, y, case$theta)
    n_hess <- fd_hessian_ref(case$d, y, case$theta)

    expect_setequal(names(a_hess), names(n_hess))
    for (p in names(n_hess)) {
      expect_equal(a_hess[[p]], n_hess[[p]],
        tolerance = 1e-4,
        label = sprintf("%s: analytic hessian '%s'", nm, p),
        expected.label = "finite differences"
      )
    }
  }
})

test_that("expected Hessians equal the expectation of the observed Hessians", {
  for (nm in names(all_distrib_cases())) {
    case <- all_distrib_cases()[[nm]]
    d <- case$d
    theta <- case$theta

    a_exp <- distrib_expected_hessian(d, 0, theta) # y is ignored, scalar is fine
    obs_names <- names(distrib_hessian(d, distrib_rng(d, 1, theta), theta))
    expect_setequal(names(a_exp), obs_names)

    for (p in names(a_exp)) {
      e_num <- expectation(
        d,
        function(y, theta) distrib_hessian(d, y, theta)[[p]],
        theta
      )
      expect_equal(a_exp[[p]][1], e_num,
        tolerance = 1e-4,
        label = sprintf("%s: analytic expected hessian '%s'", nm, p),
        expected.label = "numerical expectation of observed hessian"
      )
    }
  }
})

test_that("derivatives are vectorized over theta consistently", {
  d <- gaussian_distrib()
  y <- c(-1, 0, 3)
  theta_vec <- list(mu = c(0, 1, 2), sigma = c(1, 2, 3))

  g_vec <- distrib_gradient(d, y, theta_vec)
  for (i in 1:3) {
    g_i <- distrib_gradient(d, y[i], list(mu = theta_vec$mu[i], sigma = theta_vec$sigma[i]))
    expect_equal(g_vec$mu[i], g_i$mu[1])
    expect_equal(g_vec$sigma[i], g_i$sigma[1])
  }
})
