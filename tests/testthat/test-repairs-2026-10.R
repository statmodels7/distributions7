# Regression tests of the repairs of distributions7 0.84.0: negbin2's
# derivatives in theta at y = 0 and y = 1, negbin1's series at large theta,
# and the Poisson-inverse Gaussian at large counts. The exact values are
# mpmath references (misure_gas_congiunto/c/nb2_small_ref.py, nb_ref.py).

test_that("negbin2's derivatives in theta at y = 1 hold at a small mean", {
  d <- negbin2_distrib()
  th <- list(mu = 1e-8, theta = 50)
  # orders one to four in theta at y = 1, 60 digits
  exact <- c(3.9999999792000000054933333321e-12,
             -1.59999999152000000332799999901e-13,
             9.59999994816000002687999999002e-15,
             -7.67999995776000002713599998794e-16)
  got <- c(distrib_gradient(d, 1, th)$theta,
           distrib_hessian(d, 1, th)$theta_theta,
           distrib_deriv3(d, 1, th)$theta_theta_theta,
           distrib_deriv4(d, 1, th)$theta_theta_theta_theta)
  expect_lte(max(abs(got / exact - 1)), 1e-12)
})

test_that("negbin1's expected information comes from its series at large theta", {
  d <- negbin1_distrib()
  th <- list(mu = 1, theta = 50)
  # E[l_mu,mu], E[l_mu,theta], E[l_theta,theta] and their derivatives, 25 digits
  E <- distrib_expected_hessian(d, 0, th)
  expect_lte(abs(E$mu_mu / -0.07565578418997869659984071 - 1), 1e-13)
  expect_lte(abs(E$mu_theta / 0.001120958821054475892781128 - 1), 1e-13)
  expect_lte(abs(E$theta_theta / -0.00002241917642108951785562256 - 1), 1e-13)
  dE <- distrib_dexpected_hessian(d, 0, th)
  expect_lte(abs(dE[[dexpected_key(d@params, 1, 1, 1)]] /
                   0.0785286679988877983582472 - 1), 1e-13)
  expect_lte(abs(dE[[dexpected_key(d@params, 2, 2, 2)]] /
                   0.000001180153592665565645372124 - 1), 1e-12)
  # a row at theta = 1e8 cost 11 s in the support sums; the series is finite
  # and immediate
  t <- system.time(E8 <- distrib_expected_hessian(d, 0, list(mu = 1, theta = 1e8)))
  expect_true(all(is.finite(unlist(E8))))
  expect_lt(t[["elapsed"]], 1)
})

test_that("the Poisson-inverse Gaussian's derivatives hold at large counts", {
  d <- pig2_distrib()
  th <- list(mu = 3000, alpha = 0.4)
  # a count of 1e12 threw std::length_error when the sum held one term per count
  expect_true(all(is.finite(unlist(distrib_gradient(d, 1e12, list(mu = 1e12, alpha = 0.4))))))
  expect_true(all(is.finite(unlist(distrib_hessian(d, 1e12, list(mu = 1e12, alpha = 0.4))))))
  # orders one to four in alpha at y = 2000, against mpmath at 60 digits
  # (misure_gas_congiunto/c/fase2/pig_test/pig2_ref.py): the old sum's raw
  # moments left 2.9e-6 at order three and 3.7e-3 at order four
  gens <- list(distrib_gradient, distrib_hessian, distrib_deriv3,
               distrib_deriv4)
  for (cs in list(
    list(th = th, exact = c(2.4999888138330068323, -6.2500279654154321068,
                            31.250000000015381094, -234.37499999996154727)),
    list(th = list(mu = 1500, alpha = 40),
         exact = c(0.023882385252430571133, -0.00065289030714207751769,
                   0.000031253750809678121884, -2.3436567051416649393e-6)))) {
    got <- vapply(1:4, function(o)
      gens[[o]](d, 2000, cs$th)[[paste(rep("alpha", o), collapse = "_")]], 0)
    expect_lte(max(abs(got / cs$exact - 1)), 1e-11)
  }
})
