# fixed() on a multivariate family: the link-scale derivatives of the free
# parameters equal the parent's, also when a free parameter carries a link
# other than the identity (where a second application of the link shows).

fixed_mv_link_gap <- function(d, fix) {
  e <- stats::setNames(0.3 + 0.1 * seq_along(d@params), d@params)
  th_full <- stats::setNames(lapply(d@params, function(p)
    linkfunctions7::linkinv(d@link_params[[p]], e[[p]])), d@params)
  for (p in names(fix)) th_full[[p]] <- fix[[p]]
  z <- do.call(fixed, c(list(d), fix))
  free <- z@params
  th <- th_full[free]
  set.seed(11)
  y <- distrib_rng(d, 4, th_full)
  gap <- function(gen, keys) {
    a <- unlist(gen(z, y, th, scale = "link")[keys])
    b <- unlist(gen(d, y, th_full, scale = "link")[keys])
    max(abs(a - b) / pmax(1, abs(b)))
  }
  out <- c(
    gradient = gap(distrib_gradient, free),
    hessian = gap(distrib_hessian, hess_names(free)),
    expected_hessian = gap(distrib_expected_hessian, hess_names(free)),
    deriv3 = gap(distrib_deriv3, deriv_names(free, 3L)),
    deriv4 = gap(distrib_deriv4, deriv_names(free, 4L)))
  if (has_mv_grad_y(d)) {
    out <- c(out,
      cross_y = gap(distrib_cross_y, free),
      cross2_y = gap(distrib_cross2_y, free),
      grad_y_hess = gap(distrib_grad_y_hess, hess_names(free)),
      hess_y_hess = gap(distrib_hess_y_hess, hess_names(free)))
  }
  out
}

test_that("fixed() multivariate applies the link once (nu of the t free)", {
  d <- mvstudent_t1_distrib(2)
  g <- fixed_mv_link_gap(d, list(mu1 = 0, mu2 = 0, sigma_log_L1 = 0.1,
                                 sigma_log_L2 = -0.2, sigma_L2.1 = 0.3))
  expect_lt(max(g), 1e-12)
  g <- fixed_mv_link_gap(d, list(mu1 = 0, mu2 = 0))
  expect_lt(max(g), 1e-12)
})

test_that("fixed() multivariate applies the link once (phi of the Dirichlet free)", {
  g <- fixed_mv_link_gap(dirichlet_distrib(3),
                         list(mean_alr1 = 0.2, mean_alr2 = -0.1))
  expect_lt(max(g), 1e-12)
})
