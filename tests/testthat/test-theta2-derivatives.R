# One or two derivatives in the response, two in the parameters.

test_that("the gaussian's closed forms agree with Richardson", {
  skip_if_not_installed("numDeriv")
  y <- c(-2.1, -0.4, 0.3, 1.7, 3.2)
  d <- gaussian1_distrib()
  th <- list(mu = 0.4, sigma = 1.3)
  prs <- list(mu_mu = c("mu", "mu"), sigma_sigma = c("sigma", "sigma"),
              mu_sigma = c("mu", "sigma"))

  for (case in list(list(g = distrib_grad_y_hess(d, y, th),
                         f = function(yy, tt) distrib_grad_y(d, yy, tt)),
                    list(g = distrib_hess_y_hess(d, y, th),
                         f = function(yy, tt) distrib_hess_y(d, yy, tt)))) {
    for (nm in names(prs)) {
      pr <- prs[[nm]]
      ref <- vapply(seq_along(y), function(i) {
        h <- function(v) {
          t2 <- th
          t2[[pr[1L]]] <- v
          numDeriv::grad(function(w) {
            t3 <- t2
            t3[[pr[2L]]] <- w
            case$f(y[i], t3)
          }, t2[[pr[2L]]])
        }
        numDeriv::grad(h, th[[pr[1L]]])
      }, numeric(1))
      expect_equal(case$g[[nm]], ref, tolerance = 1e-5, info = nm)
    }
  }
})

test_that("the gaussian's closed forms are the written-out ones", {
  y <- c(-1.2, 0.5, 2.3)
  mu <- 0.4
  s <- 1.3
  g <- distrib_grad_y_hess(gaussian1_distrib(), y, list(mu = mu, sigma = s))
  h <- distrib_hess_y_hess(gaussian1_distrib(), y, list(mu = mu, sigma = s))
  expect_equal(g$mu_mu, rep(0, 3))
  expect_equal(g$mu_sigma, rep(-2 / s^3, 3), tolerance = 1e-14)
  expect_equal(g$sigma_sigma, -6 * (y - mu) / s^4, tolerance = 1e-14)
  # the curvature carries no location, so only one of its three components
  # survives
  expect_equal(h$mu_mu, rep(0, 3))
  expect_equal(h$mu_sigma, rep(0, 3))
  expect_equal(h$sigma_sigma, rep(-6 / s^4, 3), tolerance = 1e-14)
})

test_that("the fallback agrees with the closed form where both exist", {
  d <- gaussian1_distrib()
  y <- c(-1.3, 0.2, 2.4)
  th <- list(mu = 0.1, sigma = 1.7)
  expect_equal(distrib_grad_y_hess(d, y, th),
               numerical_theta2_y(d, y, th,
                                  function(t) distrib_cross_y(d, y, t)),
               tolerance = 1e-7)
  expect_equal(distrib_hess_y_hess(d, y, th),
               numerical_theta2_y(d, y, th,
                                  function(t) distrib_cross2_y(d, y, t)),
               tolerance = 1e-7)
})

test_that("a family with no closed form still answers", {
  y <- c(0.4, 1.2, 2.9)
  d <- gamma1_distrib()
  th <- list(mu = 1.5, phi = 0.8)
  for (got in list(distrib_grad_y_hess(d, y, th),
                   distrib_hess_y_hess(d, y, th))) {
    expect_named(got, hess_names(d@params))
    expect_true(all(vapply(got, function(v) all(is.finite(v)), logical(1))))
    expect_true(all(lengths(got) == length(y)))
  }
})

test_that("the link scale is the diagonal chain rule at second order", {
  skip_if_not_installed("numDeriv")
  d <- gaussian1_distrib()
  y <- c(-0.7, 0.5, 1.9)
  th <- list(mu = 0.2, sigma = 1.4)
  lnk <- distrib_grad_y_hess(d, y, th, scale = "link")
  # sigma carries a log link: h(e) = exp(e), h' = h'' = sigma, so the diagonal
  # component gains sigma^2 times itself plus sigma times the first order
  par <- distrib_grad_y_hess(d, y, th)
  cr <- distrib_cross_y(d, y, th)
  expect_equal(lnk$sigma_sigma,
               par$sigma_sigma * th$sigma^2 + cr$sigma * th$sigma,
               tolerance = 1e-12)
  expect_equal(lnk$mu_mu, par$mu_mu, tolerance = 1e-12)
  expect_equal(lnk$mu_sigma, par$mu_sigma * th$sigma, tolerance = 1e-12)

  # and against Richardson on the link-scale response gradient, which shares
  # nothing with the chain rule above
  ref <- vapply(seq_along(y), function(i) {
    numDeriv::hessian(function(e) {
      distrib_grad_y(d, y[i], list(mu = e[1], sigma = exp(e[2])))
    }, c(th$mu, log(th$sigma)))[2, 2]
  }, numeric(1))
  expect_equal(lnk$sigma_sigma, ref, tolerance = 1e-5)
})

test_that("fixed() delegates them, subset to the free parameters", {
  d <- gaussian1_distrib()
  f <- fixed(d, mu = 0)
  y <- c(-1.2, 0.3, 2.1)
  g <- distrib_grad_y_hess(f, y, list(sigma = 1.6))
  h <- distrib_hess_y_hess(f, y, list(sigma = 1.6))
  expect_named(g, "sigma_sigma")
  expect_equal(g$sigma_sigma, -6 * y / 1.6^4, tolerance = 1e-14)
  expect_equal(h$sigma_sigma, rep(-6 / 1.6^4, 3), tolerance = 1e-14)
})
