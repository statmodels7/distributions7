# Derivatives on the link (real) scale: scale = "link" applies the diagonal
# Faa di Bruno chain rule to the parameter-scale derivatives.

test_that("identity links leave the derivatives unchanged at every order", {
  set.seed(3)
  d <- gaussian_distrib(link_mu = linkfunctions7::identity_link(),
                        link_sigma = linkfunctions7::identity_link())
  th <- list(mu = 1.5, sigma = 2)
  y <- distrib_rng(d, 8, th)

  expect_equal(distrib_gradient(d, y, th, scale = "link"), distrib_gradient(d, y, th))
  expect_equal(distrib_hessian(d, y, th, scale = "link"), distrib_hessian(d, y, th))
  expect_equal(distrib_deriv3(d, y, th, scale = "link"), distrib_deriv3(d, y, th))
  expect_equal(distrib_deriv4(d, y, th, scale = "link"), distrib_deriv4(d, y, th))
  expect_equal(distrib_expected_hessian(d, y, th, scale = "link"),
               distrib_expected_hessian(d, y, th))
})

test_that("log link reproduces the closed-form chain rule exactly", {
  set.seed(4)
  d <- gaussian_distrib(link_mu = linkfunctions7::identity_link(),
                        link_sigma = linkfunctions7::log_link())
  th <- list(mu = 1.5, sigma = 2)
  y <- distrib_rng(d, 8, th)
  t <- th$sigma

  # theta = exp(eta) => h' = h'' = h''' = h'''' = theta, so the Bell polynomials give
  l1 <- distrib_gradient(d, y, th)$sigma
  l2 <- distrib_hessian(d, y, th)$sigma_sigma
  l3 <- distrib_deriv3(d, y, th)$sigma_sigma_sigma
  l4 <- distrib_deriv4(d, y, th)$sigma_sigma_sigma_sigma

  expect_equal(distrib_gradient(d, y, th, scale = "link")$sigma, l1 * t)
  expect_equal(distrib_hessian(d, y, th, scale = "link")$sigma_sigma,
               l1 * t + l2 * t^2)
  expect_equal(distrib_deriv3(d, y, th, scale = "link")$sigma_sigma_sigma,
               l1 * t + 3 * l2 * t^2 + l3 * t^3)
  expect_equal(distrib_deriv4(d, y, th, scale = "link")$sigma_sigma_sigma_sigma,
               l1 * t + 7 * l2 * t^2 + 6 * l3 * t^3 + l4 * t^4)
})

test_that("link-scale derivatives agree with finite differences in eta", {
  set.seed(5)
  eta_of <- function(d, th) {
    vapply(seq_along(d@params),
           function(i) linkfunctions7::linkfun(d@link_params[[d@params[i]]], th[[i]]), numeric(1))
  }
  theta_of <- function(d, eta) {
    out <- lapply(seq_along(d@params),
                  function(i) linkfunctions7::linkinv(d@link_params[[d@params[i]]], eta[i]))
    names(out) <- d@params
    out
  }
  ev <- function(i, n) { v <- numeric(n); v[i] <- 1; v }

  for (case in list(
    list(d = gaussian_distrib(), th = list(mu = 1.5, sigma = 2)),
    list(d = bernoulli_distrib(), th = list(mu = 0.35)),
    list(d = beta_distrib(), th = list(mu = 0.4, phi = 6))
  )) {
    d <- case$d; th <- case$th
    y <- distrib_rng(d, 5, th)
    eta <- eta_of(d, th)
    p <- length(d@params)
    L <- function(e) distrib_pdf(d, y, theta_of(d, e), log = TRUE)
    fdm <- function(e, idx, h) {
      if (!length(idx)) return(L(e))
      i <- idx[1]
      (fdm(e + h * ev(i, p), idx[-1], h) - fdm(e - h * ev(i, p), idx[-1], h)) / (2 * h)
    }

    g <- distrib_gradient(d, y, th, scale = "link")
    for (nm in names(g)) {
      i <- match(nm, d@params)
      expect_equal(g[[nm]], fdm(eta, i, 1e-4), tolerance = 1e-5,
                   label = paste(d@distrib_name, "link grad", nm))
    }

    h2 <- distrib_hessian(d, y, th, scale = "link")
    for (nm in names(h2)) {
      idx <- match(strsplit(nm, "_", fixed = TRUE)[[1]], d@params)
      expect_equal(h2[[nm]], fdm(eta, idx, 5e-3), tolerance = 1e-3,
                   label = paste(d@distrib_name, "link hess", nm))
    }
  }
})

test_that("expected information transforms as the congruence diag(h') I diag(h')", {
  for (case in list(
    list(d = gaussian_distrib(), th = list(mu = 1.5, sigma = 2)),
    list(d = gamma_distrib(), th = list(mu = 3, sigma2 = 2)),
    list(d = negbin_distrib(), th = list(mu = 4, theta = 1.7))
  )) {
    d <- case$d; th <- case$th
    y <- rep(1, 3)
    eta <- vapply(seq_along(d@params),
                  function(i) linkfunctions7::linkfun(d@link_params[[d@params[i]]], th[[i]]), numeric(1))
    hp <- vapply(seq_along(d@params),
                 function(i) linkfunctions7::linkinvderiv(d@link_params[[d@params[i]]], eta[i], order = 1),
                 numeric(1))
    eh_par <- distrib_expected_hessian(d, y, th)
    eh_lnk <- distrib_expected_hessian(d, y, th, scale = "link")
    for (i in seq_along(d@params)) {
      for (j in i:length(d@params)) {
        nm <- paste(d@params[c(i, j)], collapse = "_")
        expect_equal(eh_lnk[[nm]], eh_par[[nm]] * hp[i] * hp[j],
                     label = paste(d@distrib_name, nm))
      }
    }
  }
})

test_that("scale defaults to the parameter scale and rejects unknown values", {
  d <- gaussian_distrib()
  th <- list(mu = 1, sigma = 2)
  y <- c(0.5, 1.5)
  expect_equal(distrib_gradient(d, y, th), distrib_gradient(d, y, th, scale = "parameter"))
  expect_error(distrib_gradient(d, y, th, scale = "nonsense"))
})
