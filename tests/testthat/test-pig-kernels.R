## One compiled kernel per generic, as every other compiled family has.
##
## The two PIG families were the only ones carrying a single fifteen-column
## block, inherited from the jet the first version propagated: a jet computes
## every partial together by construction, and when the jet was replaced by
## explicit closed forms the block's shape stayed. So the value cost four
## orders, and an R helper stood between every method and its kernel doing
## the recycling, the support mask and the naming that every other family's
## kernel does for itself.
##
## The split moved statements and rewrote none, so the arithmetic must be
## unchanged BIT FOR BIT rather than to a tolerance: each kernel is asked to
## reproduce the block's own columns for its order. That is the whole
## verification of the refactor -- a statement placed in the wrong stage
## either fails to compile or fails this.

## the block columns each kernel returns, IN THE ORDER IT RETURNS THEM. The
## hessian is a permutation (mu_mu, p2_p2, mu_p2 against d20, d11, d02), so a
## mapping written by position rather than by name would transpose two entries
## of the information and this is what would catch it.
pig_layout <- list(pdf = 1L, gradient = 2:3, hessian = c(4L, 6L, 5L),
                   deriv3 = 7:10, deriv4 = 11:15)

pig_flat <- function(x) if (is.list(x)) unlist(x, use.names = FALSE) else
  as.vector(x)

test_that("each kernel is the block's own columns for its order", {
  set.seed(11)
  n <- 400
  y <- distrib_rng(pig1_distrib(), n, list(mu = 4, sigma = 0.9))
  mu <- stats::runif(n, 0.05, 60)
  al <- stats::runif(n, 0.05, 40)
  sg <- stats::runif(n, 0.02, 20)

  for (fam in c("pig1", "pig2")) {
    p2 <- if (fam == "pig1") sg else al
    blk <- get(paste0(fam, "_hd_cpp"))(y, mu, p2, 1L)
    for (g in names(pig_layout)) {
      got <- get(paste0(fam, "_", g, "_cpp"))(y, mu, p2, 1L)
      ref <- as.vector(blk[, pig_layout[[g]], drop = FALSE])
      expect_identical(pig_flat(got), ref, info = paste(fam, g))
    }
  }
})

test_that("the kernel tests the support itself", {
  ## Every other compiled family guards the support inside the kernel --
  ## betabinom.cpp writes R_NegInf for a fractional response -- and pig now
  ## does the same, so there is no R-side mask between the method and the
  ## kernel. The mass reads -Inf and a derivative reads NaN, there being no
  ## derivative to report at a point with no mass.
  y <- c(-1, 1.5, 2)
  m <- rep(3, 3)
  s <- rep(0.8, 3)
  v <- pig1_pdf_cpp(y, m, s, 1L)
  expect_identical(v[1:2], c(-Inf, -Inf))
  expect_true(is.finite(v[3]))
  for (g in c("gradient", "hessian", "deriv3", "deriv4")) {
    got <- get(paste0("pig1_", g, "_cpp"))(y, m, s, 1L)
    expect_true(all(vapply(got, function(u) all(is.nan(u[1:2])), logical(1))),
                info = g)
    expect_true(all(vapply(got, function(u) is.finite(u[3]), logical(1))),
                info = g)
  }
  expect_identical(distrib_pdf(pig1_distrib(), c(-1, 1.5),
                               list(mu = 3, sigma = 0.8)), c(0, 0))
})

test_that("a scalar parameter is recycled inside the kernel", {
  ## the R helper used to do this. The kernels read the scalar branch INSIDE
  ## the loop rather than hoisting it into a variable the workers write, which
  ## is the data race d7_par.h warns of and which shows only where the
  ## parameter varies by observation -- so the thread comparison below is run
  ## with a parameter per observation as well.
  y <- as.numeric(0:6)
  a <- pig1_gradient_cpp(y, 3, 0.8, 1L)
  expect_identical(a, pig1_gradient_cpp(y, rep(3, 7), rep(0.8, 7), 1L))
  expect_identical(a, pig1_gradient_cpp(y, 3, 0.8, 2L))
})

test_that("the derivative kernels return a named list, as every family's do", {
  y <- as.numeric(0:6)
  m <- rep(3, 7)
  s <- rep(0.8, 7)
  expect_identical(names(pig1_gradient_cpp(y, m, s, 1L)), c("mu", "sigma"))
  expect_identical(names(pig1_hessian_cpp(y, m, s, 1L)),
                   c("mu_mu", "sigma_sigma", "mu_sigma"))
  expect_identical(names(pig2_hessian_cpp(y, m, s, 1L)),
                   c("mu_mu", "alpha_alpha", "mu_alpha"))
  expect_identical(names(pig2_deriv4_cpp(y, m, s, 1L)),
                   c("mu_mu_mu_mu", "mu_mu_mu_alpha", "mu_mu_alpha_alpha",
                     "mu_alpha_alpha_alpha", "alpha_alpha_alpha_alpha"))
  ## and the method is the kernel, with nothing between them
  expect_identical(pig1_hessian_cpp(y, m, s, 1L),
                   distrib_hessian(pig1_distrib(), y,
                                   list(mu = 3, sigma = 0.8)))
})

test_that("the value does not pay for four orders", {
  ## what the split is FOR. The Bessel sum is shared by every order and is
  ## most of the cost, so the gain is modest and the check is structural: the
  ## mass comes back as one vector, so it cannot be computing the fourth
  ## derivative and discarding it.
  y <- as.numeric(0:6)
  m <- rep(3, 7)
  s <- rep(0.8, 7)
  expect_null(dim(pig1_pdf_cpp(y, m, s, 1L)))
  expect_identical(length(pig1_pdf_cpp(y, m, s, 1L)), 7L)
  expect_identical(length(pig2_pdf_cpp(y, m, s, 1L)), 7L)
  expect_identical(length(pig1_gradient_cpp(y, m, s, 1L)), 2L)
  expect_identical(length(pig1_hessian_cpp(y, m, s, 1L)), 3L)
  expect_identical(length(pig1_deriv3_cpp(y, m, s, 1L)), 4L)
  expect_identical(length(pig1_deriv4_cpp(y, m, s, 1L)), 5L)
})

test_that("the kernels are identical at any thread count", {
  set.seed(12)
  n <- 3000
  y <- distrib_rng(pig2_distrib(), n, list(mu = 5, alpha = 2))
  mu <- stats::runif(n, 0.5, 20)
  al <- stats::runif(n, 0.5, 10)
  for (g in names(pig_layout)) {
    f <- get(paste0("pig2_", g, "_cpp"))
    expect_identical(f(y, mu, al, 1L), f(y, mu, al, 2L), info = g)
  }
})
