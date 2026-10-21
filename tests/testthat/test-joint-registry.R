# The joint entries of the C registry: every component of the gradient, the
# Hessian, the expected information and its derivative at one observation.
# They call the same per-component functions as the vector kernels, so they
# are compared with the R generics on the scale of the expected vector, at
# 1e-10, as test-ccallable.R compares the scalar entries.

joint_families <- list(
  Gaussian1Distrib = function() gaussian1_distrib(),
  StudentT1Distrib = function() student_t1_distrib()
)

test_that("the covered classes are the ones listed here", {
  expect_setequal(d7_joint_classes_covered(), names(joint_families))
})

test_that("every joint entry is the generic's component", {
  for (cls in names(joint_families)) {
    d <- joint_families[[cls]]()
    P <- length(d@params)
    set.seed(11)
    n <- 60
    th <- lapply(d@params, function(p)
      linkfunctions7::linkinv(d@link_params[[p]], runif(n, -1.5, 1.5)))
    names(th) <- d@params
    y <- vapply(seq_len(n), function(i) distrib_rng(d, 1, lapply(th, `[`, i)), 0)
    # far in the tails, where the components are formed with care
    y[1:3] <- c(1e40, -1e40, 0)
    pr <- d7_joint_probe(cls, y, do.call(cbind, th))
    expect_identical(pr$covered, 1L)
    g <- distrib_gradient(d, y, th)
    h <- distrib_hessian(d, y, th)
    E <- distrib_expected_hessian(d, y, th)
    dE <- distrib_dexpected_hessian(d, y, th)
    for (a in seq_len(P)) {
      expect_agrees_on_scale(pr$grad[, a], g[[d@params[a]]],
                             label = paste(cls, "gradient", a))
      for (b in seq_len(P)) {
        key <- hess_names(d@params)[match(TRUE, hess_names(d@params) %in% c(
          paste(d@params[a], d@params[b], sep = "_"),
          paste(d@params[b], d@params[a], sep = "_")))]
        expect_agrees_on_scale(pr$hess[, a, b], h[[key]],
                               label = paste(cls, "Hessian", a, b))
        expect_agrees_on_scale(pr$info[, a, b], E[[key]],
                               label = paste(cls, "expected", a, b))
        for (c in seq_len(P)) {
          expect_agrees_on_scale(pr$dinfo[, a, b, c],
                                 dE[[dexpected_key(d@params, a, b, c)]],
                                 label = paste(cls, "dexpected", a, b, c))
        }
      }
    }
  }
})

test_that("a family the joint registry does not cover answers 0 and NaN", {
  pr <- d7_joint_probe("Gamma1Distrib", c(1, 2), cbind(c(1, 1), c(1, 1)))
  expect_identical(pr$covered, 0L)
  expect_identical(dim(pr$grad), c(2L, 0L))
  pr <- d7_joint_probe("NoSuchDistrib", c(1, 2), cbind(c(1, 1), c(1, 1)))
  expect_identical(pr$covered, -1L)
})
