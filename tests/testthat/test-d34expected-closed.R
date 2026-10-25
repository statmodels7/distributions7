# The closed third and fourth derivatives of the expected information against
# references computed with mpmath at 30 to 250 digits (misure_d34e, private),
# normwise: the error of each component over the largest component of the
# same order at the same point. The stencil these kernels replace was off by
# up to 0.6 at order 3 and above 1 at order 4 on the same points.

# Two families are allowed a wider tolerance, each at the points named: negbin1
# at mu = 2242, theta = 0.028 (2.2e-9; the subtraction of mu q' on the route
# for small theta, which the second derivative carries too), and pig1 at
# mu = 1e6 with (1 + mu) sigma between 0.1 and 0.3 (1.1e-8; the sum over the
# support, above the cut of the series).
d34_tolerance <- c(negbin1 = 5e-9, pig1 = 2e-8)

test_that("the closed third and fourth derivatives match the high-precision references", {
  ref <- readRDS(test_path("fixtures", "d34expected_ref.rds"))
  ns <- asNamespace("distributions7")
  set.seed(1)
  for (fam in unique(ref$family)) {
    rf <- ref[ref$family == fam, ]
    ctor <- get(paste0(fam, "_distrib"), envir = ns)
    for (pt in unique(rf$point)) {
      rp <- rf[rf$point == pt, ]
      d <- if (is.na(rp$size[1])) ctor() else ctor(size = rp$size[1])
      th <- as.list(rp$theta[[1]])[d@params]
      # the expected information does not depend on y; any point of the
      # support serves
      y <- distrib_rng(d, 1, th)
      got <- list(`3` = distrib_d3expected_hessian(d, y, th),
                  `4` = distrib_d4expected_hessian(d, y, th))
      tol <- if (fam %in% names(d34_tolerance)) d34_tolerance[[fam]] else 1e-9
      for (j in seq_len(nrow(rp))) {
        val <- got[[as.character(rp$order[j])]][[rp$component[j]]][1]
        # an order that vanishes at the point (odd orders of the bernoulli at
        # mu = 1/2) is compared absolutely
        err <- abs(val - rp$value[j]) / max(rp$scale[j], 1e-300)
        expect_lt(err, tol, label = paste(fam, pt, rp$component[j]))
      }
    }
  }
})

test_that("every shipped univariate family registers its own third and fourth derivatives", {
  ref <- readRDS(test_path("fixtures", "d34expected_ref.rds"))
  ns <- asNamespace("distributions7")
  for (fam in unique(ref$family)) {
    ctor <- get(paste0(fam, "_distrib"), envir = ns)
    rp <- ref[ref$family == fam, ][1, ]
    d <- if (is.na(rp$size)) ctor() else ctor(size = rp$size)
    for (g in list(distrib_d3expected_hessian, distrib_d4expected_hessian)) {
      m <- S7::method(g, S7::S7_class(d))
      expect_false(identical(attr(m, "signature")[[1]], distrib), label = fam)
    }
  }
})

test_that("the closed derivatives agree with the stencil at random parameters", {
  # the stencil is good to about 1e-8 at order 3 and 1e-5 at order 4 at
  # ordinary parameters; this guards the regions the references do not sample
  ref <- readRDS(test_path("fixtures", "d34expected_ref.rds"))
  ns <- asNamespace("distributions7")
  st <- get("numerical_dexpected_higher", envir = ns)
  set.seed(20261024)
  for (fam in unique(ref$family)) {
    ctor <- get(paste0(fam, "_distrib"), envir = ns)
    rp <- ref[ref$family == fam, ][1, ]
    d <- if (is.na(rp$size)) ctor() else ctor(size = rp$size)
    for (r in 1:3) {
      th <- generate_random_theta(d)
      y <- distrib_rng(d, 1, th)
      for (k in 3:4) {
        g <- if (k == 3L) distrib_d3expected_hessian else distrib_d4expected_hessian
        a <- unlist(g(d, y, th))
        s <- unlist(st(d, y, th, "parameter", k))[names(a)]
        err <- max(abs(a - s)) / max(abs(c(a, s)))
        expect_lt(err, if (k == 3L) 1e-6 else 1e-3,
                  label = paste(fam, k, paste(signif(unlist(th), 3), collapse = ",")))
      }
    }
  }
})
