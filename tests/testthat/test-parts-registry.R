# The parts route of the C registry (pt_parts.h): a buffer filled once per
# observation with the quantities that cost, then one entry per order over
# it. Its gradient, Hessian and expected information are compared with the
# joint entries, and its orders three and four with distrib_deriv3() and
# distrib_deriv4(), on the scale of the expected values at 1e-10 (the
# reading of test-ccallable.R, for the fused multiply-adds of arm64 macOS).
# The families and parts_twin() are in helper-parts-registry.R.

test_that("the parts route covers every class of the joint registry", {
  expect_setequal(d7_parts_classes_covered(), names(joint_families))
})

test_that("every parts entry is its source's component (part 1 of 2)", {
  for (cls in test_part(names(joint_families), 1, 2)) parts_twin(cls, 21)
})

test_that("an order above the buffer's level, or an uncovered family, is NaN", {
  d <- gamma1_distrib()
  pt <- parts_points(d, 23, c(-1, 1))
  M <- do.call(cbind, pt$th)
  pr <- d7_parts_probe("Gamma1Distrib", pt$y, M, level = 2)
  expect_true(all(is.finite(pr$grad[-(1:2), ])) && all(is.finite(pr$hess[-(1:2), , ])))
  expect_true(all(is.nan(pr$d3)) && all(is.nan(pr$d4)))
  pr1 <- d7_parts_probe("Gamma1Distrib", pt$y, M, level = 1)
  expect_true(all(is.nan(pr1$hess)) && all(is.finite(pr1$grad[-(1:2), ])))
  expect_identical(d7_parts_probe("NoSuchDistrib", 1, matrix(2))$covered, 0L)
})
