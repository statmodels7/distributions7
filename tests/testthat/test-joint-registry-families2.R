# The second half of the families of test-joint-registry.R.

test_that("every joint entry is the generic's component (part 2 of 2)", {
  for (cls in test_part(names(joint_families), 2, 2)) {
    joint_twin(cls, 11)
    # the skew t below nu = 0.2 is beyond its quadrature on both routes, the R
    # one more (the identity d_sigma E_sigma,sigma = -2 E_sigma,sigma/sigma
    # fails by 4e-6 in C and 1e-4 in R at nu = 0.083, by 2.3e-2 and 1.0e-1 at
    # nu = 0.019): it is read to nu = exp(-1.5)
    joint_twin(cls, 12, if (cls == "SkewTDistrib") c(-1.5, 1.5) else c(-4, 4))
  }
})
