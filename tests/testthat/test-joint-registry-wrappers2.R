# Part 2 of 3 of the wrappers of the joint registry (helper-joint-registry.R).

test_that("every joint entry of a wrapper is the generic's component (part 2 of 3)", {
  for (nm in test_part(names(joint_wrappers), 2, 3)) {
    d <- joint_wrappers[[nm]][[1]]()
    full <- joint_wrappers[[nm]][[2]]
    rg <- if (length(joint_wrappers[[nm]]) > 2L) joint_wrappers[[nm]][[3]] else
      list(c(-1.5, 1.5), c(-4, 4))
    # a folded family is not read at y = 1e40: its R Hessian is R2 - R1^2 of
    # log_deriv(), two terms of order 1e160 whose difference is of order
    # 1e80, so the R value there is the rounding of the subtraction (it
    # differs between platforms by 1e144), where the registry has the form
    # without the cancellation
    far <- !grepl("^fo_|folded", nm)
    joint_twin(nm, 11, rg[[1]], d = d, full = full, far = far)
    joint_twin(nm, 12, rg[[2]], d = d, full = full, far = far)
  }
})
