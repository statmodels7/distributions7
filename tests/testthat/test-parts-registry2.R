# The second half of the families of test-parts-registry.R.

test_that("every parts entry is its source's component (part 2 of 2)", {
  for (cls in test_part(names(joint_families), 2, 2)) parts_twin(cls, 21)
})
