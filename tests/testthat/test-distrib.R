test_that("distrib class hierarchy is strictly enforced", {
  # Check Continuous inheritance
  d_cont <- gaussian_distrib()
  expect_true(S7::S7_inherits(d_cont, distributions7::distrib))
  expect_true(S7::S7_inherits(d_cont, distributions7::continuous_distrib))
  expect_false(S7::S7_inherits(d_cont, distributions7::discrete_distrib))
  
  # Check Discrete inheritance
  d_disc <- poisson_distrib()
  expect_true(S7::S7_inherits(d_disc, distributions7::distrib))
  expect_true(S7::S7_inherits(d_disc, distributions7::discrete_distrib))
  expect_false(S7::S7_inherits(d_disc, distributions7::continuous_distrib))
})
