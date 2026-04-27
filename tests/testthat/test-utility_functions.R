test_that("check_params_dim validates correctly", {
  # Implicit max length
  expect_invisible(check_params_dim(list(mu = 1, sigma = 2)))
  expect_invisible(check_params_dim(list(mu = 1:5, sigma = 1)))
  
  # Should error on incompatible lengths (e.g., 2 vs 3)
  expect_error(check_params_dim(list(mu = 1:2, sigma = 1:3)), "Parameter dimension mismatch")
  
  # Explicit target length `n`
  expect_invisible(check_params_dim(list(mu = 1:5, sigma = 1), n = 5))
  
  # Should error if vector length does not match explicit n
  expect_error(check_params_dim(list(mu = 1:3, sigma = 1), n = 5), "Parameter dimension mismatch")
})

test_that("expand_params expands scalars correctly", {
  # Expands the scalar to match the maximum length
  res <- expand_params(list(mu = 1, sigma = 1:3))
  expect_equal(res$mu, c(1, 1, 1))
  expect_equal(res$sigma, 1:3)
  
  # Expands scalars based on explicit target length `n`
  res2 <- expand_params(list(mu = 2, sigma = 3), n = 2)
  expect_equal(res2$mu, c(2, 2))
  expect_equal(res2$sigma, c(3, 3))
})

test_that("transpose_params swaps rows and columns and is symmetric", {
  theta_cols <- list(mu = c(10, 20), sigma = c(1, 2))
  theta_rows <- transpose_params(theta_cols)
  
  expect_type(theta_rows, "list")
  expect_length(theta_rows, 2)
  expect_equal(theta_rows[[1]], c(mu = 10, sigma = 1))
  expect_equal(theta_rows[[2]], c(mu = 20, sigma = 2))
  
  # Inverse transformation (symmetry check)
  expect_equal(transpose_params(theta_rows), theta_cols)
})