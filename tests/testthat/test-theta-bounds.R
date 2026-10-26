# Parameter values are validated against params_bounds (open intervals) by every
# generic, so out-of-domain input errors clearly instead of yielding NaN.

test_that("out-of-domain parameters are rejected with an informative message", {
  d <- gaussian1_distrib()

  expect_error(distrib_pdf(d, 0, list(mu = 0, sigma = -1)), "outside its domain")
  expect_error(distrib_pdf(d, 0, list(mu = 0, sigma = -1)), "sigma")
  expect_error(distrib_pdf(d, 0, list(mu = 0, sigma = 0)), "\\(0, Inf\\)")
  expect_error(distrib_pdf(d, 0, list(mu = NA, sigma = 1)), "mu")
  expect_error(distrib_pdf(d, 0, list(mu = Inf, sigma = 1)), "mu")

  # bounded parameters: open interval on both sides
  b <- bernoulli_distrib()
  expect_error(distrib_pdf(b, 1, list(mu = 0)), "outside its domain")
  expect_error(distrib_pdf(b, 1, list(mu = 1)), "outside its domain")
  expect_error(distrib_pdf(b, 1, list(mu = 1.2)), "outside its domain")
  expect_silent(distrib_pdf(b, 1, list(mu = 0.5)))
})

test_that("the check applies to every generic, not just the density", {
  d <- gaussian1_distrib()
  bad <- list(mu = 0, sigma = -2)

  expect_error(distrib_cdf(d, 0, bad), "outside its domain")
  expect_error(distrib_quantile(d, 0.5, bad), "outside its domain")
  expect_error(distrib_rng(d, 1, bad), "outside its domain")
  expect_error(distrib_gradient(d, 0, bad), "outside its domain")
  expect_error(distrib_hessian(d, 0, bad), "outside its domain")
  expect_error(distrib_expected_hessian(d, 0, bad), "outside its domain")
  expect_error(distrib_deriv3(d, 0, bad), "outside its domain")
  expect_error(distrib_deriv4(d, 0, bad), "outside its domain")
  expect_error(distrib_grad_y(d, 0, bad), "outside its domain")
})

test_that("check_theta_bounds is usable directly and reports every offender", {
  d <- student_t1_distrib()
  expect_invisible(check_theta_bounds(d, list(mu = 0, sigma = 1, nu = 5)))

  err <- tryCatch(
    check_theta_bounds(d, list(mu = 0, sigma = -1, nu = -2)),
    error = function(e) conditionMessage(e)
  )
  expect_match(err, "sigma")
  expect_match(err, "nu")
  expect_match(err, "student t")
})

test_that("vectorized parameters are validated element-wise", {
  d <- gaussian1_distrib()
  expect_silent(distrib_pdf(d, c(0, 1), list(mu = c(0, 1), sigma = c(1, 2))))
  expect_error(
    distrib_pdf(d, c(0, 1), list(mu = c(0, 1), sigma = c(1, -2))),
    "outside its domain"
  )
})

# The compiled pass and the R expression it replaces, on values at and
# around every kind of bound, and the message against the function as it
# was written before the pass (kept here as the reference).
test_that("the compiled bound test agrees with the R expression", {
  bounds <- list(c(-Inf, Inf), c(0, Inf), c(0, 1), c(2, 5), c(-Inf, 0))
  special <- c(NA, NaN, Inf, -Inf, 0, -0, 1, 2, 5, -1, .Machine$double.xmin,
               4.9e-324, -4.9e-324, 1 - .Machine$double.eps / 2,
               1 + .Machine$double.eps, .Machine$double.xmax, 3)
  set.seed(7)
  for (b in bounds) {
    for (v in c(as.list(special), list(special, runif(50, -10, 10),
                                       numeric(0)))) {
      ref <- all(is.finite(v) & v > b[1] & v < b[2])
      expect_identical(bounds_ok_cpp(v, b[1], b[2]), ref,
                       info = paste(format(b), collapse = ","))
    }
  }
})

test_that("check_bounds_fast() reports as it did before the compiled pass", {
  reference <- function(params, param_bounds, theta, distrib_name) {
    problems <- character()
    for (i in seq_along(params)) {
      p <- params[i]
      b <- param_bounds[[p]]
      if (is.null(b) || i > length(theta)) next
      v <- theta[[i]]
      if (!is.numeric(v)) {
        problems <- c(problems, sprintf("  '%s' must be numeric (got %s)", p,
                                        class(v)[1]))
        next
      }
      ok <- is.finite(v) & v > b[1] & v < b[2]
      if (!all(ok)) {
        offend <- unique(v[!ok])
        shown <- paste(format(utils::head(offend, 3), trim = TRUE),
                       collapse = ", ")
        if (length(offend) > 3) shown <- paste0(shown, ", ...")
        problems <- c(problems, sprintf(
          "  '%s' = %s is outside its domain (%s, %s)",
          p, shown, format(b[1]), format(b[2])))
      }
    }
    if (length(problems) > 0) {
      stop("Invalid parameter value(s) for the '", distrib_name,
           "' distribution:\n", paste(problems, collapse = "\n"), call. = FALSE)
    }
    invisible(NULL)
  }
  msg <- function(f, ...) tryCatch({ f(...); "no error" },
                                   error = function(e) conditionMessage(e))
  d <- student_t1_distrib()
  cases <- list(
    list(mu = 0, sigma = 1, nu = 5),
    list(mu = c(0, NA, 1), sigma = c(1, -1, 0, -2, -3), nu = 5),
    list(mu = 0, sigma = 1L, nu = 0L),
    list(mu = Sys.Date(), sigma = 1, nu = 5),
    list(mu = "a", sigma = c(1, 2), nu = c(Inf, 3)),
    list(mu = 1, sigma = numeric(0), nu = NaN))
  for (th in cases) {
    expect_identical(
      msg(check_bounds_fast, d@params, d@params_bounds, th, d@distrib_name),
      msg(reference, d@params, d@params_bounds, th, d@distrib_name))
  }
})
