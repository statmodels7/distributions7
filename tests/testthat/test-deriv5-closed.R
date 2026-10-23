# The closed fifth derivatives against references computed with mpmath at 60
# digits (misure_gas_congiunto/c/fase4/D5, private), normwise: the error of
# each component over the largest component at the same point. The stencil
# these kernels replace was off by up to 0.43 on the same points.

test_that("the closed fifth derivatives match the 60-digit references", {
  ref <- readRDS(test_path("fixtures", "deriv5_ref.rds"))
  for (fam in unique(ref$family)) {
    rf <- ref[ref$family == fam, ]
    ctor <- get(paste0(fam, "_distrib"), envir = asNamespace("distributions7"))
    for (pt in unique(rf$point)) {
      rp <- rf[rf$point == pt, ]
      d <- if (is.na(rp$size[1])) ctor() else ctor(size = rp$size[1])
      th <- as.list(rp$theta[[1]])[d@params]
      got <- distrib_deriv5(d, rp$y[1], th)
      for (j in seq_len(nrow(rp))) {
        err <- abs(got[[rp$component[j]]][1] - rp$value[j]) / rp$scale[j]
        expect_lt(err, 1e-9, label = paste(fam, pt, rp$component[j]))
      }
    }
  }
})

test_that("the closed fifth derivatives are registered for the families", {
  ex5 <- get("has_exact_deriv5", envir = asNamespace("distributions7"))
  ref <- readRDS(test_path("fixtures", "deriv5_ref.rds"))
  for (fam in unique(ref$family)) {
    ctor <- get(paste0(fam, "_distrib"), envir = asNamespace("distributions7"))
    rp <- ref[ref$family == fam, ][1, ]
    d <- if (is.na(rp$size)) ctor() else ctor(size = rp$size)
    expect_true(ex5(d), label = fam)
  }
})
