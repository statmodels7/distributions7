# The compiled expected Hessian by simulation and by integration
# (d7_expect(), src/d7_expect.cpp) against the R route,
# distrib_expected_hessian(method = expect_mc() / expect_integrate()), with
# identical(): the same doubles and, for the simulation, the same state of
# the random stream afterwards. The families are those of
# helper-joint-registry.R.

expect_twin_one <- function(d, th, method, nsim = 50, seed = 1) {
  route <- distrib_scalar_route(d)
  id <- match(route$name, d7_parts_classes_covered()) - 1L
  n <- max(lengths(th))
  M <- do.call(cbind, c(lapply(th, rep_len, n), lapply(route$constants, rep, length.out = n)))
  mc <- method %in% 0:1
  ident <- if (method %in% c(1L, 3L)) "score" else "hessian"
  em <- if (mc) expect_mc(nsim = nsim, identity = ident) else expect_integrate(identity = ident)
  tol <- if (mc) rep(0, 6) else c(em@quad_atol, em@quad_rtol, em@max_panels,
                                  em@series_atol, em@series_rtol, em@max_terms)
  y <- rep(distrib_rng(d, 1, lapply(th, `[`, 1)), n)
  set.seed(seed)
  ref <- tryCatch(distrib_expected_hessian(d, y, th,
                                           method = em, scale = "parameter"),
                  error = function(e) e)
  after_r <- stats::runif(1)
  set.seed(seed)
  got <- d7_expect_probe(id, method, M, d@bounds, nsim, tol)
  after_c <- stats::runif(1)
  lab <- sprintf("%s method %d", route$name, method)
  if (is.integer(got)) {
    # handed back: R either signals an error or, for code 2, rescues
    if (got == 1L) expect_true(inherits(ref, "error") || !mc, label = paste(lab, "hand-back"))
    return(invisible(got))
  }
  expect_false(inherits(ref, "error"), label = paste(lab, "R route"))
  if (inherits(ref, "error")) return(invisible(NULL))
  ref <- do.call(cbind, ref[hess_names(d@params)])
  expect_identical(unname(got), unname(ref), label = lab)
  if (mc) expect_identical(after_c, after_r, label = paste(lab, "stream"))
  invisible(got)
}

expect_twin_family <- function(cls, methods = 0:3, nsim = 50) {
  d <- joint_families[[cls]]()
  set.seed(19)
  n <- 4
  th <- lapply(d@params, function(p)
    linkfunctions7::linkinv(d@link_params[[p]], stats::runif(n, -1, 1)))
  names(th) <- d@params
  # a repeated setting, and a scalar setting
  th <- lapply(th, function(v) { v[3] <- v[1]; v })
  th1 <- lapply(th, `[`, 2)
  for (m in methods) {
    expect_twin_one(d, th, m, nsim)
    expect_twin_one(d, th1, m, nsim)
  }
}

test_that("the compiled expected Hessian is the R route's (part 1 of 3)", {
  for (cls in test_part(names(joint_families), 1, 3)) expect_twin_family(cls)
})

test_that("the compiled expected Hessian is the R route's (part 2 of 3)", {
  for (cls in test_part(names(joint_families), 2, 3)) expect_twin_family(cls)
})

test_that("the compiled expected Hessian is the R route's (part 3 of 3)", {
  for (cls in test_part(names(joint_families), 3, 3)) expect_twin_family(cls)
})
