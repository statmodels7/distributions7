# Third and fourth derivatives of the two simplex-valued families.
#
# The licence for the orders where nothing independent exists is that the SAME
# assembly, run at orders one and two, reproduces the hand-written gradient and
# Hessian: those were derived separately and are already validated by
# check_distrib(). Above them the references are single stencils applied to the
# analytic Hessian -- one differentiation of an analytic quantity, never a
# difference of a difference.

mvh_cases <- function() {
  list(
    list(name = "dirichlet", d = dirichlet_distrib(3),
         th = list(mean_alr1 = 0.3, mean_alr2 = -0.4, phi = 6)),
    list(name = "dirichlet4", d = dirichlet_distrib(4),
         th = list(mean_alr1 = 0.2, mean_alr2 = -0.3, mean_alr3 = 0.5,
                   phi = 9)),
    list(name = "multinomial", d = multinomial_distrib(3, size = 7),
         th = list(probs_alr1 = 0.2, probs_alr2 = -0.5)),
    list(name = "mvstudent_t", d = mvstudent_t1_distrib(2),
         th = list(mu1 = 0.2, mu2 = -0.3, sigma_log_L1 = 0.1,
                   sigma_log_L2 = -0.2, sigma_L2.1 = 0.3, nu = 7))
  )
}

mvh_assembly <- function(name) {
  switch(name,
    multinomial = distributions7:::multinomial_higher,
    mvstudent_t = distributions7:::mvt_higher,
    distributions7:::dirichlet_higher
  )
}

# the pair -> Hessian component name lookup, built from the package's own two
# enumerations rather than by taking a name apart
hess_lookup <- function(params) {
  hn <- hess_names(params)
  hp <- distributions7:::hess_pairs(params)
  stats::setNames(hn, vapply(hp, function(p) paste(sort(p), collapse = ","),
                             character(1)))
}

test_that("the same assembly reproduces the hand-written score and information", {
  for (cs in mvh_cases()) {
    set.seed(7)
    y <- distrib_rng(cs$d, 40, cs$th)
    asm <- mvh_assembly(cs$name)
    g <- distrib_gradient(cs$d, y, cs$th)
    expect_equal(asm(cs$d, y, cs$th, 1L), g, tolerance = 1e-12,
                 info = cs$name)
    h <- distrib_hessian(cs$d, y, cs$th)
    a2 <- asm(cs$d, y, cs$th, 2L)
    expect_equal(a2, h[names(a2)], tolerance = 1e-12, info = cs$name)
  }
})

test_that("the third derivative agrees with one stencil on the Hessian", {
  skip_if_not_installed("numDeriv")
  for (cs in mvh_cases()) {
    set.seed(7)
    y <- distrib_rng(cs$d, 30, cs$th)
    params <- cs$d@params
    lookup <- hess_lookup(params)
    got <- distrib_deriv3(cs$d, y, cs$th)
    idx <- deriv_indices(params, 3L)
    for (i in seq_along(idx)) {
      t <- idx[[i]]
      from <- lookup[[paste(sort(params[t[1:2]]), collapse = ",")]]
      ref <- numDeriv::grad(function(v) {
        t2 <- cs$th
        t2[[params[t[3]]]] <- v
        sum(distrib_hessian(cs$d, y, t2)[[from]])
      }, cs$th[[params[t[3]]]])
      expect_lt(abs(sum(got[[i]]) - ref) / max(1, abs(ref)), 1e-6,
        label = sprintf("%s %s", cs$name, names(got)[i])
      )
    }
  }
})

test_that("the fourth derivative agrees with one stencil on the Hessian", {
  # a second-order stencil, in one variable or mixed across two, applied to the
  # analytic Hessian: a single stencil on an analytic quantity in both cases
  for (cs in mvh_cases()) {
    set.seed(3)
    y <- distrib_rng(cs$d, 20, cs$th)
    params <- cs$d@params
    lookup <- hess_lookup(params)
    got <- distrib_deriv4(cs$d, y, cs$th)
    idx <- deriv_indices(params, 4L)
    hess_at <- function(nm, shifts) {
      t2 <- cs$th
      for (p in names(shifts)) t2[[p]] <- t2[[p]] + shifts[[p]]
      sum(distrib_hessian(cs$d, y, t2)[[nm]])
    }
    for (i in seq_along(idx)) {
      t <- idx[[i]]
      from <- lookup[[paste(sort(params[t[1:2]]), collapse = ",")]]
      a <- params[t[3]]
      b <- params[t[4]]
      h <- 1e-3 * pmax(1, abs(unlist(cs$th[c(a, b)])))
      ref <- if (a == b) {
        s <- stats::setNames(list(0), a)
        (hess_at(from, stats::setNames(list(h[[1]]), a)) -
           2 * hess_at(from, s) +
           hess_at(from, stats::setNames(list(-h[[1]]), a))) / h[[1]]^2
      } else {
        hp <- stats::setNames(list(h[[1]], h[[2]]), c(a, b))
        (hess_at(from, hp) -
           hess_at(from, stats::setNames(list(h[[1]], -h[[2]]), c(a, b))) -
           hess_at(from, stats::setNames(list(-h[[1]], h[[2]]), c(a, b))) +
           hess_at(from, stats::setNames(list(-h[[1]], -h[[2]]), c(a, b)))) /
          (4 * h[[1]] * h[[2]])
      }
      expect_lt(abs(sum(got[[i]]) - ref) / max(1, abs(ref)), 1e-4,
        label = sprintf("%s %s", cs$name, names(got)[i])
      )
    }
  }
})

test_that("the components are per observation, named and finite", {
  for (cs in mvh_cases()) {
    set.seed(5)
    y <- distrib_rng(cs$d, 12, cs$th)
    for (o in 3:4) {
      v <- if (o == 3L) distrib_deriv3(cs$d, y, cs$th) else
                        distrib_deriv4(cs$d, y, cs$th)
      expect_named(v, deriv_names(cs$d@params, o))
      expect_true(all(vapply(v, length, integer(1)) == 12L), info = cs$name)
      expect_true(all(is.finite(unlist(v))), info = cs$name)
    }
  }
})

test_that("neither family is on the numerical fallback any more", {
  # the census this batch was built from, run as a test: a method is the base
  # fallback when the class it was registered on is the abstract one
  for (cs in mvh_cases()) {
    for (gen in list(distrib_deriv3, distrib_deriv4)) {
      m <- S7::method(gen, S7::S7_class(cs$d))
      owner <- attr(attr(m, "signature")[[1L]], "name")
      expect_false(owner %in% c("distrib", "multivariate_distrib"),
                   info = sprintf("%s: %s", cs$name, owner))
    }
  }
})

test_that("the expected derivatives are the observed ones averaged", {
  # sampling is the multivariate route throughout; what is checked here is that
  # the shape is right and the answer is constant across observations
  set.seed(2)
  d <- multinomial_distrib(3, size = 5)
  th <- list(probs_alr1 = 0.1, probs_alr2 = -0.2)
  y <- distrib_rng(d, 8, th)
  e <- distrib_deriv3(d, y, th, expected = TRUE, nsim = 2000)
  expect_named(e, deriv_names(d@params, 3L))
  expect_true(all(vapply(e, function(v) length(unique(v)) == 1L, logical(1))))
  expect_true(all(is.finite(unlist(e))))
})
