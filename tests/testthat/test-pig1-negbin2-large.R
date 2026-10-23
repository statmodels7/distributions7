# pig1 at large counts and negbin2's expected information at large theta,
# against references computed with mpmath at 60 and 50 digits
# (misure_gas_congiunto/c/fase3/pig1_ref.py and c/nb_ref.py, private).

test_that("pig1's derivatives at large counts match a 60-digit reference", {
  d <- pig1_distrib()
  ref <- list(
    list(y = 200, mu = 180, sigma = 0.2,
         mu_d = c(0.005509209504699214891333061, -0.0001822055907330399097382815,
                  0.000002944541199898687437588887, -6.406559134852970783893632e-8),
         sigma_d = c(-2.318241014068995996705304, 10.71831633806700263889947,
                     -98.80173086812238406405235, 1361.347221479428568446679)),
    list(y = 2000, mu = 3000, sigma = 5,
         mu_d = c(0.0001388488245480520407502676, -7.036467307698831227025076e-8,
                  5.184716067000645331044324e-11, -5.678353961585476460947576e-14),
         sigma_d = c(-0.09664818201514647986512393, 0.01865941256093690047180392,
                     -0.007195731379217394939854889, 0.004156652169281889763596144)),
    list(y = 20000, mu = 20000, sigma = 0.001,
         mu_d = c(0.00002267562432302919667507327, -0.00000238213501459933595437736,
                  3.5698892105089657309105e-10, -7.055288827498593396432438e-14),
         sigma_d = c(-476.1881107836131301765386, 453510.4075830995766121855,
                     -863825713.3919169490543939, 2468063995853.034974866542)))
  gens <- list(distrib_gradient, distrib_hessian, distrib_deriv3, distrib_deriv4)
  for (r in ref) {
    th <- list(mu = r$mu, sigma = r$sigma)
    for (o in 1:4) {
      g <- gens[[o]](d, r$y, th)
      expect_lt(abs(g[[paste(rep("mu", o), collapse = "_")]] / r$mu_d[o] - 1), 1e-10)
      expect_lt(abs(g[[paste(rep("sigma", o), collapse = "_")]] / r$sigma_d[o] - 1), 1e-10)
    }
  }
})

test_that("pig1's log mass is finite at counts past 1e10", {
  d <- pig1_distrib()
  for (y in c(1e10, 1e12)) {
    expect_true(is.finite(distrib_pdf(d, y, list(mu = y, sigma = 0.5), log = TRUE)))
  }
})

test_that("pig1's second derivative of the expected information holds at mu = 1000", {
  # the support recurrence's derivatives in w at large counts; the surface
  # was 1.2e-7 off normwise when they were cumulants of raw moments
  d <- pig1_distrib()
  ref <- c(mu_mu_mu_mu = -1.178757192457368428382951e-10,
           mu_mu_sigma_sigma = -0.01504989251584675153336061,
           mu_mu_mu_sigma = -0.000000751796429474595467175609,
           sigma_sigma_mu_mu = 0.00001693328775555317036464365,
           sigma_sigma_sigma_sigma = -440413.1782906228045864754,
           sigma_sigma_mu_sigma = 0.4890889582143867543779087,
           mu_sigma_mu_mu = 1.630267500463967619764873e-8,
           mu_sigma_sigma_sigma = 6.626795601512428664194837,
           mu_sigma_mu_sigma = 0.0001539831498490684666728612)
  out <- distrib_d2expected_hessian(d, 1, list(mu = 1000, sigma = 0.05))
  got <- vapply(names(ref), function(k) out[[k]][1], numeric(1))
  expect_lt(max(abs(got - ref)) / max(abs(ref)), 1e-10)
})

test_that("negbin2's expected information derivatives hold at theta = 1e5", {
  d <- negbin2_distrib()
  p <- d@params
  pairs <- list(c(1, 1), c(1, 2), c(2, 2))
  cd <- list(c(1, 1), c(1, 2), c(2, 2))
  ref <- list(
    list(mu = 1000,
         dE = list(c(0.0000009999019703950593079109901, -9.802960494069208901088138e-11),
                   c(-8.541496564010769411407586e-48, -2.118672982984406990451184e-48),
                   c(-9.705806343135834469896102e-18, 1.950862173586142737785736e-19)),
         d2E = list(c(-1.999998058819704144731109e-9, 1.941180295855288891342647e-15,
                      1.941180295855288891323609e-15),
                    c(1.477847346287571117916828e-43, -6.702719398327210715098423e-46,
                      2.273858341449799442255117e-44),
                    c(-9.417516968492367288916948e-21, 3.853484086267274613164082e-22,
                      -9.706046583556976460918463e-24))))
  for (r in ref) {
    th <- stats::setNames(list(r$mu, 1e5), p)
    dE <- distrib_dexpected_hessian(d, 0, th)
    d2E <- distrib_d2expected_hessian(d, 0, th)
    s1 <- max(abs(unlist(r$dE)))
    s2 <- max(abs(unlist(r$d2E)))
    for (i in 1:3) {
      a <- pairs[[i]][1]; b <- pairs[[i]][2]
      for (cc in 1:2) {
        ex <- r$dE[[i]][cc]
        got <- dE[[dexpected_key(p, a, b, cc)]]
        expect_lt(abs(got - ex) / max(abs(ex), 1e-16 * s1), 1e-10)
      }
      for (j in 1:3) {
        ex <- r$d2E[[i]][j]
        got <- d2E[[d2expected_key(p, a, b, cd[[j]][1], cd[[j]][2])]]
        expect_lt(abs(got - ex) / max(abs(ex), 1e-16 * s2), 5e-9)
      }
    }
  }
})

test_that("the PIG log mass holds at counts past 1e10", {
  # 80-digit references (uniform expansion of the Bessel function); the
  # assembly in terms of size y log y read -28.203125 at 1e12 and 0 at 9e15
  ref <- c(`1e10` = -23.598215872927656927, `1e12` = -28.203386058853873295,
           `9e15` = -37.3083659151716048)
  for (k in names(ref)) {
    y <- as.numeric(k)
    a <- sqrt(1 + 2 * 0.5 * y) / 0.5
    l1 <- distrib_pdf(pig1_distrib(), y, list(mu = y, sigma = 0.5), log = TRUE)
    l2 <- distrib_pdf(pig2_distrib(), y, list(mu = y, alpha = a), log = TRUE)
    expect_lt(abs(l1 / ref[[k]] - 1), 1e-13)
    expect_lt(abs(l2 / ref[[k]] - 1), 1e-13)
  }
})

test_that("pig2's second derivative of the expected information holds at mu = 1e4", {
  # the support recurrence's derivatives in alpha; 2.5e-4 off normwise when
  # they were cumulants of raw moments
  ref <- c(mu_mu_mu_mu = -2.230601660926372847108061e-13,
           mu_mu_alpha_alpha = -8.340737995891702541153359e-13,
           mu_mu_mu_alpha = 5.439611736451110352926104e-13,
           alpha_alpha_mu_mu = 1.522623076841447683587475e-15,
           alpha_alpha_alpha_alpha = -7.500024289386655434932021e-13,
           alpha_alpha_mu_alpha = 5.738850416912363975907761e-16)
  out <- distrib_d2expected_hessian(pig2_distrib(), 1, list(mu = 1e4, alpha = 2000))
  got <- vapply(names(ref), function(k) out[[k]][1], numeric(1))
  expect_lt(max(abs(got - ref)) / max(abs(ref)), 1e-9)
})

test_that("negbin1's second derivative of the expected information holds at theta = 0.05", {
  # the route chosen per order: on the W sums the surface at mu = 100 was
  # 2.9e-10 off relative to its largest component
  d <- negbin1_distrib()
  p <- d@params
  ex <- list(c(-0.00000190482986831101792745884, -0.00009061665448909045520236804,
               -0.01735068754801322028587279),
             c(4.528782139868751923181519e-8, -0.00004095420407502477399904835,
               -0.01610177533815273484199777),
             c(4.241652333129850450905146e-8, -0.00003562188883104009475132075,
               -2.492566450997503632716402))
  pairs <- list(c(1, 1), c(1, 2), c(2, 2))
  out <- distrib_d2expected_hessian(d, 0, stats::setNames(list(100, 0.05), p))
  err <- 0
  for (i in 1:3) for (j in 1:3) {
    k <- d2expected_key(p, pairs[[i]][1], pairs[[i]][2], pairs[[j]][1], pairs[[j]][2])
    err <- max(err, abs(out[[k]] - ex[[i]][j]))
  }
  expect_lt(err / max(abs(unlist(ex))), 1e-11)
})
