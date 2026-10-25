# The expected information of pig1 and pig2 at large mu, on the routes that
# replaced the sum over the support (closed forms, series and expansions; the
# sum cost 0.2 to 1.8 s per observation at mu = 1e6). References: mpmath at
# 28 to 70 digits (misure_open5/p6_pig_fast, p8_pig1_fixsig_impl).

rel_order <- function(got, ref) max(abs(got - ref)) / max(abs(ref))

test_that("pig1 and pig2 orders 0-2 at large mu match the references", {
  ker <- list(pig1 = list(pig1_expected_cpp, pig1_dexpected1_cpp, pig1_dexpected2_cpp),
              pig2 = list(pig2_expected_cpp, pig2_dexpected1_cpp, pig2_dexpected2_cpp))
  expect_lt(rel_order(unlist(ker$pig1[[1]](1, 100000, 8e-7, 1L), use.names = FALSE),
                      c(-9.2592592607645993e-06, -4.2866901595747671e+09,  2.5402608353035658e-03)), 1e-12)
  expect_lt(rel_order(unlist(ker$pig1[[2]](1, 100000, 8e-7, 1L), use.names = FALSE),
                      c( 9.9451303129360332e-11,  8.5733881333261963e-01, -7.9383146360338389e+04,  7.9383703130764562e+14,  4.5160189817012703e-08,  5.6450208071781281e+03)), 1e-12)
  expect_lt(rel_order(unlist(ker$pig1[[3]](1, 100000, 8e-7, 1L), use.names = FALSE),
                      c(-1.9991871159243720e-15, -1.5876647111192072e+05, -1.2701316992038133e-06, -6.1742434084062048e-01, -2.2050943926376949e+20,  2.2050987522344051e+10,  2.9897710595442298e-13,  4.6715034246835699e+09,  9.3822320092123976e-02)), 1e-12)
  expect_lt(rel_order(unlist(ker$pig1[[1]](1, 1000000, 5e-7, 1L), use.names = FALSE),
                      c(-6.6666667283950126e-07, -2.2222204526766513e+11,  3.7037007544610860e-02)), 1e-12)
  expect_lt(rel_order(unlist(ker$pig1[[2]](1, 1000000, 5e-7, 1L), use.names = FALSE),
                      c( 8.8888888477366665e-13,  4.4444441152267139e-01, -2.9629603155038151e+05,  2.9629647187866080e+17,  3.7037002743526636e-08,  7.4073946502259379e+04)), 1e-12)
  expect_lt(rel_order(unlist(ker$pig1[[3]](1, 1000000, 5e-7, 1L), use.names = FALSE),
                      c(-1.9259259190672213e-18, -5.9259268038389529e+05, -2.9629630727021227e-07,  2.4691257430596083e-08, -5.9259319615952447e+23,  5.9259264197448157e+11, -2.4691339277560661e-14, -9.8765631457752686e+10,  2.4691258344995960e-02)), 1e-12)
  expect_lt(rel_order(unlist(ker$pig2[[1]](1, 1000000, 1000, 1L), use.names = FALSE),
                      c(-4.9999962500031249e-13, -1.9999520044378176e-06,  0.0000000000000000e+00)), 1e-12)
  expect_lt(rel_order(unlist(ker$pig2[[2]](1, 1000000, 1000, 1L), use.names = FALSE),
                      c( 1.4999981250021875e-18, -9.9999850000187508e-16, -2.0997900007416073e-17,  3.9998500224248064e-09,  0.0000000000000000e+00,  0.0000000000000000e+00)), 1e-12)
  expect_lt(rel_order(unlist(ker$pig2[[3]](1, 1000000, 1000, 1L), use.names = FALSE),
                      c(-5.9999887500174997e-24, -9.9999550000937500e-19,  2.9999925000131251e-21,  4.7994913738579982e-23, -1.1999310143705963e-11,  2.9993373596530131e-20,  0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00)), 1e-12)
  expect_lt(rel_order(unlist(ker$pig2[[1]](1, 100000, 1000000, 1L), use.names = FALSE),
                      c(-9.0049628097900101e-06, -5.9846872498401897e-15,  0.0000000000000000e+00)), 1e-12)
  expect_lt(rel_order(unlist(ker$pig2[[2]](1, 100000, 1000000, 1L), use.names = FALSE),
                      c( 9.9901481466315840e-11, -9.8518533684157342e-13, -1.2923352369870858e-19,  2.4892718930585585e-20,  0.0000000000000000e+00,  0.0000000000000000e+00)), 1e-12)
  expect_lt(rel_order(unlist(ker$pig2[[3]](1, 100000, 1000000, 1L), use.names = FALSE),
                      c(-2.0009559224060443e-15,  1.9411077428858722e-18,  2.9262930797274458e-19, -1.5693236859568364e-24, -1.2914138847710100e-25,  5.4463274735700042e-25,  0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00)), 1e-12)
})

test_that("pig1 at large mu and moderate sigma reads the expansion at fixed sigma", {
  d <- pig1_distrib()
  e <- distrib_expected_hessian(d, 1, list(mu = 10000, sigma = 0.5))
  expect_lt(rel_order(unlist(e)[c("mu_mu", "sigma_sigma", "mu_sigma")], c(-2.4986403242764036e-08, -1.9969595214069082e+00,  9.9828010468251761e-05)), 1e-13)
  e <- distrib_expected_hessian(d, 1, list(mu = 10000, sigma = 1))
  expect_lt(rel_order(unlist(e)[c("mu_mu", "sigma_sigma", "mu_sigma")], c(-1.4989544638048999e-08, -4.9915426968720555e-01,  4.9910435925128041e-05)), 1e-13)
  e <- distrib_expected_hessian(d, 1, list(mu = 1e+06, sigma = 0.5))
  expect_lt(rel_order(unlist(e)[c("mu_mu", "sigma_sigma", "mu_sigma")], c(-2.4999863752840768e-12, -1.9999695009583118e+00,  9.9998275051365503e-07)), 1e-13)
  e <- distrib_expected_hessian(d, 1, list(mu = 1e+06, sigma = 1))
  expect_lt(rel_order(unlist(e)[c("mu_mu", "sigma_sigma", "mu_sigma")], c(-1.4999895004515714e-12, -4.9999150043207224e-01,  4.9999100044107181e-07)), 1e-13)
  # the sum over the support took about 0.3 s per row here
  th <- list(mu = rep(1e6, 200), sigma = rep(0.5, 200))
  t0 <- proc.time()[["elapsed"]]
  distrib_d2expected_hessian(d, rep(1, 200), th)
  expect_lt(proc.time()[["elapsed"]] - t0, 5)
})
