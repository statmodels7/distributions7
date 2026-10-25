# Regression tests for the precision repairs of 2026-10-25 (misure_open5):
# lognormal2's observed Hessian at small v/m^2, the folded/truncated skew t
# and folded Student t expected information at small nu, the skew t's d4E
# component d^4 E_{alpha,nu}/d nu^4 at nu = 0.05, and negbin1's expected
# information and its derivatives at large mu and small theta.
# References: mpmath, 30 to 90 digits, rounded to 17 significant digits.

# largest |got - ref| over the entries of ref, relative to the largest |ref|
rel_err <- function(got, ref) {
  got <- unlist(got)[names(ref)]
  max(abs(got - ref)) / max(abs(ref))
}

test_that("lognormal2 second derivatives are accurate at v/m^2 = 1e-8", {
  d <- lognormal2_distrib()
  cases <- list(
  list(m = 10000, v = 1, y = 9998.0001999916658, fn = "p2_",
       ref = c(mean_mean = -1.0006000599974998, mean_var = 2.0003999949982885, var_var = -3.4999999800021535)),
  list(m = 10000, v = 1, y = 9998.0001999916658, fn = "grad_y_hess_",
       ref = c(mean_mean = 0.00030010001350085005, mean_var = -1.0006001000078335, var_var = 4.00080007000241)),
  list(m = 10000, v = 1, y = 9998.0001999916658, fn = "hess_y_hess_",
       ref = c(mean_mean = -5.0024005550773403e-08, mean_var = 0.00030016004000600057, var_var = -2.0012003200503385)),
  list(m = 10000, v = 1, y = 10000, fn = "p2_",
       ref = c(mean_mean = -1.0000000199999999, mean_var = -2.499999966666667e-13, var_var = 0.5)),
  list(m = 10000, v = 1, y = 10000, fn = "grad_y_hess_",
       ref = c(mean_mean = 0.00029999999950000003, mean_var = -1, var_var = 0)),
  list(m = 10000, v = 1, y = 10000, fn = "hess_y_hess_",
       ref = c(mean_mean = -4.999999995e-08, mean_var = 0.00029999999999999997, var_var = -2)),
  list(m = 10000, v = 1, y = 10003.000450037502, fn = "p2_",
       ref = c(mean_mean = -0.99910011000374932, mean_var = -2.999099992505275, var_var = -8.4999999550031529)),
  list(m = 10000, v = 1, y = 10003.000450037502, fn = "grad_y_hess_",
       ref = c(mean_mean = 0.00029985003099647519, mean_var = -0.99910022497075179, var_var = -5.9982002549830504)),
  list(m = 10000, v = 1, y = 10003.000450037502, fn = "hess_y_hess_",
       ref = c(mean_mean = -4.9964012547240406e-08, mean_var = 0.0002997600899790033, var_var = -1.9982007198245284))
  )
  gen <- list(p2_ = distrib_hessian, grad_y_hess_ = distrib_grad_y_hess,
              hess_y_hess_ = distrib_hess_y_hess)
  for (cs in cases) {
    th <- list(mean = cs$m, var = cs$v)
    got <- gen[[cs$fn]](d, cs$y, th, scale = "parameter")
    expect_lt(rel_err(got, cs$ref), 5e-14)
  }
})

# expected Hessian (minus the information) from the stored information, upper triangle by rows
info_ref <- function(params, vals) {
  p <- length(params)
  nm <- unlist(lapply(seq_len(p), function(a) paste(params[a], params[a:p], sep = "_")))
  setNames(-vals, nm)
}
sk_par <- c("mu", "sigma", "alpha", "nu")
t1_par <- c("mu", "sigma", "nu")
ref3 <- list(
  fold_skewt_0.05 = c(0.33246553590111944, 0.017252523113298814, 0.063089970676464241, 0.1281495705280139, 0.029067257899008294, 0.0022721283667064386, -0.65246560848478607, 0.014205940938138812, 0.0079842745643933406, 383.27171921409416),
  trunc_skewt_0.05 = c(0.56465097910323592, 0.021977656628514669, 0.11480434228926814, 0.54183457630903764, 0.04689202747355159, 0.02042364745883991, -0.83044085729339745, 0.047903756715246291, -0.85323813779854407, 373.79679338171087),
  fold_t1_0.05 = c(0.25765997541848557, 0.0084528874939854282, 0.051430542887059845, 0.027106621470032358, -0.67033904429014957, 383.11619579067678),
  fold_skewt_0.1 = c(0.33077304356136655, 0.029945776205122714, 0.078790657715826218, 0.080370717840235056, 0.055360044781139581, 0.0063667023932479265, -0.61495812426943264, 0.020847636923648844, 0.0062600397244649475, 91.752187388817617),
  trunc_skewt_0.1 = c(0.54996212619917118, 0.033254347018649044, 0.14560408410306516, 0.37559145493390883, 0.086715872244148284, 0.034823156482770484, -0.71825958887473473, 0.070308910268399388, -0.59289314342317456, 87.403897164014538),
  fold_t1_0.1 = c(0.24136109247793192, 0.019767422393386625, 0.046990430584055692, 0.051748557953860339, -0.63150448080473176, 91.68850570220043)
)
th_sk <- function(nu) list(mu = 0.5, sigma = 1, alpha = 0.5, nu = nu)
th_t1 <- function(nu) list(mu = 0.5, sigma = 1, nu = nu)

test_that("folded and truncated skew t expected information at small nu", {
  for (nu in c(0.05, 0.1)) {
    key <- as.character(nu)
    d <- folded(skewt_distrib())
    got <- distrib_expected_hessian(d, 1, th_sk(nu), scale = "parameter")
    expect_lt(rel_err(got, info_ref(sk_par, ref3[[paste0("fold_skewt_", key)]])), 5e-14)
    d <- truncated(skewt_distrib(), lower = -0.5)
    got <- distrib_expected_hessian(d, 1, th_sk(nu), scale = "parameter")
    expect_lt(rel_err(got, info_ref(sk_par, ref3[[paste0("trunc_skewt_", key)]])), 5e-14)
  }
})

test_that("folded Student t expected information at small nu", {
  for (nu in c(0.05, 0.1)) {
    d <- folded(student_t1_distrib())
    got <- distrib_expected_hessian(d, 1, th_t1(nu), scale = "parameter")
    expect_lt(rel_err(got, info_ref(t1_par, ref3[[paste0("fold_t1_", nu)]])), 5e-14)
  }
})

test_that("folded skew t information derivative matches differences of the information", {
  d <- folded(skewt_distrib())
  th <- th_sk(0.05)
  P <- names(th)
  up <- unlist(distrib_dexpected_hessian(d, 1, th, scale = "parameter"))
  ev <- function(c, s) {
    t2 <- th
    t2[[c]] <- th[[c]] + s
    unlist(distrib_expected_hessian(d, 1, t2, scale = "parameter"))
  }
  err <- 0
  for (c in seq_along(P)) {
    h <- 1e-3 * min(abs(th[[c]]), if (P[c] %in% c("sigma", "nu")) th[[c]] else Inf)
    D1 <- (ev(c, h) - ev(c, -h)) / (2 * h)
    D2 <- (ev(c, h / 2) - ev(c, -h / 2)) / h
    R <- (4 * D2 - D1) / 3
    for (nm in names(R)) {
      s <- strsplit(nm, "_", fixed = TRUE)[[1]]
      tu <- c(s[-(1:2)], P[c])
      kk <- paste(c(s[1:2], P[sort(match(tu, P))]), collapse = "_")
      expect_true(kk %in% names(up))
      err <- max(err, abs(R[[nm]] - up[[kk]]))
    }
  }
  expect_lt(err / max(abs(up)), 1e-10)
})

test_that("skew t d4E component in (alpha, nu, nu, nu, nu) at nu = 0.05", {
  d <- skewt_distrib()
  th <- list(mu = 0, sigma = 1, alpha = -0.9, nu = 0.05)
  got <- distrib_d4expected_hessian(d, 1, th, scale = "parameter")
  ref <- 0.95092178628762332
  key <- grep("^alpha_nu_", names(got), value = TRUE)
  key <- key[vapply(strsplit(key, "_"), function(s) all(s[-(1:2)] == "nu") && length(s) == 6L, NA)]
  expect_length(key, 1L)
  expect_lt(abs(got[[key]] - ref) / ref, 2e-8)
})

ref5 <- list(
  list(mu = 1e5, theta = 0.01,
       ref = c(mu_mu = -9.9009901039113811e-06, mu_mu_theta = 9.8029595234792045e-06, mu_mu_theta_theta = -1.9411897134598582e-05, mu_mu_theta_theta_theta = 5.765938584220771e-05, mu_mu_theta_theta_theta_theta = -0.00022835513381538959, mu_mu_mu = 9.9009901088128601e-11, mu_mu_mu_theta = -9.8029585528893438e-11, mu_mu_mu_theta_theta = 1.9411991310616203e-10, mu_mu_mu_theta_theta_theta = -5.7659951018066746e-10, mu_mu_mu_mu = -1.9801980227428679e-15, mu_mu_mu_mu_theta = 1.9605915164599251e-15, mu_mu_mu_mu_theta_theta = -3.8824170973211516e-15, mu_mu_mu_mu_mu = 5.9405940711694907e-20, mu_mu_mu_mu_mu_theta = -5.8817739670260314e-20, mu_mu_mu_mu_mu_mu = -2.3762376296441508e-24, theta_theta = -0.49014797633572277, theta_theta_theta = 0.97059482485790416, theta_theta_theta_theta = -2.8829723370697673, theta_theta_theta_theta_theta = 11.417784387636575, theta_theta_theta_theta_theta_theta = -56.524048671305927, theta_theta_mu = -4.8367732954121579e-13, theta_theta_mu_theta = -4.6769293462030917e-11, theta_theta_mu_theta_theta = 3.1303704359618388e-10, theta_theta_mu_theta_theta_theta = -1.9613795115631325e-09, theta_theta_mu_mu = 9.6735456490638457e-18, theta_theta_mu_mu_theta = 9.3538568654025426e-16, theta_theta_mu_mu_theta_theta = -6.2607574801292912e-15, theta_theta_mu_mu_mu = -2.9020634121910393e-22, theta_theta_mu_mu_mu_theta = -2.8061565115197466e-20, theta_theta_mu_mu_mu_mu = 1.1608252518651807e-26, mu_theta = 4.9014797633572277e-08, mu_theta_theta = 4.8044202808714376e-06, mu_theta_theta_theta = -1.9123599263451107e-05, mu_theta_theta_theta_theta = 8.5347391673329371e-05, mu_theta_theta_theta_theta_theta = -0.00045105897063833241, mu_theta_mu = -4.9014792796798987e-13, mu_theta_mu_theta = -4.8044193295011729e-11, mu_theta_mu_theta_theta = 1.9123689671667595e-10, mu_theta_mu_theta_theta_theta = -8.5348311170665034e-10, mu_theta_mu_mu = 9.8029575920052319e-18, mu_theta_mu_mu_theta = 9.608836756262095e-16, mu_theta_mu_mu_theta_theta = -3.8247560159715013e-15, mu_theta_mu_mu_mu = -2.940886987395228e-22, mu_theta_mu_mu_mu_theta = -2.882650456056636e-20, mu_theta_mu_mu_mu_mu = 1.1763546788755662e-26)),
  list(mu = 1e6, theta = 0.001,
       ref = c(mu_mu = -9.9900099900149799e-07, mu_mu_theta = 9.9800299500798792e-07, mu_mu_theta_theta = -1.9940129740479151e-06, mu_mu_theta_theta_theta = 5.9760658443323254e-06, mu_mu_theta_theta_theta_theta = -2.3880394922614049e-05, mu_mu_mu = 9.9900099900199695e-13, mu_mu_mu_theta = -9.9800299401098205e-13, mu_mu_mu_theta_theta = 1.9940139680658691e-12, mu_mu_mu_theta_theta_theta = -5.9760718084579466e-12, mu_mu_mu_mu = -1.9980019980049918e-18, mu_mu_mu_mu_theta = 1.996005986027952e-18, mu_mu_mu_mu_theta_theta = -3.9880299241676408e-18, mu_mu_mu_mu_mu = 5.9940059940179697e-24, mu_mu_mu_mu_mu_theta = -5.9880179521018206e-24, mu_mu_mu_mu_mu_mu = -2.3976023976083855e-29, theta_theta = -0.49900149750416017, theta_theta_theta = 0.99700648669212077, theta_theta_theta_theta = -2.9880332525109647, theta_theta_theta_theta_theta = 11.940200437404719, theta_theta_theta_theta_theta_theta = -59.641399552423856, theta_theta_mu = -4.9833682634846908e-16, theta_theta_mu_theta = -4.9667714073035034e-13, theta_theta_mu_theta_theta = 3.3124071146044597e-12, theta_theta_mu_theta_theta_theta = -2.0856584014067022e-11, theta_theta_mu_mu = 9.9667365170292007e-22, theta_theta_mu_mu_theta = 9.9335427947862889e-19, theta_theta_mu_mu_theta_theta = -6.6248161934241043e-18, theta_theta_mu_mu_mu = -2.9900209521267067e-27, theta_theta_mu_mu_mu_theta = -2.9800628324896713e-24, theta_theta_mu_mu_mu_mu = 1.1960083796578611e-32, mu_theta = 4.9900149750416013e-10, mu_theta_theta = 4.9800449101746807e-07, mu_theta_theta_theta = -1.9910249401317305e-06, mu_theta_theta_theta_theta = 8.9521595570954887e-06, mu_theta_theta_theta_theta_theta = -4.7701160350066453e-05, mu_theta_mu = -4.9900149700582336e-16, mu_theta_mu_theta = -4.9800449002245412e-13, mu_theta_mu_theta_theta = 1.9910259301736047e-12, mu_theta_mu_theta_theta_theta = -8.9521694734602501e-12, mu_theta_mu_mu = 9.9800299301497303e-22, mu_theta_mu_mu_theta = 9.9600897805488029e-19, mu_theta_mu_mu_theta_theta = -3.9820538404309524e-18, mu_theta_mu_mu_mu = -2.9940089760548981e-27, mu_theta_mu_mu_mu_theta = -2.9880269281945569e-24, mu_theta_mu_mu_mu_mu = 1.197603589225951e-32)),
  list(mu = 1e6, theta = 0.01,
       ref = c(mu_mu = -9.9009900995000501e-07, mu_mu_theta = 9.802960397010196e-07, mu_mu_theta_theta = -1.9411812376159985e-06, mu_mu_theta_theta_theta = 5.7658877186057021e-06, mu_mu_theta_theta_theta_theta = -2.2835210190135592e-05, mu_mu_mu = 9.9009900999901965e-13, mu_mu_mu_theta = -9.8029602999511845e-13, mu_mu_mu_theta_theta = 1.9411821793766797e-12, mu_mu_mu_theta_theta_theta = -5.7658933703171271e-12, mu_mu_mu_mu = -1.9801980200960689e-18, mu_mu_mu_mu_theta = 1.9605920405784347e-18, mu_mu_mu_mu_theta_theta = -3.8823662422746668e-18, mu_mu_mu_mu_mu = 5.9405940605822957e-24, mu_mu_mu_mu_mu_theta = -5.8817760634998985e-24, mu_mu_mu_mu_mu_mu = -2.376237624350554e-29, theta_theta = -0.49014801986668627, theta_theta_theta = 0.97059061562075266, theta_theta_theta_theta = -2.8829441638031068, theta_theta_theta_theta_theta = 11.417607864772943, theta_theta_theta_theta_theta_theta = -56.522847569118056, theta_theta_mu = -4.8367741429966893e-15, theta_theta_mu_theta = -4.6769309905066825e-13, theta_theta_mu_theta_theta = 3.1303554885710343e-12, theta_theta_mu_theta_theta_theta = -1.9613507913779604e-11, theta_theta_mu_mu = 9.6735481918173071e-21, theta_theta_mu_mu_theta = 9.35386179831293e-19, theta_theta_mu_mu_theta_theta = -6.2607126379639086e-18, theta_theta_mu_mu_mu = -2.9020644292923711e-26, theta_theta_mu_mu_mu_theta = -2.80615848468375e-24, theta_theta_mu_mu_mu_mu = 1.1608257604158201e-31, mu_theta = 4.9014801986668623e-09, mu_theta_theta = 4.8044211371047871e-07, mu_theta_theta_theta = -1.9123517896034742e-06, mu_theta_theta_theta_theta = 8.5346564127615924e-06, mu_theta_theta_theta_theta_theta = -4.5105202983400592e-05, mu_theta_mu = -4.9014801502991212e-15, mu_theta_mu_theta = -4.8044210419677356e-13, mu_theta_mu_theta_theta = 1.9123526936861174e-12, mu_theta_mu_theta_theta_theta = -8.5346656076929787e-12, mu_theta_mu_mu = 9.8029602038627609e-21, mu_theta_mu_mu_theta = 9.6088418936613728e-19, mu_theta_mu_mu_theta_theta = -3.8247071955374684e-18, mu_theta_mu_mu_mu = -2.940888032138184e-26, mu_theta_mu_mu_mu_theta = -2.8826525110161826e-24, mu_theta_mu_mu_mu_mu = 1.1763552012470158e-31)),
  list(mu = 1e7, theta = 0.1,
       ref = c(mu_mu = -9.0909090950413218e-08, mu_mu_theta = 8.2644627347858756e-08, mu_mu_theta_theta = -1.5026296564442304e-07, mu_mu_theta_theta_theta = 4.0980810674879438e-07, mu_mu_theta_theta_theta_theta = -1.4902113650052266e-06, mu_mu_mu = 9.090909099173554e-15, mu_mu_mu_theta = -8.2644626596543974e-15, mu_mu_mu_theta_theta = 1.5026297110853036e-14, mu_mu_mu_theta_theta_theta = -4.0980814027854689e-14, mu_mu_mu_mu = -1.818181820661157e-21, mu_mu_mu_mu_theta = 1.6528925169045843e-21, mu_mu_mu_mu_theta_theta = -3.005259531452751e-21, mu_mu_mu_mu_mu = 5.4545454644628101e-28, mu_mu_mu_mu_mu_theta = -4.9586775056348667e-28, mu_mu_mu_mu_mu_mu = -2.1818181867768595e-34, theta_theta = -0.41322313686451295, theta_theta_theta = 0.7513148260592396, theta_theta_theta_theta = -2.0490405463693691, theta_theta_theta_theta_theta = 7.4510569610643316, theta_theta_theta_theta_theta_theta = -33.868442722651423, theta_theta_mu = -3.6313548163832166e-16, theta_theta_mu_theta = -2.5157661515089428e-15, theta_theta_mu_theta_theta = 1.8027415687445532e-14, theta_theta_mu_theta_theta_theta = -1.0843545187138705e-13, theta_theta_mu_mu = 7.2627095781253566e-23, theta_theta_mu_mu_theta = 5.0315322272654845e-22, theta_theta_mu_mu_theta_theta = -3.6054831316185799e-21, theta_theta_mu_mu_mu = -2.1788128570452846e-29, theta_theta_mu_mu_mu_theta = -1.5094596454539252e-28, theta_theta_mu_mu_mu_mu = 8.7152513626118474e-36, mu_theta = 4.1322313686451295e-09, mu_theta_theta = 3.3809165425858895e-08, mu_theta_theta_theta = -1.2977255974815423e-07, mu_theta_theta_theta_theta = 5.402015943001675e-07, mu_theta_theta_theta_theta_theta = -2.6417383571992181e-06, mu_theta_mu = -4.1322313323315812e-16, mu_theta_mu_theta = -3.3809164811146798e-15, mu_theta_mu_theta_theta = 1.2977256297694498e-14, mu_theta_mu_theta_theta_theta = -5.4020163753886936e-14, mu_theta_mu_mu = 8.2644625920360661e-23, mu_theta_mu_mu_theta = 6.7618328392869418e-22, mu_theta_mu_mu_theta_theta = -2.5954513241147126e-21, mu_theta_mu_mu_mu = -2.4793387558226914e-29, mu_theta_mu_mu_mu_theta = -2.0285498149033573e-28, mu_theta_mu_mu_mu_mu = 9.9173549361382521e-36))
)

test_that("negbin1 expected information and derivatives at large mu, small theta", {
  d <- negbin1_distrib()
  gen <- list(distrib_expected_hessian, distrib_dexpected_hessian, distrib_d2expected_hessian,
              distrib_d3expected_hessian, distrib_d4expected_hessian)
  for (cs in ref5) {
    th <- list(mu = cs$mu, theta = cs$theta)
    ks <- lengths(strsplit(names(cs$ref), "_")) - 2L
    grp <- vapply(strsplit(names(cs$ref), "_"), function(s) paste(s[1:2], collapse = "_"), "")
    for (k in 0:4) {
      got <- unlist(gen[[k + 1]](d, 1, th, scale = "parameter"))
      for (g in unique(grp)) {
        ref <- cs$ref[ks == k & grp == g]
        expect_true(all(names(ref) %in% names(got)))
        expect_lt(rel_err(got, ref), 1e-12)
      }
    }
  }
})

test_that("expectation() of a truncated skew t keeps the adaptive quadrature", {
  # the compiled rule serves the derivatives of the expected information only:
  # with an indicator it read 0.41923 for 0.41469
  d <- truncated(skewt_distrib(), lower = -0.5)
  th <- list(mu = 0.5, sigma = 1.2, alpha = 0.7, nu = 3)
  got <- expectation(d, function(y, theta, ...) as.numeric(y < 1), th)
  expect_equal(got, 0.414685791513054, tolerance = 1e-8)
})
