# The skew t's derivatives in the degrees of freedom are exact: the
# derivatives of the t distribution function in its degrees of freedom are
# integrals, taken by quadrature (pt_skewt_exact.h). Values from mpmath:
# the ratios (d_m^j T_m(w))/T_m(w) by mpmath.diff on the regularized incomplete
# beta at 150 digits, and the derivatives of the log-density by mpmath.diff at
# 50 digits.

test_that("the derivatives of the t distribution function in its degrees of freedom", {
  ref <- rbind(
    c(-40, 1.3, -3.31953115260994011137, -42.1218816277547929276, -664.814153025691129680),
    c(-3, 1.3, -0.828573546810343331933, -2.06125089857357092402, -13.9158502614536456303),
    c(-0.5, 1.3, -0.0646803651580502829729, -0.140950228444827371605, -0.999509738633082261204),
    c(0.3, 1.3, 0.0214592663187527485563, 0.0472263974581798125754, 0.338545011507710437034),
    c(-40, 5, -2.47717827407196558441, -16.0852866294510464983, -113.090045293867866759),
    c(-3, 5, -0.251119473907996779625, -0.0780209204057726010582, -0.0651537038601256786366),
    c(-0.5, 5, -0.00637932831222599899309, -0.00138758149518837591254, -0.000988322343097450808514),
    c(0.3, 5, 0.00189851164944845955603, 0.000415087025738214689510, 0.000296848103536429768462))
  for (i in seq_len(nrow(ref))) {
    got <- skewt_tdf_ratio_cpp(ref[i, 1], ref[i, 2])
    expect_equal(got[1, c(1, 3, 5)], ref[i, 3:5], tolerance = 1e-13,
                 label = sprintf("w = %g, m = %g", ref[i, 1], ref[i, 2]))
  }
})

test_that("the components in nu agree with mpmath at every order", {
  d <- skewt_distrib()
  th <- list(mu = 0.2, sigma = 1.3, alpha = 0.7, nu = 4)
  g <- c(distrib_gradient(d, 0.7, th), distrib_hessian(d, 0.7, th),
         distrib_deriv3(d, 0.7, th), distrib_deriv4(d, 0.7, th),
         distrib_deriv5(d, 0.7, th))
  ref <- c(nu = 0.0176062963276277110798, mu_nu = -0.0158839182194376117145,
           alpha_nu = -0.00131074010163531792124, nu_nu = -0.00857895090531822452542,
           sigma_nu_nu = 0.00295064129548065838924, nu_nu_nu = 0.00624075687400320957757,
           mu_mu_nu_nu = -0.0156172627207254317112, nu_nu_nu_nu = -0.00602957566654267134509,
           alpha_nu_nu_nu_nu = 0.000449069703790849527796,
           nu_nu_nu_nu_nu = 0.00725948644291969396617)
  for (k in names(ref)) expect_equal(g[[k]], ref[[k]], tolerance = 1e-12, label = k)
  # the stencils they replace read the score in nu 1e-13 to 1.6e-10 out and
  # nu_nu_nu_nu 2e-4 out; far in the tail of z
  g3 <- c(distrib_gradient(d, -40, th), distrib_hessian(d, -40, th),
          distrib_deriv3(d, -40, th), distrib_deriv4(d, -40, th))
  expect_equal(g3[["nu"]], -2.36750662410320496777, tolerance = 1e-13)
  expect_equal(g3[["nu_nu"]], 0.0960561785096225106370, tolerance = 1e-13)
  expect_equal(g3[["sigma_alpha_nu"]], 0.00433685387297816966467, tolerance = 1e-12)
  expect_equal(g3[["nu_nu_nu_nu"]], 0.000538678320866201302611, tolerance = 1e-11)
})

test_that("the components hold in both tails of the tilting argument", {
  d <- skewt_distrib()
  # w = alpha z sqrt((nu+1)/(nu+z^2)) near -78, where q = t/T is nearly
  # d_w log t and (d_m T)/T nearly (d_m t)/t; expanded in them the derivatives
  # in nu lost up to 1e5
  th <- list(mu = 0, sigma = 1, alpha = 40, nu = 60)
  expect_equal(distrib_hessian(d, -2, th)$nu_nu, 0.00817031428662349746107,
               tolerance = 1e-12)
  expect_equal(distrib_deriv3(d, -2, th)$nu_nu_nu, -0.000138405202943803224326,
               tolerance = 1e-10)
  expect_equal(distrib_deriv4(d, -2, th)$alpha_alpha_alpha_alpha,
               0.000129379432448811628152, tolerance = 1e-8)
  # and near w = 9.3, where T is one to within 1e-7
  th2 <- list(mu = 0, sigma = 1, alpha = 5, nu = 4)
  h <- distrib_hessian(d, 3, th2)
  expect_equal(h$nu_nu, 0.0239935912984086975367, tolerance = 1e-13)
  expect_equal(h$alpha_nu, -0.000134025699859189170172, tolerance = 1e-13)
})
