#' @include reparametrize.R
NULL

# The second parametrizations obtained through reparametrize(), each
# supplying its hand-written map derivatives (reparam_maps.R), rather than
# written out. Each is a map of a few lines; everything else -- the density,
# the derivatives to fourth order observed and expected, the moments, the
# validator, the fit -- comes from the parent through the partition sum.
#
# Where the literature already numbers a parametrization that number is used,
# so the Weibull in the mean is weibull3 and not weibull2: in gamlss WEI2 is a
# different parametrization, and a number that means one thing there and
# another here would mislead exactly the reader who knows the field.


#' Lognormal Distribution in the Mean and Variance of Y, Obtained
#'
#' @description
#' The same family as [lognormal2_distrib()], obtained through
#' [reparametrize()] on [lognormal1_distrib()] rather than written out.
#'
#' @details
#' This exists as a reference for the tests and is not exported:
#' [lognormal2_distrib()] carries its own kernels, so the two are independent
#' implementations of one law.
#'
#' @param link_mean Link function for the mean. Defaults to the log.
#' @param link_var Link function for the variance. Defaults to the log.
#'
#' @return A reparametrized distribution object.
#'
#' @seealso [lognormal2_distrib()]
#'
#' @keywords internal
lognormal2_by_reparam <- function(link_mean = log_link(), link_var = log_link()) {
  reparametrize(
    lognormal1_distrib(),
    map = function(psi) {
      s2 <- log(1 + psi$var / psi$mean^2)
      list(mu = log(psi$mean) - s2 / 2, sigma2 = s2)
    },
    params = c("mean", "var"),
    bounds = list(mean = c(0, Inf), var = c(0, Inf)),
    links = list(mean = link_mean, var = link_var),
    map_derivs = md_lognormal2,
    interpretation = c(mean = "mean", var = "variance"),
    name = "lognormal2"
  )
}


#' Weibull Distribution in the Mean, Obtained
#'
#' @description
#' The same family as [weibull3_distrib()], obtained through
#' [reparametrize()] on [weibull1_distrib()] rather than written out.
#'
#' @details
#' This exists as a reference for the tests and is not exported:
#' [weibull3_distrib()] carries its own kernels, so the two are independent
#' implementations of one law.
#'
#' @param link_mean Link function for the mean. Defaults to the log.
#' @param link_sigma Link function for the shape. Defaults to the log.
#'
#' @return A reparametrized distribution object.
#'
#' @seealso [weibull3_distrib()]
#'
#' @keywords internal
weibull3_by_reparam <- function(link_mean = log_link(), link_sigma = log_link()) {
  reparametrize(
    weibull1_distrib(),
    map = function(psi) {
      list(mu = psi$mean / gamma(1 + 1 / psi$sigma), sigma = psi$sigma)
    },
    params = c("mean", "sigma"),
    bounds = list(mean = c(0, Inf), sigma = c(0, Inf)),
    links = list(mean = link_mean, sigma = link_sigma),
    map_derivs = md_weibull3,
    interpretation = c(mean = "mean", sigma = "shape"),
    name = "weibull3"
  )
}


#' Student t Distribution in the Standard Deviation, Obtained
#'
#' @description
#' The same family as [student_t2_distrib()], obtained through
#' [reparametrize()] on [student_t1_distrib()] rather than written out.
#'
#' @details
#' This exists as a reference for the tests and is not exported:
#' [student_t2_distrib()] carries its own kernels, so the two are independent
#' implementations of one law.
#'
#' @param link_mu Link function for the location. Defaults to the identity.
#' @param link_sigma Link function for the standard deviation. Defaults to the
#'   log.
#' @param link_nu Link function for the degrees of freedom. Defaults to a link
#'   bounded below at two.
#'
#' @return A reparametrized distribution object.
#'
#' @seealso [student_t2_distrib()]
#'
#' @keywords internal
student_t2_by_reparam <- function(link_mu = identity_link(),
                                  link_sigma = log_link(),
                                  link_nu = bounded_link(lwr = 2)) {
  reparametrize(
    student_t1_distrib(),
    map = function(psi) {
      list(mu = psi$mu,
           sigma = psi$sigma * sqrt((psi$nu - 2) / psi$nu),
           nu = psi$nu)
    },
    params = c("mu", "sigma", "nu"),
    bounds = list(mu = c(-Inf, Inf), sigma = c(0, Inf), nu = c(2, Inf)),
    links = list(mu = link_mu, sigma = link_sigma, nu = link_nu),
    map_derivs = md_student_t2,
    interpretation = c(mu = "location", sigma = "standard deviation",
                       nu = "degrees of freedom"),
    name = "student t2"
  )
}


#' Generalized Gamma Distribution in the Mean, Obtained
#'
#' @description
#' The same family as [gengamma2_distrib()], obtained through
#' [reparametrize()] on [gengamma1_distrib()] rather than written out.
#'
#' @details
#' This exists as a reference for the tests and is not exported:
#' [gengamma2_distrib()] carries its own kernels, so the two are independent
#' implementations of one law.
#'
#' @param link_mean Link function for the mean. Defaults to the log.
#' @param link_d Link function for the shape. Defaults to the log.
#' @param link_p Link function for the power. Defaults to the log.
#'
#' @return A reparametrized distribution object.
#'
#' @seealso [gengamma2_distrib()]
#'
#' @keywords internal
gengamma2_by_reparam <- function(link_mean = log_link(), link_d = log_link(),
                                 link_p = log_link()) {
  reparametrize(
    gengamma1_distrib(),
    map = function(psi) {
      list(a = psi$mean * gamma(psi$d / psi$p) / gamma((psi$d + 1) / psi$p),
           d = psi$d, p = psi$p)
    },
    params = c("mean", "d", "p"),
    bounds = list(mean = c(0, Inf), d = c(0, Inf), p = c(0, Inf)),
    links = list(mean = link_mean, d = link_d, p = link_p),
    map_derivs = md_gengamma2,
    interpretation = c(mean = "mean", d = "shape", p = "power"),
    name = "gengamma2"
  )
}


#' Inverse Gaussian Distribution in the Mean and Shape, Obtained
#'
#' @description
#' The same family as [invgauss2_distrib()], obtained through
#' [reparametrize()] rather than written out.
#'
#' @details
#' This exists as a check rather than as a second way of doing the same thing.
#' [invgauss2_distrib()] carries its own kernels, so the two are
#' independent implementations of one object and their agreement needs no
#' tolerance to be chosen. It is not exported for that reason.
#'
#' @return A reparametrized distribution object.
#'
#' @seealso [invgauss2_distrib()]
#'
#' @keywords internal
invgauss2_by_reparam <- function() {
  reparametrize(
    invgauss1_distrib(),
    map = function(psi) list(mu = psi$mu, phi = 1 / psi$lambda),
    params = c("mu", "lambda"),
    bounds = list(mu = c(0, Inf), lambda = c(0, Inf)),
    links = list(mu = log_link(), lambda = log_link()),
    map_derivs = md_invgauss2,
    interpretation = c(mu = "mean", lambda = "shape"),
    name = "invgauss2 (reparametrized)"
  )
}
