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


#' Lognormal Distribution in the Mean and Variance of Y
#'
#' @description
#' Creates a lognormal distribution object parametrized by the mean and the
#' variance of \eqn{Y} itself, rather than of \eqn{\log Y}.
#'
#' @details
#' The parameters of \code{\link{lognormal1_distrib}} describe \eqn{\log Y}, so
#' neither of them is a moment of \eqn{Y}. Here they are, through
#' \deqn{\mu_{\log} = \log\dfrac{m^2}{\sqrt{v + m^2}}, \qquad
#'       \sigma^2_{\log} = \log\left(1 + \dfrac{v}{m^2}\right)}
#' which is the parametrization a regression on the mean wants.
#'
#' Built with \code{\link{reparametrize}}, so every derivative to fourth order,
#' observed and expected, is exact.
#'
#' @section The distribution:
#' \deqn{f(y) = \frac{1}{y\sqrt{2\pi s^{2}}}\exp\!\left\{-\frac{(\log y - m)^{2}}{2s^{2}}\right\}, \quad s^{2} = \log\!\left(1+\frac{v}{\mu^{2}}\right)\!, \; m = \log\mu - \frac{s^{2}}{2}}
#' on \eqn{y \in (0, \infty)}.
#'
#' \deqn{\mathbb{E}[Y] = \mu, \qquad \operatorname{Var}(Y) = v}
#'
#' @param link_mean Link function for the mean. Defaults to the log.
#' @param link_var Link function for the variance. Defaults to the log.
#'
#' @return A reparametrized distribution object.
#'
#' @seealso \code{\link{lognormal1_distrib}}, \code{\link{reparametrize}}
#'
#' @examples
#' d <- lognormal2_distrib()
#' theta <- list(mean = 3, var = 2)
#' c(mean = mean(d, theta), variance = variance(d, theta))
#'
#' @export
lognormal2_distrib <- function(link_mean = log_link(), link_var = log_link()) {
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


#' Weibull Distribution in the Mean
#'
#' @description
#' Creates a Weibull distribution object parametrized by its mean and its
#' shape.
#'
#' @details
#' The first parameter of \code{\link{weibull1_distrib}} is the scale and not
#' the mean: the mean is \eqn{\mu\,\Gamma(1 + 1/\sigma)}. Inverting that gives
#' the map used here,
#' \deqn{\mu = \dfrac{m}{\Gamma(1 + 1/\sigma)},}
#' so every derivative becomes a derivative of the gamma function, which is
#' why \code{\link{weibull1_distrib}} is not written this way.
#'
#' The number follows gamlss, where the Weibull in the mean is \code{WEI3}.
#' Leaving \code{weibull2} unused is deliberate: it names a different
#' parametrization there.
#'
#' @section The distribution:
#' \deqn{f(y) = \frac{\sigma}{b}\left(\frac{y}{b}\right)^{\sigma-1}e^{-(y/b)^{\sigma}}, \qquad b = \frac{\mu}{\Gamma(1+1/\sigma)}}
#' on \eqn{y \in (0, \infty)}.
#'
#' \deqn{\mathbb{E}[Y] = \mu, \qquad \operatorname{Var}(Y) = b^{2}\left[\Gamma(1+2/\sigma) - \Gamma(1+1/\sigma)^{2}\right]}
#'
#' @param link_mean Link function for the mean. Defaults to the log.
#' @param link_sigma Link function for the shape. Defaults to the log.
#'
#' @return A reparametrized distribution object.
#'
#' @seealso \code{\link{weibull1_distrib}}, \code{\link{reparametrize}}
#'
#' @examples
#' d <- weibull3_distrib()
#' theta <- list(mean = 4, sigma = 1.7)
#' mean(d, theta)
#'
#' @references
#' Rigby, R. A. and Stasinopoulos, D. M. (2005). Generalized additive models
#' for location, scale and shape. \emph{Journal of the Royal Statistical
#' Society, Series C} 54, 507-554.
#'
#' @export
weibull3_distrib <- function(link_mean = log_link(), link_sigma = log_link()) {
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


#' Student t Distribution in the Standard Deviation
#'
#' @description
#' Creates a Student t distribution object whose second parameter is the
#' standard deviation rather than the scale.
#'
#' @details
#' The scale of \code{\link{student_t1_distrib}} is not the standard
#' deviation: the two differ by \eqn{\sqrt{\nu/(\nu-2)}}. Here the map is
#' \deqn{\sigma_{\text{scale}} = \sigma\sqrt{\dfrac{\nu-2}{\nu}},}
#' which exists only for \eqn{\nu > 2}, and the constructor bounds \eqn{\nu}
#' there rather than letting the map return a complex number several frames
#' down. This is \code{TF2} in gamlss.
#'
#' The restriction is the point rather than a limitation: a family
#' parametrized by a standard deviation is a family whose standard deviation
#' exists.
#'
#' @section The distribution:
#' \deqn{f(y) = \frac{1}{s_0}\,t_{\nu}\!\left(\frac{y-\mu}{s_0}\right), \qquad s_0 = \sigma\sqrt{\frac{\nu-2}{\nu}}}
#' on \eqn{y \in \mathbb{R}}.
#'
#' \deqn{\mathbb{E}[Y] = \mu, \qquad \operatorname{Var}(Y) = \sigma^{2}}
#'
#' @param link_mu Link function for the location. Defaults to the identity.
#' @param link_sigma Link function for the standard deviation. Defaults to the
#'   log.
#' @param link_nu Link function for the degrees of freedom. Defaults to a link
#'   bounded below at two.
#'
#' @return A reparametrized distribution object.
#'
#' @seealso \code{\link{student_t1_distrib}}, \code{\link{reparametrize}}
#'
#' @examples
#' d <- student_t2_distrib()
#' theta <- list(mu = 0, sigma = 2, nu = 8)
#' variance(d, theta)
#'
#' @export
student_t2_distrib <- function(link_mu = identity_link(),
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


#' Generalized Gamma Distribution in the Mean
#'
#' @description
#' Creates a generalized gamma distribution object whose first parameter is the
#' mean.
#'
#' @details
#' The Stacy parametrization of \code{\link{gengamma1_distrib}} carries a
#' scale, a shape and a power, and exposes no mean at all, which is awkward for
#' a family a regression would put a linear predictor on. Since
#' \eqn{\mathbb{E}[Y] = a\,\Gamma((d+1)/p)/\Gamma(d/p)}, the map is
#' \deqn{a = m\,\dfrac{\Gamma(d/p)}{\Gamma((d+1)/p)}.}
#'
#' @section The distribution:
#' \deqn{f(y) = \frac{p\,y^{d-1}}{a^{d}\,\Gamma(d/p)}\,e^{-(y/a)^{p}}, \qquad a = \mu\,\frac{\Gamma(d/p)}{\Gamma((d+1)/p)}}
#' on \eqn{y \in (0, \infty)}.
#'
#' \deqn{\mathbb{E}[Y] = \mu}
#'
#' @param link_mean Link function for the mean. Defaults to the log.
#' @param link_d Link function for the shape. Defaults to the log.
#' @param link_p Link function for the power. Defaults to the log.
#'
#' @return A reparametrized distribution object.
#'
#' @seealso \code{\link{gengamma1_distrib}}, \code{\link{reparametrize}}
#'
#' @examples
#' d <- gengamma2_distrib()
#' theta <- list(mean = 5, d = 3, p = 1.5)
#' mean(d, theta)
#'
#' @export
gengamma2_distrib <- function(link_mean = log_link(), link_d = log_link(),
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
#' The same family as \code{\link{invgauss2_distrib}}, obtained through
#' \code{\link{reparametrize}} rather than written out.
#'
#' @details
#' This exists as a check rather than as a second way of doing the same thing.
#' \code{\link{invgauss2_distrib}} carries its own kernels, so the two are
#' independent implementations of one object and their agreement needs no
#' tolerance to be chosen. It is not exported for that reason.
#'
#' @return A reparametrized distribution object.
#'
#' @seealso \code{\link{invgauss2_distrib}}
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
