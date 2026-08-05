#' S7 Class for Probability Distributions
#'
#' @description
#' The base S7 class for probability distributions.
#' It carries the name, the parameters with their domains and links, and the support.
#'
#' @param distrib_name A single character string specifying the name of the distribution (e.g., \code{"student t"}).
#' @param dimension A character string indicating the dimensionality (\code{"univariate"} or \code{"multivariate"}).
#' @param bounds A numeric vector of length 2 defining the overall support of the distribution \code{c(lower, upper)}.
#' @param params A character vector containing the names of the distribution parameters (e.g., \code{c("mu", "sigma")}).
#' @param params_interpretation A character vector (typically named) providing the statistical interpretation of each parameter (e.g., \code{c(mu = "location")}).
#' @param n_params A numeric value specifying the total number of parameters.
#' @param params_bounds A list of numeric vectors of length 2, specifying the valid mathematical domain for each individual parameter.
#' @param link_params A list of link function objects corresponding to each parameter, primarily used to map parameters to the unconstrained real line for optimization algorithms.
#' @param params_smooth An optional named logical vector flagging, for each parameter, whether the log-likelihood is differentiable with respect to it. Defaults to all \code{TRUE} (leave empty). Set an entry to \code{FALSE} for parameters at which the log-likelihood has a kink (e.g. the location of a Laplace distribution): the observed Hessian is then degenerate and the expected information must be obtained from the score variance rather than from \eqn{-\mathbb{E}[H]} (see \code{\link{distrib_expected_hessian}}).
#'
#' @import S7
#' @section Methods:
#' Registered on the base class, so every distribution inherits them unless it
#' registers something more specific. Those that compute derivatives do so by
#' finite differences, which is why a subclass that implements nothing but
#' \code{\link{distrib_pdf}} is already fully functional. Note that
#' \code{distrib_pdf} itself has no default: the density is the one thing a
#' distribution must supply.
#'
#'   \code{\link[=distrib_deriv3.distrib]{distrib_deriv3()}},
#'   \code{\link[=distrib_deriv4.distrib]{distrib_deriv4()}},
#'   \code{\link[=distrib_expected_hessian.distrib]{distrib_expected_hessian()}},
#'   \code{\link[=distrib_gradient.distrib]{distrib_gradient()}},
#'   \code{\link[=distrib_hessian.distrib]{distrib_hessian()}},
#'   \code{\link[=generate_random_theta.distrib]{generate_random_theta()}},
#'   \code{\link[=kurtosis]{kurtosis()}},
#'   \code{\link[=mean.distrib]{mean()}},
#'   \code{\link[=print.distrib]{print()}},
#'   \code{\link[=skewness]{skewness()}},
#'   \code{\link[=std_dev]{std_dev()}},
#'   \code{\link[=variance]{variance()}}
#'
#' @return An object of class \code{distrib}.
#'
#' @examples
#' d <- gaussian1_distrib()
#' S7::S7_inherits(d, distrib)
#' d@params
#' d@params_bounds
#'
#' @export
distrib <- S7::new_class("distrib",
  properties = list(
    distrib_name          = S7::class_character,
    dimension             = S7::class_character,
    bounds                = S7::class_numeric,
    params                = S7::class_character,
    params_interpretation = S7::class_character,
    n_params              = S7::class_numeric,
    params_bounds         = S7::class_list,
    link_params           = S7::class_list,
    params_smooth         = S7::class_logical
  ),
  validator = function(self) {
    errors <- character()

    if (!self@dimension %in% c("univariate", "multivariate")) {
      errors <- c(errors, "Property 'dimension' must be 'univariate' or 'multivariate'.")
    }
    if (length(self@bounds) != 2) {
      errors <- c(errors, "Property 'bounds' must be a numeric vector of length 2 [lower, upper].")
    }
    if (self@bounds[1] >= self@bounds[2]) {
      errors <- c(errors, "Lower bound must be strictly less than the upper bound.")
    }
    if (self@n_params != length(self@params)) {
      errors <- c(errors, "Property 'n_params' must match the length of 'params'.")
    }

    missing_bounds <- setdiff(self@params, names(self@params_bounds))
    if (length(missing_bounds) > 0) {
      errors <- c(errors, paste0(
        "Property 'params_bounds' is missing an entry for: ",
        paste(missing_bounds, collapse = ", "), "."
      ))
    }

    missing_links <- setdiff(self@params, names(self@link_params))
    if (length(missing_links) > 0) {
      errors <- c(errors, paste0(
        "Property 'link_params' is missing an entry for: ",
        paste(missing_links, collapse = ", "), "."
      ))
    }

    bad_links <- names(which(!vapply(
      self@link_params,
      function(l) S7::S7_inherits(l, linkfunctions7::link),
      logical(1)
    )))
    if (length(bad_links) > 0) {
      errors <- c(errors, paste0(
        "Property 'link_params' must contain 'linkfunctions7::link' objects. Not a link: ",
        paste(bad_links, collapse = ", "), "."
      ))
    }

    if (length(self@params_smooth) > 0) {
      missing_smooth <- setdiff(self@params, names(self@params_smooth))
      if (length(missing_smooth) > 0) {
        errors <- c(errors, paste0(
          "Property 'params_smooth', when supplied, must be named and cover every parameter. Missing: ",
          paste(missing_smooth, collapse = ", "), "."
        ))
      }
    }

    if (length(errors) > 0) errors else NULL
  }
)

#' @title Default Atoms: None
#' @name distrib_atoms.distrib
#' @description
#' Ordinary distributions have nothing to declare here: a continuous one has no
#' atoms and a discrete one is made of nothing else. Only a mixed distribution
#' registers something more specific.
#' @param distrib A \code{distrib} object.
#' @param theta A named list of parameters.
#' @return A list with empty \code{y} and \code{p}.
#' @keywords internal
S7::method(distrib_atoms, distrib) <- function(distrib, theta) {
  list(y = numeric(0), p = numeric(0))
}

#' Per-Parameter Smoothness of the Log-Likelihood
#'
#' @description
#' Returns a named logical vector indicating, for each parameter of a distribution,
#' whether the log-likelihood is differentiable with respect to it. This reads the
#' \code{params_smooth} property, defaulting to all \code{TRUE} when the property was
#' left empty.
#'
#' @param distrib An object inheriting from class \code{"distrib"}.
#' @return A named logical vector, one entry per parameter.
#' @examples
#' param_smoothness(gaussian1_distrib())
#'
#' # the Laplace location is a kink, so it is not smooth
#' param_smoothness(laplace_distrib())
#'
#' @export
param_smoothness <- function(distrib) {
  ps <- distrib@params_smooth
  if (length(ps) == 0) {
    stats::setNames(rep(TRUE, distrib@n_params), distrib@params)
  } else {
    ps[distrib@params]
  }
}

#' S7 Class for Continuous Distributions
#'
#' @description A subclass of \code{distrib} specifically for continuous probability distributions.
#' @inheritParams distrib
#'
#' @section Methods:
#' Defaults for continuous distributions, built from the density alone: the cdf by
#' quadrature, the quantile function by root finding, and the generator by
#' Generalized Ratio-of-Uniforms (\code{\link{rng_grou}}) or inverse transform when
#' an analytical quantile is available.
#'
#'   \code{\link[=distrib_cdf.continuous_distrib]{distrib_cdf()}},
#'   \code{\link[=distrib_grad_y.continuous_distrib]{distrib_grad_y()}},
#'   \code{\link[=distrib_hess_y.continuous_distrib]{distrib_hess_y()}},
#'   \code{\link[=distrib_quantile.continuous_distrib]{distrib_quantile()}},
#'   \code{\link[=distrib_rng.continuous_distrib]{distrib_rng()}},
#'   \code{\link[=expectation]{expectation()}},
#'   \code{\link[=plot.continuous_distrib]{plot()}}
#'
#' Everything else is inherited from \code{\link{distrib}}.
#'
#' @return An object of class \code{continuous_distrib}.
#'
#' @examples
#' S7::S7_inherits(gaussian1_distrib(), continuous_distrib)
#' S7::S7_inherits(poisson_distrib(), continuous_distrib)
#'
#' @export
continuous_distrib <- S7::new_class("continuous_distrib",
  parent = distrib
)

#' S7 Class for Discrete Distributions
#'
#' @description A subclass of \code{distrib} specifically for discrete probability distributions.
#' @inheritParams distrib
#'
#' @section Methods:
#' Defaults for discrete distributions, built from the probability mass function
#' alone by accumulating it over the support; they require a finite lower bound,
#' which every standard count distribution has.
#'
#'   \code{\link[=distrib_cdf.discrete_distrib]{distrib_cdf()}},
#'   \code{\link[=distrib_quantile.discrete_distrib]{distrib_quantile()}},
#'   \code{\link[=distrib_rng.discrete_distrib]{distrib_rng()}},
#'   \code{\link[=expectation]{expectation()}},
#'   \code{\link[=plot.discrete_distrib]{plot()}}
#'
#' Everything else is inherited from \code{\link{distrib}}.
#'
#' @return An object of class \code{discrete_distrib}.
#'
#' @examples
#' S7::S7_inherits(poisson_distrib(), discrete_distrib)
#' S7::S7_inherits(gaussian1_distrib(), discrete_distrib)
#'
#' @export
discrete_distrib <- S7::new_class("discrete_distrib",
  parent = distrib
)
