#' S7 Class for Probability Distributions
#'
#' @description
#' A base class for probability distributions using the S7 object-oriented system.
#' It encapsulates the essential metadata, parameter definitions, and boundaries of a distribution.
#'
#' @param distrib_name A single character string specifying the name of the distribution (e.g., \code{"student t"}).
#' @param dimension A character string indicating the dimensionality (\code{"univariate"} or \code{"multivariate"}).
#' @param bounds A numeric vector of length 2 defining the overall support of the distribution \code{c(lower, upper)}.
#' @param params A character vector containing the names of the distribution parameters (e.g., \code{c("mu", "sigma")}).
#' @param params_interpretation A character vector (typically named) providing the statistical interpretation of each parameter (e.g., \code{c(mu = "location")}).
#' @param n_params A numeric value specifying the total number of parameters.
#' @param params_bounds A list of numeric vectors of length 2, specifying the valid mathematical domain for each individual parameter.
#' @param link_params A list of link function objects corresponding to each parameter, primarily used to map parameters to the unconstrained real line for optimization algorithms.
#'
#' @import S7
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
    link_params           = S7::class_list
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
    
    if (length(errors) > 0) errors else NULL
  }
)

#' S7 Class for Continuous Distributions
#'
#' @description A subclass of \code{distrib} specifically for continuous probability distributions.
#' @inheritParams distrib
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
#' @export
discrete_distrib <- S7::new_class("discrete_distrib",
  parent = distrib
)
