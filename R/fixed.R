#' @include distrib.R generics.R utility_functions.R moments.R cross_derivatives.R cross2_derivatives.R cross3_derivatives.R cross_theta2_derivatives.R mv_summary.R y_higher.R
NULL

#' @title S7 Class for Distributions With Fixed Parameters (Continuous)
#' @name FixedContinuousDistrib
#'
#' @section Methods:
#' The following methods are registered on this class, and each is linked to
#' the page that documents it:
#' - the distribution: [`distrib_pdf()`][distrib_pdf.FixedContinuousDistrib],
#'   [`distrib_cdf()`][distrib_cdf.FixedContinuousDistrib],
#'   [`distrib_quantile()`][distrib_quantile.FixedContinuousDistrib],
#'   [`distrib_rng()`][distrib_rng.FixedContinuousDistrib] and
#'   [`distrib_atoms()`][distrib_atoms.FixedContinuousDistrib];
#' - the derivatives in the parameters:
#'   [`distrib_gradient()`][distrib_gradient.FixedContinuousDistrib],
#'   [`distrib_hessian()`][distrib_hessian.FixedContinuousDistrib],
#'   [`distrib_expected_hessian()`][distrib_expected_hessian.FixedContinuousDistrib],
#'   [`distrib_deriv3()`][distrib_deriv3.FixedContinuousDistrib],
#'   [`distrib_deriv4()`][distrib_deriv4.FixedContinuousDistrib] and
#'   [`distrib_deriv5()`][distrib_deriv5.FixedContinuousDistrib];
#' - the derivatives of the expected information: [distrib_dexpected_hessian()]
#'   and [distrib_d2expected_hessian()];
#' - the derivatives in the response:
#'   [`distrib_grad_y()`][distrib_grad_y.FixedContinuousDistrib],
#'   [`distrib_hess_y()`][distrib_hess_y.FixedContinuousDistrib],
#'   [`distrib_deriv3_y()`][distrib_deriv3_y.FixedContinuousDistrib] and
#'   [`distrib_deriv4_y()`][distrib_deriv4_y.FixedContinuousDistrib];
#' - the mixed derivatives in the response and the parameters:
#'   [`distrib_cross_y()`][distrib_cross_y.FixedContinuousDistrib],
#'   [`distrib_cross2_y()`][distrib_cross2_y.FixedContinuousDistrib],
#'   [`distrib_cross3_y()`][distrib_cross3_y.FixedContinuousDistrib],
#'   [`distrib_grad_y_hess()`][distrib_grad_y_hess.FixedContinuousDistrib] and
#'   [`distrib_hess_y_hess()`][distrib_hess_y_hess.FixedContinuousDistrib];
#' - the derivatives of the distribution function:
#'   [`distrib_grad_cdf()`][distrib_grad_cdf.FixedContinuousDistrib] and
#'   [`distrib_hess_cdf()`][distrib_hess_cdf.FixedContinuousDistrib];
#' - the moments: [`variance()`][variance.FixedContinuousDistrib],
#'   [`std_dev()`][std_dev.FixedContinuousDistrib],
#'   [`skewness()`][skewness.FixedContinuousDistrib] and
#'   [`kurtosis()`][kurtosis.FixedContinuousDistrib];
#' - the exactness flags: [expected_hessian_exact()] and [has_exact_deriv4()];
#' - the fitting helpers: [distrib_scalar_route()].
#'
#' Every other generic is inherited from [continuous_distrib()].
#'
#' @aliases distrib_atoms.FixedContinuousDistrib
#' @aliases distrib_cdf.FixedContinuousDistrib
#' @aliases distrib_cross2_y.FixedContinuousDistrib
#' @aliases distrib_cross3_y.FixedContinuousDistrib
#' @aliases distrib_cross_y.FixedContinuousDistrib
#' @aliases distrib_deriv3.FixedContinuousDistrib
#' @aliases distrib_deriv3_y.FixedContinuousDistrib
#' @aliases distrib_deriv4.FixedContinuousDistrib
#' @aliases distrib_deriv4_y.FixedContinuousDistrib
#' @aliases distrib_expected_hessian.FixedContinuousDistrib
#' @aliases distrib_grad_cdf.FixedContinuousDistrib
#' @aliases distrib_grad_y.FixedContinuousDistrib
#' @aliases distrib_grad_y_hess.FixedContinuousDistrib
#' @aliases distrib_gradient.FixedContinuousDistrib
#' @aliases distrib_hess_cdf.FixedContinuousDistrib
#' @aliases distrib_hess_y.FixedContinuousDistrib
#' @aliases distrib_hess_y_hess.FixedContinuousDistrib
#' @aliases distrib_hessian.FixedContinuousDistrib
#' @aliases distrib_pdf.FixedContinuousDistrib
#' @aliases distrib_quantile.FixedContinuousDistrib
#' @aliases distrib_rng.FixedContinuousDistrib
#' @aliases kurtosis.FixedContinuousDistrib
#' @aliases skewness.FixedContinuousDistrib
#' @aliases std_dev.FixedContinuousDistrib
#' @aliases variance.FixedContinuousDistrib
#'
#' @description
#' This is the S7 class of a CONTINUOUS distribution in which some parameters of
#' the wrapped distribution are held at known values. It is the only wrapper in
#' this package that REMOVES parameters: the law is the parent's, evaluated at a
#' `theta` short by however many were fixed. Build one with [fixed()], which
#' picks this class or one of its two siblings from the parent.
#'
#' @details
#' # How every method works
#'
#' The free parameters are the parent's minus the fixed ones, IN THE PARENT'S
#' ORDER. Every method splices the fixed values back into `theta` at their
#' positions and delegates to the parent, so the parent's closed forms are used
#' wherever they exist, and a derivative method then keeps only the components
#' in which every index is a free parameter. Every method of this class
#' delegates its computation to the parent, so the speed and the accuracy of
#' each method are those of the parent's method.
#'
#' # Why the components can be subset by name
#'
#' The free set preserves the parent's ORDER, so a combination of free
#' parameters produces the same name string under the parent's enumeration as
#' under the wrapper's. Subsetting by name therefore cannot pair a component
#' with the wrong indices. Re-parsing a name by splitting on the underscore
#' makes this mistake for a parameter whose own name contains an underscore.
#'
#' # The methods carry no documentation pages
#'
#' The continuous and discrete branches register their methods inside a loop
#' over the two classes, so roxygen has no top-level assignment to which it can
#' attach a block. The multivariate branch registers its methods at the top
#' level and leaves them undocumented to match, because each of them is the
#' same delegation.
#' Read the parent's page for what a method computes and this page for what the
#' wrapper does to it.
#'
#' @param parent_distrib The wrapped `continuous_distrib` object.
#' @param fixed_params A named list of the fixed values, one single finite
#'   number each, strictly inside the corresponding parameter's open domain.
#' @inheritParams distrib
#'
#' @return An S7 object of class `FixedContinuousDistrib`, inheriting from
#'   `continuous_distrib` and from `distrib`. It carries `parent_distrib` and
#'   `fixed_params` beside the parent's properties. For an object built by
#'   [fixed()], `params` is the parent's less the fixed names, `n_params` falls
#'   by as many, `bounds` is the parent's unchanged, and `distrib_name` is the
#'   parent's with the held values in brackets.
#'
#' @seealso [fixed()] to build one, [FixedDiscreteDistrib] and
#'   [FixedMultivariateDistrib] for the two siblings, and [folded()], whose
#'   half-normal is a fixed folded gaussian.
#'
#' @examples
#' # A gaussian centered at zero has one free parameter instead of two.
#' d <- fixed(gaussian1_distrib(), mu = 0)
#' d@params
#' d@fixed_params
#' d@distrib_name
#'
#' # The law is the parent's at the spliced theta.
#' all.equal(distrib_pdf(d, c(-1, 0, 1), list(sigma = 2)),
#'           dnorm(c(-1, 0, 1), 0, 2))
#'
#' # A derivative keeps only the components among the free parameters.
#' set.seed(1)
#' y <- rnorm(20, 0, 2)
#' names(distrib_gradient(d, y, list(sigma = 2)))
#' names(distrib_deriv3(d, y, list(sigma = 2)))
#'
#' # Fixing every parameter is allowed and gives a fully known law.
#' d0 <- fixed(gaussian1_distrib(), mu = 0, sigma = 1)
#' c(n_params = d0@n_params, density = distrib_pdf(d0, 0, list()))
FixedContinuousDistrib <- S7::new_class("FixedContinuousDistrib",
  parent = continuous_distrib,
  properties = list(
    parent_distrib = distrib,
    fixed_params = S7::class_list
  )
)

#' @title S7 Class for Distributions With Fixed Parameters (Discrete)
#' @name FixedDiscreteDistrib
#'
#' @section Methods:
#' The following methods are registered on this class, and each is linked to
#' the page that documents it:
#' - the distribution: [`distrib_pdf()`][distrib_pdf.FixedDiscreteDistrib],
#'   [`distrib_cdf()`][distrib_cdf.FixedDiscreteDistrib],
#'   [`distrib_quantile()`][distrib_quantile.FixedDiscreteDistrib],
#'   [`distrib_rng()`][distrib_rng.FixedDiscreteDistrib] and
#'   [`distrib_atoms()`][distrib_atoms.FixedDiscreteDistrib];
#' - the derivatives in the parameters:
#'   [`distrib_gradient()`][distrib_gradient.FixedDiscreteDistrib],
#'   [`distrib_hessian()`][distrib_hessian.FixedDiscreteDistrib],
#'   [`distrib_expected_hessian()`][distrib_expected_hessian.FixedDiscreteDistrib],
#'   [`distrib_deriv3()`][distrib_deriv3.FixedDiscreteDistrib],
#'   [`distrib_deriv4()`][distrib_deriv4.FixedDiscreteDistrib] and
#'   [`distrib_deriv5()`][distrib_deriv5.FixedDiscreteDistrib];
#' - the derivatives of the expected information: [distrib_dexpected_hessian()]
#'   and [distrib_d2expected_hessian()];
#' - the derivatives in the response:
#'   [`distrib_grad_y()`][distrib_grad_y.FixedDiscreteDistrib] and
#'   [`distrib_hess_y()`][distrib_hess_y.FixedDiscreteDistrib];
#' - the mixed derivatives in the response and the parameters:
#'   [`distrib_cross_y()`][distrib_cross_y.FixedDiscreteDistrib],
#'   [`distrib_cross2_y()`][distrib_cross2_y.FixedDiscreteDistrib],
#'   [`distrib_grad_y_hess()`][distrib_grad_y_hess.FixedDiscreteDistrib] and
#'   [`distrib_hess_y_hess()`][distrib_hess_y_hess.FixedDiscreteDistrib];
#' - the derivatives of the distribution function:
#'   [`distrib_grad_cdf()`][distrib_grad_cdf.FixedDiscreteDistrib] and
#'   [`distrib_hess_cdf()`][distrib_hess_cdf.FixedDiscreteDistrib];
#' - the moments: [`variance()`][variance.FixedDiscreteDistrib],
#'   [`std_dev()`][std_dev.FixedDiscreteDistrib],
#'   [`skewness()`][skewness.FixedDiscreteDistrib] and
#'   [`kurtosis()`][kurtosis.FixedDiscreteDistrib];
#' - the exactness flags: [expected_hessian_exact()] and [has_exact_deriv4()];
#' - the fitting helpers: [distrib_scalar_route()].
#'
#' Every other generic is inherited from [discrete_distrib()].
#'
#' @aliases distrib_atoms.FixedDiscreteDistrib
#' @aliases distrib_cdf.FixedDiscreteDistrib
#' @aliases distrib_cross2_y.FixedDiscreteDistrib
#' @aliases distrib_cross_y.FixedDiscreteDistrib
#' @aliases distrib_deriv3.FixedDiscreteDistrib
#' @aliases distrib_deriv4.FixedDiscreteDistrib
#' @aliases distrib_expected_hessian.FixedDiscreteDistrib
#' @aliases distrib_grad_cdf.FixedDiscreteDistrib
#' @aliases distrib_grad_y.FixedDiscreteDistrib
#' @aliases distrib_grad_y_hess.FixedDiscreteDistrib
#' @aliases distrib_gradient.FixedDiscreteDistrib
#' @aliases distrib_hess_cdf.FixedDiscreteDistrib
#' @aliases distrib_hess_y.FixedDiscreteDistrib
#' @aliases distrib_hess_y_hess.FixedDiscreteDistrib
#' @aliases distrib_hessian.FixedDiscreteDistrib
#' @aliases distrib_pdf.FixedDiscreteDistrib
#' @aliases distrib_quantile.FixedDiscreteDistrib
#' @aliases distrib_rng.FixedDiscreteDistrib
#' @aliases kurtosis.FixedDiscreteDistrib
#' @aliases skewness.FixedDiscreteDistrib
#' @aliases std_dev.FixedDiscreteDistrib
#' @aliases variance.FixedDiscreteDistrib
#'
#' @description
#' This is the S7 class of a DISCRETE distribution in which some parameters of
#' the wrapped distribution are held at known values. It behaves exactly as
#' [FixedContinuousDistrib] does; the split into three classes exists so that
#' the wrapper inherits the right base class, and with it the right defaults for
#' anything that the parent does not register. Build one with [fixed()], which
#' picks the class from the parent.
#'
#' @details
#' # How every method works
#'
#' The free parameters are the parent's minus the fixed ones, IN THE PARENT'S
#' ORDER. Every method splices the fixed values back into `theta` at their
#' positions and delegates to the parent, so the parent's closed forms are used
#' wherever they exist, and a derivative method then keeps only the components
#' in which every index is a free parameter. Every method of this class
#' delegates its computation to the parent, so the speed and the accuracy of
#' each method are those of the parent's method.
#'
#' # Why the components can be subset by name
#'
#' The free set preserves the parent's ORDER, so a combination of free
#' parameters produces the same name string under the parent's enumeration as
#' under the wrapper's. Subsetting by name therefore cannot pair a component
#' with the wrong indices. Re-parsing a name by splitting on the underscore
#' makes this mistake for a parameter whose own name contains an underscore.
#'
#' # The methods carry no documentation pages
#'
#' The continuous and discrete branches register their methods inside a loop
#' over the two classes, so roxygen has no top-level assignment to which it can
#' attach a block. The multivariate branch registers its methods at the top
#' level and leaves them undocumented to match, because each of them is the
#' same delegation.
#' Read the parent's page for what a method computes and this page for what the
#' wrapper does to it.
#'
#' @param parent_distrib The wrapped `discrete_distrib` object.
#' @param fixed_params A named list of the fixed values, one single finite
#'   number each, strictly inside the corresponding parameter's open domain.
#' @inheritParams distrib
#'
#' @return An S7 object of class `FixedDiscreteDistrib`, inheriting from
#'   `discrete_distrib` and from `distrib`. It carries `parent_distrib` and
#'   `fixed_params` beside the parent's properties.
#'
#' @seealso [fixed()] to build one, [FixedContinuousDistrib] and
#'   [FixedMultivariateDistrib] for the two siblings, and [zero_inflated()],
#'   one of whose own parameters can usefully be held.
#'
#' @examples
#' # A Poisson with its mean known has no free parameter at all.
#' d <- fixed(poisson_distrib(), mu = 3)
#' c(n_params = d@n_params, class = class(d)[1])
#' all.equal(distrib_pdf(d, 0:3, list()), dpois(0:3, 3))
#'
#' # A wrapper's OWN parameter can be held, so a zero-inflated model can have
#' # a known inflation rate.
#' zi <- fixed(zero_inflated(poisson_distrib()), zi = 0.3)
#' zi@params
#' all.equal(distrib_pdf(zi, 0:3, list(mu = 3)),
#'           0.3 * (0:3 == 0) + 0.7 * dpois(0:3, 3))
FixedDiscreteDistrib <- S7::new_class("FixedDiscreteDistrib",
  parent = discrete_distrib,
  properties = list(
    parent_distrib = distrib,
    fixed_params = S7::class_list
  )
)

#' @title S7 Class for Distributions With Fixed Parameters (Multivariate)
#' @name FixedMultivariateDistrib
#'
#' @section Methods:
#' The following methods are registered on this class, and each is linked to
#' the page that documents it:
#' - the distribution: [`distrib_pdf()`][distrib_pdf.FixedMultivariateDistrib]
#'   and [`distrib_rng()`][distrib_rng.FixedMultivariateDistrib];
#' - the derivatives in the parameters:
#'   [`distrib_gradient()`][distrib_gradient.FixedMultivariateDistrib],
#'   [`distrib_hessian()`][distrib_hessian.FixedMultivariateDistrib],
#'   [`distrib_expected_hessian()`][distrib_expected_hessian.FixedMultivariateDistrib],
#'   [`distrib_deriv3()`][distrib_deriv3.FixedMultivariateDistrib] and
#'   [`distrib_deriv4()`][distrib_deriv4.FixedMultivariateDistrib];
#' - the derivatives in the response:
#'   [`distrib_grad_y()`][distrib_grad_y.FixedMultivariateDistrib],
#'   [`distrib_hess_y()`][distrib_hess_y.FixedMultivariateDistrib],
#'   [distrib_deriv3_y()] and [distrib_deriv4_y()];
#' - the mixed derivatives in the response and the parameters:
#'   [`distrib_cross_y()`][distrib_cross_y.FixedMultivariateDistrib],
#'   [`distrib_cross2_y()`][distrib_cross2_y.FixedMultivariateDistrib],
#'   [distrib_cross3_y()],
#'   [`distrib_grad_y_hess()`][distrib_grad_y_hess.FixedMultivariateDistrib] and
#'   [`distrib_hess_y_hess()`][distrib_hess_y_hess.FixedMultivariateDistrib];
#' - the moments: [`variance()`][variance.FixedMultivariateDistrib];
#' - the multivariate summaries:
#'   [`mv_location()`][mv_location.FixedMultivariateDistrib],
#'   [`mv_sigma()`][mv_sigma.FixedMultivariateDistrib],
#'   [`mv_marginal()`][mv_marginal.FixedMultivariateDistrib],
#'   [`mv_support()`][mv_support.FixedMultivariateDistrib],
#'   [`mv_derived()`][mv_derived.FixedMultivariateDistrib] and
#'   [`mv_reference_draw()`][mv_reference_draw.FixedMultivariateDistrib];
#' - the exactness flags: [expected_hessian_exact()] and [has_exact_deriv4()].
#'
#' Every other generic is inherited from [multivariate_distrib()].
#'
#' @aliases distrib_cross2_y.FixedMultivariateDistrib
#' @aliases distrib_cross_y.FixedMultivariateDistrib
#' @aliases distrib_deriv3.FixedMultivariateDistrib
#' @aliases distrib_deriv4.FixedMultivariateDistrib
#' @aliases distrib_expected_hessian.FixedMultivariateDistrib
#' @aliases distrib_grad_y.FixedMultivariateDistrib
#' @aliases distrib_grad_y_hess.FixedMultivariateDistrib
#' @aliases distrib_gradient.FixedMultivariateDistrib
#' @aliases distrib_hess_y.FixedMultivariateDistrib
#' @aliases distrib_hess_y_hess.FixedMultivariateDistrib
#' @aliases distrib_hessian.FixedMultivariateDistrib
#' @aliases distrib_pdf.FixedMultivariateDistrib
#' @aliases distrib_rng.FixedMultivariateDistrib
#' @aliases mv_derived.FixedMultivariateDistrib
#' @aliases mv_location.FixedMultivariateDistrib
#' @aliases mv_marginal.FixedMultivariateDistrib
#' @aliases mv_reference_draw.FixedMultivariateDistrib
#' @aliases mv_sigma.FixedMultivariateDistrib
#' @aliases mv_support.FixedMultivariateDistrib
#' @aliases variance.FixedMultivariateDistrib
#'
#' @description
#' This is the S7 class of a MULTIVARIATE distribution in which some parameters
#' of the wrapped distribution are held at known values. It follows
#' [FixedContinuousDistrib] in splicing the held values into `theta` and
#' delegating to the parent. It has no distribution function, quantile or atoms,
#' and adds the multivariate generics: [mv_location()], [mv_sigma()],
#' [mv_marginal()], [mv_support()], [mv_reference_draw()] and [mv_derived()] all
#' delegate to the parent at the spliced `theta`, and [mv_derived()] reports the
#' parent's quantities with the Jacobian columns of the fixed parameters
#' removed.
#'
#' The motivating case is a CENTERED PRIOR: holding the mean components of a
#' multivariate family at zero leaves only the matrix parameter, and a random
#' effect is distributed according to a law with this matrix as its only
#' parameter.
#'
#' @details
#' # How every method works
#'
#' The free parameters are the parent's minus the fixed ones, IN THE PARENT'S
#' ORDER. Every method splices the fixed values back into `theta` at their
#' positions and delegates to the parent, so the parent's closed forms are used
#' wherever they exist, and a derivative method then keeps only the components
#' in which every index is a free parameter. Every method of this class
#' delegates its computation to the parent, so the speed and the accuracy of
#' each method are those of the parent's method.
#'
#' # Why the components can be subset by name
#'
#' The free set preserves the parent's ORDER, so a combination of free
#' parameters produces the same name string under the parent's enumeration as
#' under the wrapper's. Subsetting by name therefore cannot pair a component
#' with the wrong indices. Re-parsing a name by splitting on the underscore
#' makes this mistake for a parameter whose own name contains an underscore.
#'
#' # The methods carry no documentation pages
#'
#' The continuous and discrete branches register their methods inside a loop
#' over the two classes, so roxygen has no top-level assignment to which it can
#' attach a block. The multivariate branch registers its methods at the top
#' level and leaves them undocumented to match, because each of them is the
#' same delegation.
#' Read the parent's page for what a method computes and this page for what the
#' wrapper does to it.
#'
#' # Generics that are rejected
#'
#' The multivariate branch sits BESIDE the continuous and discrete ones rather
#' than under either, so the generics that a multivariate family rejects by
#' design, the distribution function and the quantile, are inherited
#' unregistered and go on rejecting through this wrapper.
#'
#' @param parent_distrib The wrapped `multivariate_distrib` object.
#' @param fixed_params A named list of the fixed values, one single finite
#'   number each, strictly inside the corresponding parameter's open domain.
#' @param n_dim The number of coordinates, carried from the parent: fixing a
#'   parameter removes it from the parameter set and leaves the dimension of
#'   the response unchanged.
#' @inheritParams distrib
#'
#' @return An S7 object of class `FixedMultivariateDistrib`, inheriting from
#'   `multivariate_distrib` and from `distrib`. It carries `parent_distrib` and
#'   `fixed_params` beside the parent's properties, `n_dim` included.
#'
#' @seealso [fixed()] to build one, [FixedContinuousDistrib] and
#'   [FixedDiscreteDistrib] for the two siblings, [mvgaussian1_distrib()] for a
#'   parent, and [mv_summary()] for the quantities that a fit of one reports.
#'
#' @examples
#' # A centered two-dimensional gaussian has only the matrix parameter, as the
#' # law of a random effect does.
#' d <- fixed(mvgaussian1_distrib(2), mu1 = 0, mu2 = 0)
#' d@params
#' d
#'
#' theta <- list(sigma_log_L1 = 0, sigma_log_L2 = 0, sigma_L2.1 = 0.5)
#'
#' # The law is the parent's at the spliced theta.
#' y <- rbind(c(0, 0), c(1, -1))
#' all.equal(distrib_pdf(d, y, theta, log = TRUE),
#'           distrib_pdf(mvgaussian1_distrib(2), y,
#'                       c(list(mu1 = 0, mu2 = 0), theta), log = TRUE))
#'
#' # The multivariate generics work through the wrapper.
#' mv_sigma(d, theta)
#' names(mv_derived(d, theta)$value)
#' mv_marginal(d, theta, 1)$distrib@n_dim
#'
#' # And so do the errors that the base class registers.
#' try(distrib_cdf(d, y, theta))
FixedMultivariateDistrib <- S7::new_class("FixedMultivariateDistrib",
  parent = multivariate_distrib,
  properties = list(
    parent_distrib = distrib,
    fixed_params = S7::class_list
  )
)

#' @title Is This a Fixed-Parameter Wrapper
#'
#' @description
#' Returns `TRUE` for a distribution produced by [fixed()], in any of its three
#' forms. [fixed()] uses it to COLLAPSE a `fixed()` of a `fixed()` into one
#' wrapper. The two describe the same law either way, and nesting would leave
#' the inner object's free set to be reconstructed at every call.
#'
#' @param distrib A `distrib` object.
#'
#' @return `TRUE` for a `FixedContinuousDistrib`, a `FixedDiscreteDistrib` or a
#'   `FixedMultivariateDistrib`, `FALSE` otherwise.
#'
#' @seealso [fixed()], which consults this, and [FixedContinuousDistrib] for
#'   the class.
#'
#' @examples
#' c(plain = distributions7:::is_fixed(gaussian1_distrib()),
#'   fixed = distributions7:::is_fixed(fixed(gaussian1_distrib(), mu = 0)),
#'   multivariate = distributions7:::is_fixed(
#'     fixed(mvgaussian1_distrib(2), mu1 = 0)))
#'
#' # Two calls therefore collapse into one wrapper that holds both values.
#' d <- fixed(fixed(gaussian1_distrib(), mu = 0), sigma = 1)
#' c(class = class(d)[1], held = paste(names(d@fixed_params), collapse = ", "))
#'
#' @keywords internal
is_fixed <- function(distrib) {
  S7::S7_inherits(distrib, FixedContinuousDistrib) ||
    S7::S7_inherits(distrib, FixedDiscreteDistrib) ||
    S7::S7_inherits(distrib, FixedMultivariateDistrib)
}

#' @title Splice the Fixed Values Back Into a Full Parameter List
#'
#' @description
#' Combines the wrapper's free `theta` with its fixed values into the full
#' parameter list that the PARENT expects, in the parent's own order. Every
#' method of every `Fixed*` class begins with this call and then delegates, so
#' it is the single point at which the two parameter sets are reconciled.
#'
#' @param distrib A `FixedContinuousDistrib`, `FixedDiscreteDistrib` or
#'   `FixedMultivariateDistrib` object.
#' @param theta A named list of the FREE parameters, already aligned. A
#'   fully-fixed wrapper takes `list()`.
#'
#' @return A named list of the parent's parameters, complete and in the
#'   parent's order.
#'
#' @seealso [fixed()] for the wrapper and [FixedContinuousDistrib] for what
#'   the methods do with the result.
#'
#' @examples
#' d <- fixed(gaussian1_distrib(), mu = 0)
#' str(distributions7:::fixed_full_theta(d, list(sigma = 2)))
#'
#' # The parent is then called at exactly this list.
#' all.equal(distrib_pdf(d, c(-1, 1), list(sigma = 2)),
#'           distrib_pdf(gaussian1_distrib(), c(-1, 1),
#'                       distributions7:::fixed_full_theta(d, list(sigma = 2))))
#'
#' @keywords internal
fixed_full_theta <- function(distrib, theta) {
  theta <- align_theta(distrib, theta)
  free <- distrib@params
  theta <- theta[seq_len(distrib@n_params)]
  names(theta) <- free

  parent <- distrib@parent_distrib
  out <- vector("list", parent@n_params)
  names(out) <- parent@params
  out[free] <- theta
  fp <- distrib@fixed_params
  out[names(fp)] <- fp
  out
}

# ---------------------------------------------------------------------------
# Method registration.
#
# Every method has the same one-line behavior -- splice and delegate -- so
# they are registered in a loop over the two classes. The functions do not
# read the loop variable, so no closure capture is involved. The derivative
# methods subset the parent's result by the names generated from the free
# parameter set: a combination of free parameters produces the same name
# string under the parent's enumeration as under the wrapper's, because the
# free set preserves the parent's order, so subsetting by name cannot pair a
# component with the wrong indices -- the mistake that re-parsing names by
# splitting on the underscore commits for a parameter whose own name contains
# one.
# ---------------------------------------------------------------------------

# The third and fourth response derivatives delegate like the orders below.
# Without these a fixed family fell back to the base class's stencil on the
# log-density, so fixed(student_t1_distrib(), mu = 0) -- a heavy-tailed prior --
# took a difference where its parent has the closed form. Continuous only: a
# discrete family has no response derivative to delegate to.
S7::method(distrib_deriv3_y, FixedContinuousDistrib) <- function(distrib, y,
                                                                 theta, ...) {
  distrib_deriv3_y(distrib@parent_distrib, y, fixed_full_theta(distrib, theta),
                   ...)
}

S7::method(distrib_deriv4_y, FixedContinuousDistrib) <- function(distrib, y,
                                                                 theta, ...) {
  distrib_deriv4_y(distrib@parent_distrib, y, fixed_full_theta(distrib, theta),
                   ...)
}

# and the mixed third-response derivative, subset to the free parameters as
# cross_y and cross2_y are; a heavy-tailed prior reads it through penalties7
S7::method(distrib_cross3_y, FixedContinuousDistrib) <- function(distrib, y,
                                                                 theta,
                                                                 scale = c("parameter", "link"),
                                                                 ...) {
  res <- distrib_cross3_y(distrib@parent_distrib, y,
                          fixed_full_theta(distrib, theta),
                          scale = "parameter", ...)
  res[distrib@params]
}

for (.fixed_cls in list(FixedContinuousDistrib, FixedDiscreteDistrib)) {
  S7::method(distrib_pdf, .fixed_cls) <- function(distrib, y, theta, log = FALSE, ...) {
    distrib_pdf(distrib@parent_distrib, y, fixed_full_theta(distrib, theta),
      log = log
    )
  }

  S7::method(distrib_cdf, .fixed_cls) <- function(distrib, q, theta, ...) {
    distrib_cdf(distrib@parent_distrib, q, fixed_full_theta(distrib, theta), ...)
  }

  S7::method(distrib_quantile, .fixed_cls) <- function(distrib, p, theta, ...) {
    distrib_quantile(distrib@parent_distrib, p, fixed_full_theta(distrib, theta), ...)
  }

  S7::method(distrib_rng, .fixed_cls) <- function(distrib, n, theta, ...) {
    distrib_rng(distrib@parent_distrib, n, fixed_full_theta(distrib, theta), ...)
  }

  S7::method(distrib_atoms, .fixed_cls) <- function(distrib, theta, ...) {
    distrib_atoms(distrib@parent_distrib, fixed_full_theta(distrib, theta), ...)
  }

  S7::method(distrib_grad_y, .fixed_cls) <- function(distrib, y, theta, ...) {
    distrib_grad_y(distrib@parent_distrib, y, fixed_full_theta(distrib, theta), ...)
  }

  S7::method(distrib_hess_y, .fixed_cls) <- function(distrib, y, theta, ...) {
    distrib_hess_y(distrib@parent_distrib, y, fixed_full_theta(distrib, theta), ...)
  }

  S7::method(distrib_cross_y, .fixed_cls) <- function(distrib, y, theta,
                                                      scale = c("parameter", "link"),
                                                      ...) {
    res <- distrib_cross_y(distrib@parent_distrib, y,
      fixed_full_theta(distrib, theta),
      scale = "parameter", ...
    )
    res[distrib@params]
  }

  S7::method(distrib_cross2_y, .fixed_cls) <- function(distrib, y, theta,
                                                       scale = c("parameter", "link"),
                                                       ...) {
    res <- distrib_cross2_y(distrib@parent_distrib, y,
      fixed_full_theta(distrib, theta),
      scale = "parameter", ...
    )
    res[distrib@params]
  }

  S7::method(distrib_grad_y_hess, .fixed_cls) <- function(distrib, y, theta,
                                                          scale = c("parameter", "link"),
                                                          ...) {
    res <- distrib_grad_y_hess(distrib@parent_distrib, y,
      fixed_full_theta(distrib, theta),
      scale = "parameter", ...
    )
    res[hess_names(distrib@params)]
  }

  S7::method(distrib_hess_y_hess, .fixed_cls) <- function(distrib, y, theta,
                                                          scale = c("parameter", "link"),
                                                          ...) {
    res <- distrib_hess_y_hess(distrib@parent_distrib, y,
      fixed_full_theta(distrib, theta),
      scale = "parameter", ...
    )
    res[hess_names(distrib@params)]
  }

  S7::method(distrib_gradient, .fixed_cls) <- function(distrib, y, theta,
                                                       scale = c("parameter", "link"),
                                                       ...) {
    res <- distrib_gradient(distrib@parent_distrib, y,
      fixed_full_theta(distrib, theta),
      scale = "parameter", ...
    )
    res[distrib@params]
  }

  S7::method(distrib_hessian, .fixed_cls) <- function(distrib, y, theta,
                                                      scale = c("parameter", "link"),
                                                      ...) {
    res <- distrib_hessian(distrib@parent_distrib, y,
      fixed_full_theta(distrib, theta),
      scale = "parameter", ...
    )
    res[hess_names(distrib@params)]
  }

  S7::method(distrib_expected_hessian, .fixed_cls) <- function(distrib, y, theta,
                                                               scale = c("parameter", "link"),
                                                               approx = c("opg", "bartlett", "integrate", "mc"),
                                                               nsim = 10000, ...) {
    res <- distrib_expected_hessian(distrib@parent_distrib, y,
      fixed_full_theta(distrib, theta),
      scale = "parameter", approx = approx, nsim = nsim, ...
    )
    res[hess_names(distrib@params)]
  }

  # A METHOD REGISTERED ON THIS CLASS IS NOT BY ITSELF A CLOSED FORM. The
  # default `expected_hessian_exact()` reads the OWNER of the registered
  # `distrib_expected_hessian` method, and that owner is `.fixed_cls` here
  # regardless of whether the wrapped family has one: `fixed(pig1_distrib(),
  # sigma = 0.8)` would otherwise answer TRUE and a consumer -- statmodels7's
  # `vcov()` among them -- would report standard errors from an approximation
  # believing them exact. Exactness is the PARENT's.
  S7::method(expected_hessian_exact, .fixed_cls) <- function(x, ...) {
    expected_hessian_exact(x@parent_distrib)
  }
  S7::method(has_exact_deriv4, .fixed_cls) <- function(x, ...) {
    has_exact_deriv4(x@parent_distrib)
  }


  S7::method(distrib_deriv3, .fixed_cls) <- function(distrib, y, theta,
                                                     expected = FALSE,
                                                     scale = c("parameter", "link"),
                                                     approx = c("integrate", "bartlett", "mc", "opg"),
                                                     nsim = 10000, ...) {
    res <- distrib_deriv3(distrib@parent_distrib, y,
      fixed_full_theta(distrib, theta),
      expected = expected, scale = "parameter", approx = approx, nsim = nsim, ...
    )
    res[deriv_names(distrib@params, 3L)]
  }

  S7::method(distrib_deriv4, .fixed_cls) <- function(distrib, y, theta,
                                                     expected = FALSE,
                                                     scale = c("parameter", "link"),
                                                     approx = c("integrate", "bartlett", "mc", "opg"),
                                                     nsim = 10000, ...) {
    res <- distrib_deriv4(distrib@parent_distrib, y,
      fixed_full_theta(distrib, theta),
      expected = expected, scale = "parameter", approx = approx, nsim = nsim, ...
    )
    res[deriv_names(distrib@params, 4L)]
  }

  S7::method(distrib_grad_cdf, .fixed_cls) <- function(distrib, q, theta,
                                                       lower.tail = TRUE, log = TRUE,
                                                       ...) {
    res <- distrib_grad_cdf(distrib@parent_distrib, q,
      fixed_full_theta(distrib, theta),
      lower.tail = lower.tail, log = log, ...
    )
    res[distrib@params]
  }

  S7::method(distrib_hess_cdf, .fixed_cls) <- function(distrib, q, theta,
                                                       lower.tail = TRUE, log = TRUE,
                                                       ...) {
    res <- distrib_hess_cdf(distrib@parent_distrib, q,
      fixed_full_theta(distrib, theta),
      lower.tail = lower.tail, log = log, ...
    )
    res[hess_names(distrib@params)]
  }

  # The moments delegate so that a parent with a closed form keeps it; the
  # law is the parent's law at the full parameter vector, so nothing changes
  # but where theta comes from.
  S7::method(mean, .fixed_cls) <- function(x, theta, ...) {
    mean(x@parent_distrib, fixed_full_theta(x, theta), ...)
  }

  S7::method(variance, .fixed_cls) <- function(x, theta, ...) {
    variance(x@parent_distrib, fixed_full_theta(x, theta), ...)
  }

  S7::method(std_dev, .fixed_cls) <- function(x, theta, ...) {
    std_dev(x@parent_distrib, fixed_full_theta(x, theta), ...)
  }

  S7::method(skewness, .fixed_cls) <- function(x, theta, ...) {
    skewness(x@parent_distrib, fixed_full_theta(x, theta), ...)
  }

  S7::method(kurtosis, .fixed_cls) <- function(x, theta, ...) {
    kurtosis(x@parent_distrib, fixed_full_theta(x, theta), ...)
  }

  S7::method(print, .fixed_cls) <- function(x, ...) {
    # The base print iterates over the free parameters, which works down to
    # one of them; with none it computes the width of an empty name set, so
    # that case prints its own header instead. print() is a base generic and a
    # method registered on it is an S3 method, so the parent class's print is
    # reached with NextMethod(); super() only works inside S7 generics, and
    # S7::method() only retrieves from them.
    if (x@n_params > 0L) {
      NextMethod()
    } else {
      d_name <- gsub("(^|[[:space:]])([[:alpha:]])", "\\1\\U\\2",
        x@distrib_name,
        perl = TRUE
      )
      d_type <- if (S7::S7_inherits(x, continuous_distrib)) "Continuous" else "Discrete"
      cat(sprintf("Distribution: %s\n", d_name))
      cat(sprintf("Type:         %s\n", d_type))
      cat(sprintf("Dimensions:   %s\n", x@dimension))
      cat("\nParameters:   none free\n")
    }
    fp <- x@fixed_params
    cat("\nFixed:\n")
    for (nm in names(fp)) {
      cat(sprintf("  %s = %s\n", nm, format(fp[[nm]])))
    }
    invisible(x)
  }
}
rm(.fixed_cls)

#' The Derivatives of a Fixed Family's Expected Information
#'
#' @description
#' Returns the parent family's first (`k = 1`) or second (`k = 2`)
#' derivatives of the expected information at the full parameter vector,
#' subset to the free parameters, on the parameter scale. The function is
#' registered through [register_dexpected()], which carries the derivatives to
#' the link scale with the wrapper's links.
#'
#' @param d A fixed family, from [fixed()].
#' @param y The response, read for its length.
#' @param th The free parameters.
#' @param k The order, `1L` or `2L`.
#' @param t The thread count passed to the parent's methods.
#'
#' @return A named list keyed as [dexpected_names()] (`k = 1`) or
#'   [d2expected_names()] (`k = 2`) over the free parameters.
#'
#' @keywords internal
fixed_dexpected <- function(d, y, th, k, t) {
  parent <- d@parent_distrib
  full <- fixed_full_theta(d, th)
  if (k == 1L) {
    res <- distrib_dexpected_hessian(parent, y, full, scale = "parameter",
                                     threads = t)
    return(res[dexpected_names(d@params)])
  }
  res <- distrib_d2expected_hessian(parent, y, full, scale = "parameter",
                                    threads = t)
  res[d2expected_names(d@params)]
}
register_dexpected(FixedContinuousDistrib, fixed_dexpected)
register_dexpected(FixedDiscreteDistrib, fixed_dexpected)

# ---------------------------------------------------------------------------
# The multivariate branch.
#
# Registered apart from the loop above because the two sets of generics differ:
# a multivariate family has no distribution function, no quantile and no atoms,
# and it has the mv_* accessors instead. The generics it rejects are NOT
# registered here, so the rejection is inherited from multivariate_distrib and
# keeps its message.
#
# The derivative methods subset by name from the free set, exactly as the
# univariate ones do; a multivariate family's parameters are already flattened
# into scalars, so the enumeration that generates the component names is the
# same one and no special case is involved.
# ---------------------------------------------------------------------------

S7::method(distrib_pdf, FixedMultivariateDistrib) <-
  function(distrib, y, theta, log = FALSE, ...) {
    distrib_pdf(distrib@parent_distrib, y, fixed_full_theta(distrib, theta),
                log = log)
  }

S7::method(distrib_rng, FixedMultivariateDistrib) <-
  function(distrib, n, theta, ...) {
    distrib_rng(distrib@parent_distrib, n, fixed_full_theta(distrib, theta),
                ...)
  }

S7::method(distrib_grad_y, FixedMultivariateDistrib) <-
  function(distrib, y, theta, ...) {
    distrib_grad_y(distrib@parent_distrib, y, fixed_full_theta(distrib, theta),
                   ...)
  }

S7::method(distrib_hess_y, FixedMultivariateDistrib) <-
  function(distrib, y, theta, ...) {
    distrib_hess_y(distrib@parent_distrib, y, fixed_full_theta(distrib, theta),
                   ...)
  }

# the third and fourth response derivatives and the mixed one, which a
# multivariate prior whose Hessian moves with the coefficients is read through
S7::method(distrib_deriv3_y, FixedMultivariateDistrib) <-
  function(distrib, y, theta, ...) {
    distrib_deriv3_y(distrib@parent_distrib, y,
                     fixed_full_theta(distrib, theta), ...)
  }

S7::method(distrib_deriv4_y, FixedMultivariateDistrib) <-
  function(distrib, y, theta, ...) {
    distrib_deriv4_y(distrib@parent_distrib, y,
                     fixed_full_theta(distrib, theta), ...)
  }

S7::method(distrib_cross3_y, FixedMultivariateDistrib) <-
  function(distrib, y, theta, scale = c("parameter", "link"), ...) {
    out <- distrib_cross3_y(distrib@parent_distrib, y,
                            fixed_full_theta(distrib, theta),
                            scale = "parameter", ...)
    out[distrib@params]
  }

# Every derivative below reads the parent on the PARAMETER scale, as the
# univariate methods do: the generic of this class applies the link once, so a
# parent asked for the link scale would apply it a second time. With every free
# parameter on the identity link the two coincide, which is how the defect went
# unseen; nu of the multivariate t or phi of the Dirichlet left free exposed it.
S7::method(distrib_cross_y, FixedMultivariateDistrib) <-
  function(distrib, y, theta, scale = c("parameter", "link"), ...) {
    out <- distrib_cross_y(distrib@parent_distrib, y,
                           fixed_full_theta(distrib, theta),
                           scale = "parameter", ...)
    out[distrib@params]
  }

S7::method(distrib_cross2_y, FixedMultivariateDistrib) <-
  function(distrib, y, theta, scale = c("parameter", "link"), ...) {
    out <- distrib_cross2_y(distrib@parent_distrib, y,
                            fixed_full_theta(distrib, theta),
                            scale = "parameter", ...)
    out[distrib@params]
  }

S7::method(distrib_grad_y_hess, FixedMultivariateDistrib) <-
  function(distrib, y, theta, scale = c("parameter", "link"), ...) {
    out <- distrib_grad_y_hess(distrib@parent_distrib, y,
                               fixed_full_theta(distrib, theta),
                               scale = "parameter", ...)
    out[hess_names(distrib@params)]
  }

S7::method(distrib_hess_y_hess, FixedMultivariateDistrib) <-
  function(distrib, y, theta, scale = c("parameter", "link"), ...) {
    out <- distrib_hess_y_hess(distrib@parent_distrib, y,
                               fixed_full_theta(distrib, theta),
                               scale = "parameter", ...)
    out[hess_names(distrib@params)]
  }

S7::method(distrib_gradient, FixedMultivariateDistrib) <-
  function(distrib, y, theta, scale = c("parameter", "link"), ...) {
    out <- distrib_gradient(distrib@parent_distrib, y,
                            fixed_full_theta(distrib, theta),
                            scale = "parameter", ...)
    out[distrib@params]
  }

S7::method(distrib_hessian, FixedMultivariateDistrib) <-
  function(distrib, y, theta, scale = c("parameter", "link"), ...) {
    out <- distrib_hessian(distrib@parent_distrib, y,
                           fixed_full_theta(distrib, theta),
                           scale = "parameter", ...)
    out[hess_names(distrib@params)]
  }

S7::method(distrib_expected_hessian, FixedMultivariateDistrib) <-
  function(distrib, y, theta, scale = c("parameter", "link"),
           approx = c("opg", "bartlett", "integrate", "mc"),
           nsim = 10000, ...) {
    out <- distrib_expected_hessian(distrib@parent_distrib, y,
                                    fixed_full_theta(distrib, theta),
                                    scale = "parameter", approx = approx,
                                    nsim = nsim, ...)
    out[hess_names(distrib@params)]
  }

# Exactness is the parent's, not this class's own -- see the univariate
# .fixed_cls registration above for why the default predicate gets this wrong
# without it.
S7::method(expected_hessian_exact, FixedMultivariateDistrib) <- function(x, ...) {
  expected_hessian_exact(x@parent_distrib)
}

S7::method(has_exact_deriv4, FixedMultivariateDistrib) <- function(x, ...) {
  has_exact_deriv4(x@parent_distrib)
}


S7::method(distrib_deriv3, FixedMultivariateDistrib) <-
  function(distrib, y, theta, expected = FALSE,
           scale = c("parameter", "link"),
           approx = c("integrate", "bartlett", "mc", "opg"),
           nsim = 10000, ...) {
    out <- distrib_deriv3(distrib@parent_distrib, y,
                          fixed_full_theta(distrib, theta),
                          expected = expected, scale = "parameter",
                          approx = approx, nsim = nsim, ...)
    out[deriv_names(distrib@params, 3L)]
  }

S7::method(distrib_deriv4, FixedMultivariateDistrib) <-
  function(distrib, y, theta, expected = FALSE,
           scale = c("parameter", "link"),
           approx = c("integrate", "bartlett", "mc", "opg"),
           nsim = 10000, ...) {
    out <- distrib_deriv4(distrib@parent_distrib, y,
                          fixed_full_theta(distrib, theta),
                          expected = expected, scale = "parameter",
                          approx = approx, nsim = nsim, ...)
    out[deriv_names(distrib@params, 4L)]
  }

S7::method(mv_location, FixedMultivariateDistrib) <- function(distrib, theta) {
  mv_location(distrib@parent_distrib, fixed_full_theta(distrib, theta))
}

S7::method(mv_sigma, FixedMultivariateDistrib) <- function(distrib, theta) {
  mv_sigma(distrib@parent_distrib, fixed_full_theta(distrib, theta))
}

S7::method(mv_marginal, FixedMultivariateDistrib) <-
  function(distrib, theta, which, ...) {
    mv_marginal(distrib@parent_distrib, fixed_full_theta(distrib, theta),
                which, ...)
  }

S7::method(mv_support, FixedMultivariateDistrib) <- function(distrib, theta, ...) {
  mv_support(distrib@parent_distrib, fixed_full_theta(distrib, theta), ...)
}

S7::method(mv_reference_draw, FixedMultivariateDistrib) <-
  function(distrib, theta, n, ...) {
    mv_reference_draw(distrib@parent_distrib, fixed_full_theta(distrib, theta),
                      n, ...)
  }

S7::method(mv_derived, FixedMultivariateDistrib) <-
  function(distrib, theta, ...) {
    # WITHOUT this the wrapper falls to the base method, which reports the
    # distinct entries of the covariance rather than the standard deviations
    # and correlations the family declares -- a centered prior would then be
    # read on a scale its own family does not use.
    out <- mv_derived(distrib@parent_distrib, fixed_full_theta(distrib, theta),
                      ...)
    if (is.null(out)) return(NULL)
    # the quantities are the parent's; the Jacobian keeps the columns of the
    # free parameters alone, the fixed ones having no derivative to report
    keep <- match(distrib@params, distrib@parent_distrib@params)
    out$jacobian <- out$jacobian[, keep, drop = FALSE]
    out
  }

S7::method(mean, FixedMultivariateDistrib) <- function(x, theta, ...) {
  mean(x@parent_distrib, fixed_full_theta(x, theta), ...)
}

S7::method(variance, FixedMultivariateDistrib) <- function(x, theta, ...) {
  variance(x@parent_distrib, fixed_full_theta(x, theta), ...)
}

S7::method(print, FixedMultivariateDistrib) <- function(x, ...) {
  # as in the univariate wrappers: the base print iterates over the free
  # parameters and cannot size an empty set, so that case prints its own header
  if (x@n_params > 0L) {
    NextMethod()
  } else {
    cat(sprintf("Distribution: %s\n", x@distrib_name))
    cat("Type:         Multivariate\n")
    cat(sprintf("Dimensions:   %d\n", x@n_dim))
    cat("\nParameters:   none free\n")
  }
  fp <- x@fixed_params
  cat("\nFixed:\n")
  for (nm in names(fp)) {
    cat(sprintf("  %s = %s\n", nm, format(fp[[nm]])))
  }
  invisible(x)
}

#' @title Fix Parameters of a Distribution at Known Values
#'
#' @description
#' Returns the distribution obtained by holding some parameters of `distrib` at
#' known values, leaving the others to be supplied and estimated. It is the
#' only wrapper in this package that REMOVES parameters, and it derives
#' nothing: the density is the parent's at the reassembled vector and the
#' derivatives are the parent's components among the free indices.
#'
#' @details
#' # The construction
#'
#' Splitting the parent's parameters into a fixed part \eqn{\theta_C = c} and a
#' free part \eqn{\theta_F},
#' \deqn{f_{\mathrm{fix}}(y; \theta_F) = f(y; \theta_F, c),
#'   \qquad \ell^{(i_1 \cdots i_k)}_{\mathrm{fix}} = \ell^{(i_1 \cdots i_k)},
#'   \quad i_1, \dots, i_k \in F.}
#' `theta` carries only the free parameters, every generic returns what the
#' parent returns at the full vector, and a derivative is the parent's
#' restricted to the free indices: a subvector of the score, a submatrix of the
#' Hessian, sub-arrays at orders three and four. The wrapper does not recompute
#' anything and does not change the normalizing constant, and [fit_distrib()]
#' estimates the free parameters and reports standard errors and intervals only
#' for those parameters.
#'
#' # What is accepted
#'
#' Fixed values are single finite numbers, strictly inside the OPEN domain of
#' their parameter. A multivariate family's mean is the vector parameter
#' `mu1, ..., mup`, and the name `mu` fixes all of it: `mu = 0` holds every
#' coordinate at zero, and `mu = c(0, 1, 2)` holds each at its own value. The
#' same holds for any parameter named as a stem followed by the indices
#' `1, ..., k`. Fixing a parameter of a distribution that is already a
#' fixed-parameter wrapper COLLAPSES the two into one wrapper around the
#' original parent. Fixing a WRAPPER's own parameter is allowed and useful:
#' `fixed(zero_inflated(d), zi = 0.3)` is a zero-inflated model with a known
#' inflation rate. Fixing every parameter is allowed too and gives a fully
#' known distribution with an empty parameter set. A call with no named value
#' signals an error, because the result would be the parent unchanged and the
#' missing argument would go unreported.
#'
#' # Purpose
#'
#' The wrapper builds priors. `fixed(gaussian1_distrib(), mu = 0)` is the ridge
#' penalty with a free scale, `fixed(laplace2_distrib(), mu = 0)` is the
#' lasso, and `fixed(mvgaussian1_distrib(p), mu = 0)` is the distribution of a
#' random effect. `fixed(folded(gaussian1_distrib()), mu = 0)` is the
#' half-normal.
#'
#' The per-parameter smoothness declaration travels with the free parameters,
#' so fixing the location of a Laplace leaves a distribution whose remaining
#' parameter is smooth.
#'
#' @section Notation:
#' \eqn{f} is the parent's density, \eqn{\theta_C = c} the fixed parameters,
#' \eqn{\theta_F} the free ones, \eqn{F} their index set and
#' \eqn{\ell^{(i_1\cdots i_k)}} a derivative of the log-density in the
#' parameters named.
#'
#' @param distrib The distribution whose parameters are to be fixed, inheriting
#'   from `continuous_distrib`, `discrete_distrib` or `multivariate_distrib`.
#'   A parameter whose name is a prefix of `distrib` (the shape `d` of
#'   [gengamma1_distrib()] and [gengamma2_distrib()]) is fixed only when
#'   `distrib` is named, as in `fixed(distrib = gengamma1_distrib(), d = 2)`;
#'   otherwise partial matching assigns the value to `distrib` and an error is
#'   signalled.
#' @param ... The fixed values, named after the parameters that they fix, as in
#'   `fixed(gaussian1_distrib(), mu = 0)`. Each must be a single finite number
#'   strictly inside its parameter's domain, and each name must be a parameter
#'   of `distrib`, or the stem of a vector parameter such as a multivariate
#'   mean, whose value is then one number for every coordinate or one number
#'   per coordinate. A name that is not, a value outside the domain, a value that
#'   is not a single number, and an empty `...` are each rejected with an error
#'   that names the condition that failed.
#'
#' @return An S7 object of class [FixedContinuousDistrib],
#'   [FixedDiscreteDistrib] or [FixedMultivariateDistrib], matching the
#'   parent's branch. Its `params` are the parent's less the fixed names in the
#'   parent's order, `n_params` falls by as many, `fixed_params` holds the
#'   values, and `distrib_name` is the parent's with the held values in
#'   brackets.
#'
#' @seealso [zero_inflated()], [truncated()], [transformation()] and
#'   [folded()] for the other wrappers, and [FixedContinuousDistrib] for what
#'   the methods do.
#'
#' @examples
#' # This multivariate gaussian has its whole mean held at zero.
#' fixed(mvgaussian1_distrib(3), mu = 0)@params
#'
#' # Fixing the mean of a gaussian leaves only sigma as a free parameter.
#' d <- fixed(gaussian1_distrib(), mu = 0)
#' d@params
#' theta <- list(sigma = 2)
#' distrib_pdf(d, c(-1, 0, 1), theta)
#'
#' # The score is the corresponding component of the parent's, unchanged.
#' full <- distrib_gradient(gaussian1_distrib(), c(-1, 0, 1),
#'                          list(mu = 0, sigma = 2))
#' all.equal(distrib_gradient(d, c(-1, 0, 1), theta)$sigma, full$sigma)
#'
#' # The lasso prior is a Laplace in its rate, centered at zero.
#' lasso <- fixed(laplace2_distrib(), mu = 0)
#' lasso@params
#' b <- c(-1, 0.5, 2)
#' all.equal(-distrib_pdf(lasso, b, list(lambda = 2), log = TRUE),
#'           2 * abs(b) - log(2 / 2))
#'
#' # A wrapper's own parameter can be held.
#' fixed(zero_inflated(poisson_distrib()), zi = 0.3)@params
#'
#' # Two calls collapse into one wrapper, and fixing every parameter is allowed.
#' fixed(fixed(gaussian1_distrib(), mu = 0), sigma = 1)@fixed_params
#'
#' # The four errors below each name the condition that failed.
#' try(fixed(gaussian1_distrib(), nope = 1))
#' try(fixed(gaussian1_distrib(), sigma = -1))
#' try(fixed(gaussian1_distrib(), mu = c(0, 1)))
#' try(fixed(gaussian1_distrib()))
#'
#' @export
fixed <- function(distrib, ...) {
  is_distrib <- function(x) {
    S7::S7_inherits(x, continuous_distrib) ||
      S7::S7_inherits(x, discrete_distrib) ||
      S7::S7_inherits(x, multivariate_distrib)
  }
  # a value named by a prefix of `distrib` (the shape d of gengamma1 and
  # gengamma2) is matched to `distrib` by partial matching, and the
  # distribution itself lands in `...`
  if (!is_distrib(distrib) && any(vapply(list(...), is_distrib, logical(1)))) {
    stop(paste0(
      "fixed() received a distribution in '...' and a value in place of ",
      "'distrib'.\n  A value named by a prefix of 'distrib', such as d = 2 for ",
      "the shape d of gengamma1\n  and gengamma2, is matched to 'distrib' by ",
      "partial matching. Name the distribution,\n  as in ",
      "fixed(distrib = gengamma1_distrib(), d = 2)."),
      call. = FALSE)
  }
  if (!is_distrib(distrib)) {
    stop(paste0("Input must inherit from 'discrete_distrib', ",
                "'continuous_distrib' or 'multivariate_distrib'."),
      call. = FALSE
    )
  }

  fix <- fixed_expand_stems(list(...), distrib@params)
  if (length(fix) == 0L) {
    stop(paste0(
      "fixed() needs at least one named value, as in fixed(d, mu = 0). With\n",
      "  none, the result would be the parent distribution unchanged, and the\n",
      "  missing argument would go unreported."
    ), call. = FALSE)
  }
  nms <- names(fix)
  if (is.null(nms) || any(!nzchar(nms))) {
    stop("Every fixed value must be named after the parameter that it fixes.",
      call. = FALSE
    )
  }
  if (anyDuplicated(nms)) {
    stop(sprintf(
      "Parameter '%s' is fixed more than once in the same call.",
      nms[duplicated(nms)][1L]
    ), call. = FALSE)
  }

  # Collapse a fixed() of a fixed(): one wrapper around the original parent,
  # carrying both sets of values. A parameter already fixed is caught below
  # by the membership check, since it is no longer among the free parameters.
  inherited <- list()
  if (is_fixed(distrib)) {
    unknown <- setdiff(nms, distrib@params)
    if (length(unknown)) {
      stop(sprintf(
        "'%s' is not a free parameter of '%s'. The free parameters are %s.",
        unknown[1L], distrib@distrib_name,
        paste(distrib@params, collapse = ", ")
      ), call. = FALSE)
    }
    inherited <- distrib@fixed_params
    distrib <- distrib@parent_distrib
  }

  unknown <- setdiff(nms, distrib@params)
  if (length(unknown)) {
    stop(sprintf(
      "'%s' is not a parameter of '%s'. The parameters are %s.",
      unknown[1L], distrib@distrib_name,
      paste(distrib@params, collapse = ", ")
    ), call. = FALSE)
  }

  for (nm in nms) {
    v <- fix[[nm]]
    if (!is.numeric(v) || length(v) != 1L || !is.finite(v)) {
      stop(sprintf(
        "The value fixing '%s' must be a single finite number.", nm
      ), call. = FALSE)
    }
    b <- distrib@params_bounds[[nm]]
    # The domains are open intervals, the convention align_theta() enforces
    # for a free parameter; a fixed one obeys the same rule.
    if (v <= b[1] || v >= b[2]) {
      stop(sprintf(
        "The value fixing '%s' (%s) is outside its open domain (%s, %s).",
        nm, format(v), format(b[1]), format(b[2])
      ), call. = FALSE)
    }
    fix[[nm]] <- as.numeric(v)
  }

  all_fixed <- c(inherited, fix)
  # Keep the parent's ordering in the record, so the printed name is stable
  # whatever order the calls arrived in.
  all_fixed <- all_fixed[intersect(distrib@params, names(all_fixed))]
  free <- setdiff(distrib@params, names(all_fixed))

  smooth <- param_smoothness(distrib)[free]
  if (length(free) == 0L) smooth <- logical(0)

  # No spaces in the label: the print method capitalises the first letter
  # after every space in the distribution's name, and a parameter name must
  # not come out as "Sigma". Same convention as truncated's "[lower=0]".
  label <- paste(
    vapply(
      names(all_fixed),
      function(nm) sprintf("%s=%s", nm, format(all_fixed[[nm]])),
      character(1)
    ),
    collapse = ","
  )

  common <- list(
    parent_distrib = distrib,
    fixed_params = all_fixed,
    distrib_name = sprintf("fixed %s [%s]", distrib@distrib_name, label),
    dimension = distrib@dimension,
    bounds = distrib@bounds,
    params = free,
    params_interpretation = distrib@params_interpretation[free],
    n_params = length(free),
    params_bounds = distrib@params_bounds[free],
    link_params = distrib@link_params[free],
    params_smooth = smooth
  )

  if (S7::S7_inherits(distrib, multivariate_distrib)) {
    common$n_dim <- distrib@n_dim
    do.call(FixedMultivariateDistrib, common)
  } else if (S7::S7_inherits(distrib, discrete_distrib)) {
    do.call(FixedDiscreteDistrib, common)
  } else {
    do.call(FixedContinuousDistrib, common)
  }
}


#' @title Expand the Stem of a Vector Parameter in fixed()
#'
#' @description
#' Replaces a name that is not a parameter but is the stem of the parameters
#' `stem1, ..., stemk` by those parameters, recycling a single value or taking
#' one value per coordinate. With this expansion, `fixed(mvgaussian1_distrib(3),
#' mu = 0)` holds the whole mean at zero.
#'
#' @param fix The named list of values passed to [fixed()].
#' @param params The parameter names of the distribution that is being fixed.
#'
#' @return The list with every stem expanded into its coordinates. A name that
#'   is a parameter, or matches no stem, is returned unchanged for [fixed()]
#'   to check.
#'
#' @examples
#' distributions7:::fixed_expand_stems(list(mu = 0, sigma = 1),
#'                                     c("mu1", "mu2", "sigma"))
#' distributions7:::fixed_expand_stems(list(mu = c(1, 2)), c("mu1", "mu2"))
#'
#' @keywords internal
fixed_expand_stems <- function(fix, params) {
  nms <- names(fix)
  if (is.null(nms)) return(fix)
  out <- list()
  for (i in seq_along(fix)) {
    nm <- nms[i]
    v <- fix[[i]]
    if (!nzchar(nm) || nm %in% params) {
      out[length(out) + 1L] <- list(v)
      names(out)[length(out)] <- nm
      next
    }
    k <- 0L
    while (paste0(nm, k + 1L) %in% params) k <- k + 1L
    if (k == 0L) {
      out[length(out) + 1L] <- list(v)
      names(out)[length(out)] <- nm
      next
    }
    if (!is.numeric(v) || !(length(v) %in% c(1L, k))) {
      stop(sprintf(paste0(
        "The value fixing '%s' must be one number for every coordinate or %d ",
        "numbers, one per coordinate %s1, ..., %s%d."
      ), nm, k, nm, nm, k), call. = FALSE)
    }
    v <- rep_len(v, k)
    for (j in seq_len(k)) {
      out[length(out) + 1L] <- list(v[[j]])
      names(out)[length(out)] <- paste0(nm, j)
    }
  }
  out
}

