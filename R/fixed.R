#' @include distrib.R generics.R utility_functions.R moments.R cross_derivatives.R cross2_derivatives.R cross_theta2_derivatives.R mv_summary.R
NULL

#' @title S7 Class for Distributions With Fixed Parameters (Continuous)
#' @name FixedContinuousDistrib
#'
#' @section Methods:
#' `fixed()` registers 22 methods on this class:
#' `distrib_atoms()`, `distrib_cdf()`, `distrib_cross2_y()`, `distrib_cross_y()`, `distrib_deriv3()`, `distrib_deriv4()`, `distrib_expected_hessian()`, `distrib_grad_cdf()`, `distrib_grad_y()`, `distrib_grad_y_hess()`, `distrib_gradient()`, `distrib_hess_cdf()`, `distrib_hess_y()`, `distrib_hess_y_hess()`, `distrib_hessian()`, `distrib_pdf()`, `distrib_quantile()`, `distrib_rng()`, `kurtosis()`, `skewness()`, `std_dev()`, `variance()`.
#'
#' Every one splices the held values back into `theta` and delegates to the
#' parent. The derivative methods then subset the parent's answer by the names
#' the free parameter set generates, which is safe because the free set
#' preserves the parent's order: the same combination of free parameters
#' produces the same component string under either enumeration.
#'
#' @aliases distrib_atoms.FixedContinuousDistrib
#' @aliases distrib_cdf.FixedContinuousDistrib
#' @aliases distrib_cross2_y.FixedContinuousDistrib
#' @aliases distrib_cross_y.FixedContinuousDistrib
#' @aliases distrib_deriv3.FixedContinuousDistrib
#' @aliases distrib_deriv4.FixedContinuousDistrib
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
#' The S7 class of a CONTINUOUS distribution in which some parameters of the
#' wrapped distribution are held at known values. It is the only wrapper in
#' this package that REMOVES parameters: the law is the parent's, evaluated at
#' a `theta` short by however many were fixed. Build one with [fixed()], which
#' picks this class or one of its two siblings from the parent.
#'
#' @details
#' # How every method works
#'
#' The free parameters are the parent's minus the fixed ones, IN THE PARENT'S
#' ORDER. Every method splices the fixed values back into `theta` at their
#' positions and delegates to the parent, so the parent's closed forms are used
#' wherever they exist, and a derivative method then keeps only the components
#' in which every index is a free parameter. No method of this class computes
#' anything of its own, and none is faster or more accurate than the parent's.
#'
#' # Why the components can be subset by name
#'
#' The free set preserves the parent's ORDER, so a combination of free
#' parameters produces the same name string under the parent's enumeration as
#' under the wrapper's. Subsetting by name therefore cannot pair a component
#' with the wrong indices. That is the mistake re-parsing a name by splitting
#' on the underscore commits, for a parameter whose own name contains one.
#'
#' # The methods carry no documentation pages
#'
#' The continuous and discrete branches register theirs inside a loop over the
#' two classes, so there is no top-level assignment for roxygen to attach a
#' block to. The multivariate branch registers its own at the top level and is
#' left undocumented to match, every one of them being the same delegation.
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
#' # A gaussian centered at zero: one free parameter instead of two.
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
#' # Fixing everything is legal and gives a fully known law.
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
#' `fixed()` registers 22 methods on this class:
#' `distrib_atoms()`, `distrib_cdf()`, `distrib_cross2_y()`, `distrib_cross_y()`, `distrib_deriv3()`, `distrib_deriv4()`, `distrib_expected_hessian()`, `distrib_grad_cdf()`, `distrib_grad_y()`, `distrib_grad_y_hess()`, `distrib_gradient()`, `distrib_hess_cdf()`, `distrib_hess_y()`, `distrib_hess_y_hess()`, `distrib_hessian()`, `distrib_pdf()`, `distrib_quantile()`, `distrib_rng()`, `kurtosis()`, `skewness()`, `std_dev()`, `variance()`.
#'
#' Every one splices the held values back into `theta` and delegates to the
#' parent, exactly as on [FixedContinuousDistrib()]; the two classes share
#' one registration loop. The response derivatives are inherited refusals
#' from `discrete_distrib` and are not among them.
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
#' The S7 class of a DISCRETE distribution in which some parameters of the
#' wrapped distribution are held at known values. It behaves exactly as
#' [FixedContinuousDistrib] does; the split into three classes exists so that
#' the wrapper inherits the right base class, and with it the right defaults
#' for anything the parent does not register. Build one with [fixed()], which
#' picks the class from the parent.
#'
#' @details
#' # How every method works
#'
#' The free parameters are the parent's minus the fixed ones, IN THE PARENT'S
#' ORDER. Every method splices the fixed values back into `theta` at their
#' positions and delegates to the parent, so the parent's closed forms are used
#' wherever they exist, and a derivative method then keeps only the components
#' in which every index is a free parameter. No method of this class computes
#' anything of its own, and none is faster or more accurate than the parent's.
#'
#' # Why the components can be subset by name
#'
#' The free set preserves the parent's ORDER, so a combination of free
#' parameters produces the same name string under the parent's enumeration as
#' under the wrapper's. Subsetting by name therefore cannot pair a component
#' with the wrong indices. That is the mistake re-parsing a name by splitting
#' on the underscore commits, for a parameter whose own name contains one.
#'
#' # The methods carry no documentation pages
#'
#' The continuous and discrete branches register theirs inside a loop over the
#' two classes, so there is no top-level assignment for roxygen to attach a
#' block to. The multivariate branch registers its own at the top level and is
#' left undocumented to match, every one of them being the same delegation.
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
#' # A Poisson with its mean known: no free parameter at all.
#' d <- fixed(poisson_distrib(), mu = 3)
#' c(n_params = d@n_params, class = class(d)[1])
#' all.equal(distrib_pdf(d, 0:3, list()), dpois(0:3, 3))
#'
#' # A wrapper's OWN parameter can be held, which is what makes a
#' # zero-inflated model with a known inflation rate.
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
#' `fixed()` registers 20 methods on this class:
#' `distrib_cross2_y()`, `distrib_cross_y()`, `distrib_deriv3()`, `distrib_deriv4()`, `distrib_expected_hessian()`, `distrib_grad_y()`, `distrib_grad_y_hess()`, `distrib_gradient()`, `distrib_hess_y()`, `distrib_hess_y_hess()`, `distrib_hessian()`, `distrib_pdf()`, `distrib_rng()`, `mv_derived()`, `mv_location()`, `mv_marginal()`, `mv_reference_draw()`, `mv_sigma()`, `mv_support()`, `variance()`.
#'
#' Every one splices the held values back into `theta` and delegates to the
#' parent. The set differs from the univariate classes': there is no
#' distribution function, no quantile and no atom, and the `mv_*` accessors
#' take their place. The generics a multivariate family rejects are not
#' registered here at all, so the refusal is inherited from
#' `multivariate_distrib` and keeps its own message.
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
#' The S7 class of a MULTIVARIATE distribution in which some parameters of the
#' wrapped distribution are held at known values. It behaves exactly as
#' [FixedContinuousDistrib] does and adds the multivariate contract:
#' [mv_location()], [mv_sigma()], [mv_marginal()], [mv_support()],
#' [mv_reference_draw()] and [mv_derived()] all delegate to the parent at the
#' spliced `theta`, and [mv_derived()] reports the parent's quantities with the
#' Jacobian columns of the fixed parameters removed.
#'
#' The motivating case is a CENTERED PRIOR: holding the mean components of a
#' multivariate family at zero leaves the matrix parameter alone, and a matrix
#' parameter alone is what a random effect is distributed by.
#'
#' @details
#' # How every method works
#'
#' The free parameters are the parent's minus the fixed ones, IN THE PARENT'S
#' ORDER. Every method splices the fixed values back into `theta` at their
#' positions and delegates to the parent, so the parent's closed forms are used
#' wherever they exist, and a derivative method then keeps only the components
#' in which every index is a free parameter. No method of this class computes
#' anything of its own, and none is faster or more accurate than the parent's.
#'
#' # Why the components can be subset by name
#'
#' The free set preserves the parent's ORDER, so a combination of free
#' parameters produces the same name string under the parent's enumeration as
#' under the wrapper's. Subsetting by name therefore cannot pair a component
#' with the wrong indices. That is the mistake re-parsing a name by splitting
#' on the underscore commits, for a parameter whose own name contains one.
#'
#' # The methods carry no documentation pages
#'
#' The continuous and discrete branches register theirs inside a loop over the
#' two classes, so there is no top-level assignment for roxygen to attach a
#' block to. The multivariate branch registers its own at the top level and is
#' left undocumented to match, every one of them being the same delegation.
#' Read the parent's page for what a method computes and this page for what the
#' wrapper does to it.
#'
#' # What the refusals do
#'
#' The multivariate branch sits BESIDE the continuous and discrete ones rather
#' than under either, so the generics a multivariate family rejects by design,
#' the distribution function and the quantile, are inherited unregistered and
#' go on rejecting through this wrapper.
#'
#' @param parent_distrib The wrapped `multivariate_distrib` object.
#' @param fixed_params A named list of the fixed values, one single finite
#'   number each, strictly inside the corresponding parameter's open domain.
#' @param n_dim The number of coordinates, carried from the parent: fixing a
#'   parameter removes it from the parameter set and leaves the dimension of
#'   the response alone.
#' @inheritParams distrib
#'
#' @return An S7 object of class `FixedMultivariateDistrib`, inheriting from
#'   `multivariate_distrib` and from `distrib`. It carries `parent_distrib` and
#'   `fixed_params` beside the parent's properties, `n_dim` included.
#'
#' @seealso [fixed()] to build one, [FixedContinuousDistrib] and
#'   [FixedDiscreteDistrib] for the two siblings, [mvgaussian_distrib()] for a
#'   parent, and [mv_summary()] for the quantities a fit of one reports.
#'
#' @examples
#' # A centered two-dimensional gaussian: the matrix alone, which is what a
#' # random effect is distributed by.
#' d <- fixed(mvgaussian_distrib(2), mu1 = 0, mu2 = 0)
#' d@params
#' d
#'
#' theta <- list(sigma_log_L1 = 0, sigma_log_L2 = 0, sigma_L2.1 = 0.5)
#'
#' # The law is the parent's at the spliced theta.
#' y <- rbind(c(0, 0), c(1, -1))
#' all.equal(distrib_pdf(d, y, theta, log = TRUE),
#'           distrib_pdf(mvgaussian_distrib(2), y,
#'                       c(list(mu1 = 0, mu2 = 0), theta), log = TRUE))
#'
#' # The multivariate contract survives the wrapper.
#' mv_sigma(d, theta)
#' names(mv_derived(d, theta)$value)
#' mv_marginal(d, theta, 1)$distrib@n_dim
#'
#' # And so do the refusals the base class registers.
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
#' forms, which is how `fixed()` of a `fixed()` COLLAPSES into one wrapper. The
#' two describe the same law either way, and nesting would leave the inner
#' object's free set to be reconstructed at every call.
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
#'     fixed(mvgaussian_distrib(2), mu1 = 0)))
#'
#' # Which is why two calls collapse into one wrapper holding both values.
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
#' parameter list the PARENT expects, in the parent's own order. Every method
#' of every `Fixed*` class begins with this call and then delegates, so it is
#' the single point at which the two parameter sets are reconciled.
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
#' # Which is exactly what the parent is then called at.
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
                                                               approx = c("bartlett", "integrate", "mc", "opg"),
                                                               nsim = 10000, ...) {
    res <- distrib_expected_hessian(distrib@parent_distrib, y,
      fixed_full_theta(distrib, theta),
      scale = "parameter", approx = approx, nsim = nsim, ...
    )
    res[hess_names(distrib@params)]
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

S7::method(distrib_cross_y, FixedMultivariateDistrib) <-
  function(distrib, y, theta, scale = c("parameter", "link"), ...) {
    out <- distrib_cross_y(distrib@parent_distrib, y,
                           fixed_full_theta(distrib, theta),
                           scale = scale, ...)
    out[distrib@params]
  }

S7::method(distrib_cross2_y, FixedMultivariateDistrib) <-
  function(distrib, y, theta, scale = c("parameter", "link"), ...) {
    out <- distrib_cross2_y(distrib@parent_distrib, y,
                            fixed_full_theta(distrib, theta),
                            scale = scale, ...)
    out[distrib@params]
  }

S7::method(distrib_grad_y_hess, FixedMultivariateDistrib) <-
  function(distrib, y, theta, scale = c("parameter", "link"), ...) {
    out <- distrib_grad_y_hess(distrib@parent_distrib, y,
                               fixed_full_theta(distrib, theta),
                               scale = scale, ...)
    out[hess_names(distrib@params)]
  }

S7::method(distrib_hess_y_hess, FixedMultivariateDistrib) <-
  function(distrib, y, theta, scale = c("parameter", "link"), ...) {
    out <- distrib_hess_y_hess(distrib@parent_distrib, y,
                               fixed_full_theta(distrib, theta),
                               scale = scale, ...)
    out[hess_names(distrib@params)]
  }

S7::method(distrib_gradient, FixedMultivariateDistrib) <-
  function(distrib, y, theta, scale = c("parameter", "link"), ...) {
    out <- distrib_gradient(distrib@parent_distrib, y,
                            fixed_full_theta(distrib, theta),
                            scale = scale, ...)
    out[distrib@params]
  }

S7::method(distrib_hessian, FixedMultivariateDistrib) <-
  function(distrib, y, theta, scale = c("parameter", "link"), ...) {
    out <- distrib_hessian(distrib@parent_distrib, y,
                           fixed_full_theta(distrib, theta),
                           scale = scale, ...)
    out[hess_names(distrib@params)]
  }

S7::method(distrib_expected_hessian, FixedMultivariateDistrib) <-
  function(distrib, y, theta, scale = c("parameter", "link"),
           approx = c("bartlett", "integrate", "mc", "opg"),
           nsim = 10000, ...) {
    out <- distrib_expected_hessian(distrib@parent_distrib, y,
                                    fixed_full_theta(distrib, theta),
                                    scale = scale, approx = approx,
                                    nsim = nsim, ...)
    out[hess_names(distrib@params)]
  }

S7::method(distrib_deriv3, FixedMultivariateDistrib) <-
  function(distrib, y, theta, expected = FALSE,
           scale = c("parameter", "link"),
           approx = c("integrate", "bartlett", "mc", "opg"),
           nsim = 10000, ...) {
    out <- distrib_deriv3(distrib@parent_distrib, y,
                          fixed_full_theta(distrib, theta),
                          expected = expected, scale = scale,
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
                          expected = expected, scale = scale,
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
#' `theta` carries only the free parameters, every generic answers as the
#' parent does at the full vector, and a derivative is the parent's restricted
#' to the free indices: a subvector of the score, a submatrix of the Hessian,
#' sub-arrays at orders three and four. Nothing is recomputed, no normalizing
#' constant changes, and [fit_distrib()] estimates the free parameters with
#' standard errors and intervals for them alone.
#'
#' # What is accepted
#'
#' Fixed values are single finite numbers, strictly inside the OPEN domain of
#' their parameter. Fixing a parameter of a distribution that is already a
#' fixed-parameter wrapper collapses the two into one wrapper around the
#' original parent. Fixing a WRAPPER's own parameter is allowed and useful:
#' `fixed(zero_inflated(d), zi = 0.3)` is a zero-inflated model with a known
#' inflation rate. Fixing every parameter is allowed too and gives a fully
#' known distribution with an empty parameter set. Calling with no named value
#' is an error: the result would be the parent unchanged, and returning it
#' silently would hide a missing argument.
#'
#' # What it is for
#'
#' A prior. `fixed(gaussian1_distrib(), mu = 0)` is the ridge penalty with its
#' scale free, `fixed(laplace2_distrib(), mu = 0)` is the lasso, and
#' `fixed(mvgaussian_distrib(p), mu1 = 0, ...)` is what a random effect is
#' distributed by. `fixed(folded(gaussian1_distrib()), mu = 0)` is the
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
#' @param ... The fixed values, named after the parameters they fix, as in
#'   `fixed(gaussian1_distrib(), mu = 0)`. Each must be a single finite number
#'   strictly inside its parameter's domain, and each name must be a parameter
#'   of `distrib`. A name that is not, a value outside the domain, a value that
#'   is not a single number, and an empty `...` are each rejected with an error
#'   saying which condition failed.
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
#' # A gaussian with a known mean: only sigma remains.
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
#' # The lasso prior: a Laplace in its rate, centered at zero.
#' lasso <- fixed(laplace2_distrib(), mu = 0)
#' lasso@params
#' b <- c(-1, 0.5, 2)
#' all.equal(-distrib_pdf(lasso, b, list(lambda = 2), log = TRUE),
#'           2 * abs(b) - log(2 / 2))
#'
#' # A wrapper's own parameter can be held.
#' fixed(zero_inflated(poisson_distrib()), zi = 0.3)@params
#'
#' # Two calls collapse into one wrapper, and fixing everything is legal.
#' fixed(fixed(gaussian1_distrib(), mu = 0), sigma = 1)@fixed_params
#'
#' # Four refusals, each naming the condition that failed.
#' try(fixed(gaussian1_distrib(), nope = 1))
#' try(fixed(gaussian1_distrib(), sigma = -1))
#' try(fixed(gaussian1_distrib(), mu = c(0, 1)))
#' try(fixed(gaussian1_distrib()))
#'
#' @export
fixed <- function(distrib, ...) {
  if (!S7::S7_inherits(distrib, continuous_distrib) &&
    !S7::S7_inherits(distrib, discrete_distrib) &&
    !S7::S7_inherits(distrib, multivariate_distrib)) {
    stop(paste0("Input must inherit from 'discrete_distrib', ",
                "'continuous_distrib' or 'multivariate_distrib'."),
      call. = FALSE
    )
  }

  fix <- list(...)
  if (length(fix) == 0L) {
    stop(paste0(
      "fixed() needs at least one named value, as in fixed(d, mu = 0). With\n",
      "  none, the result would be the parent distribution unchanged; returning\n",
      "  it silently would hide a missing argument rather than report it."
    ), call. = FALSE)
  }
  nms <- names(fix)
  if (is.null(nms) || any(!nzchar(nms))) {
    stop("Every fixed value must be named after the parameter it fixes.",
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
        "'%s' is not a free parameter of '%s'. Free parameters: %s.",
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
      "'%s' is not a parameter of '%s'. Parameters: %s.",
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
