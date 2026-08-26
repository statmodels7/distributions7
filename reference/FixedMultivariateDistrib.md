# S7 Class for Distributions With Fixed Parameters (Multivariate)

The S7 class of a MULTIVARIATE distribution in which some parameters of
the wrapped distribution are held at known values. It behaves exactly as
[FixedContinuousDistrib](https://statmodels7.github.io/distributions7/reference/FixedContinuousDistrib.md)
does and adds the multivariate contract:
[`mv_location()`](https://statmodels7.github.io/distributions7/reference/mv_location.md),
[`mv_sigma()`](https://statmodels7.github.io/distributions7/reference/mv_sigma.md),
[`mv_marginal()`](https://statmodels7.github.io/distributions7/reference/mv_marginal.md),
[`mv_support()`](https://statmodels7.github.io/distributions7/reference/mv_support.md),
[`mv_reference_draw()`](https://statmodels7.github.io/distributions7/reference/mv_reference_draw.md)
and
[`mv_derived()`](https://statmodels7.github.io/distributions7/reference/mv_derived.md)
all delegate to the parent at the spliced `theta`, and
[`mv_derived()`](https://statmodels7.github.io/distributions7/reference/mv_derived.md)
reports the parent's quantities with the Jacobian columns of the fixed
parameters removed.

The motivating case is a CENTERED PRIOR: holding the mean components of
a multivariate family at zero leaves the matrix parameter alone, and a
matrix parameter alone is what a random effect is distributed by.

## Usage

``` r
FixedMultivariateDistrib(
  distrib_name = character(0),
  dimension = character(0),
  bounds = integer(0),
  params = character(0),
  params_interpretation = character(0),
  n_params = integer(0),
  params_bounds = list(),
  link_params = list(),
  params_smooth = logical(0),
  n_dim = integer(0),
  parent_distrib = distrib(),
  fixed_params = list()
)
```

## Arguments

- distrib_name:

  A single character string specifying the name of the distribution
  (e.g., `"student t"`).

- dimension:

  A character string indicating the dimensionality (`"univariate"` or
  `"multivariate"`).

- bounds:

  A numeric vector of length 2 defining the overall support of the
  distribution `c(lower, upper)`.

- params:

  A character vector containing the names of the distribution parameters
  (e.g., `c("mu", "sigma")`).

- params_interpretation:

  A character vector (typically named) providing the statistical
  interpretation of each parameter (e.g., `c(mu = "location")`).

- n_params:

  A numeric value specifying the total number of parameters.

- params_bounds:

  A list of numeric vectors of length 2, specifying the valid
  mathematical domain for each individual parameter.

- link_params:

  A list of link function objects corresponding to each parameter,
  primarily used to map parameters to the unconstrained real line for
  optimization algorithms.

- params_smooth:

  An optional named logical vector flagging, for each parameter, whether
  the log-likelihood is differentiable with respect to it. Defaults to
  all `TRUE` (leave empty). Set an entry to `FALSE` for parameters at
  which the log-likelihood has a kink (e.g. the location of a Laplace
  distribution): the observed Hessian is then degenerate and the
  expected information must be obtained from the score variance rather
  than from \\-\mathbb{E}\[H\]\\ (see
  [`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.md)).

- n_dim:

  The number of coordinates, carried from the parent: fixing a parameter
  removes it from the parameter set and leaves the dimension of the
  response alone.

- parent_distrib:

  The wrapped `multivariate_distrib` object.

- fixed_params:

  A named list of the fixed values, one single finite number each,
  strictly inside the corresponding parameter's open domain.

## Value

An S7 object of class `FixedMultivariateDistrib`, inheriting from
`multivariate_distrib` and from `distrib`. It carries `parent_distrib`
and `fixed_params` beside the parent's properties, `n_dim` included.

## Methods

[`fixed()`](https://statmodels7.github.io/distributions7/reference/fixed.md)
registers 20 methods on this class:
[`distrib_cross2_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross2_y.md),
[`distrib_cross_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross_y.md),
[`distrib_deriv3()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.md),
[`distrib_deriv4()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.md),
[`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.md),
[`distrib_grad_y()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.md),
[`distrib_grad_y_hess()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y_hess.md),
[`distrib_gradient()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.md),
[`distrib_hess_y()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.md),
[`distrib_hess_y_hess()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y_hess.md),
[`distrib_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.md),
[`distrib_pdf()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.md),
[`distrib_rng()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.md),
[`mv_derived()`](https://statmodels7.github.io/distributions7/reference/mv_derived.md),
[`mv_location()`](https://statmodels7.github.io/distributions7/reference/mv_location.md),
[`mv_marginal()`](https://statmodels7.github.io/distributions7/reference/mv_marginal.md),
[`mv_reference_draw()`](https://statmodels7.github.io/distributions7/reference/mv_reference_draw.md),
[`mv_sigma()`](https://statmodels7.github.io/distributions7/reference/mv_sigma.md),
[`mv_support()`](https://statmodels7.github.io/distributions7/reference/mv_support.md),
[`variance()`](https://statmodels7.github.io/distributions7/reference/variance.md).

Every one splices the held values back into `theta` and delegates to the
parent. The set differs from the univariate classes': there is no
distribution function, no quantile and no atom, and the `mv_*` accessors
take their place. The generics a multivariate family rejects are not
registered here at all, so the refusal is inherited from
`multivariate_distrib` and keeps its own message.

## How every method works

The free parameters are the parent's minus the fixed ones, IN THE
PARENT'S ORDER. Every method splices the fixed values back into `theta`
at their positions and delegates to the parent, so the parent's closed
forms are used wherever they exist, and a derivative method then keeps
only the components in which every index is a free parameter. No method
of this class computes anything of its own, and none is faster or more
accurate than the parent's.

## Why the components can be subset by name

The free set preserves the parent's ORDER, so a combination of free
parameters produces the same name string under the parent's enumeration
as under the wrapper's. Subsetting by name therefore cannot pair a
component with the wrong indices. That is the mistake re-parsing a name
by splitting on the underscore commits, for a parameter whose own name
contains one.

## The methods carry no documentation pages

The continuous and discrete branches register theirs inside a loop over
the two classes, so there is no top-level assignment for roxygen to
attach a block to. The multivariate branch registers its own at the top
level and is left undocumented to match, every one of them being the
same delegation. Read the parent's page for what a method computes and
this page for what the wrapper does to it.

## What the refusals do

The multivariate branch sits BESIDE the continuous and discrete ones
rather than under either, so the generics a multivariate family rejects
by design, the distribution function and the quantile, are inherited
unregistered and go on rejecting through this wrapper.

## See also

[`fixed()`](https://statmodels7.github.io/distributions7/reference/fixed.md)
to build one,
[FixedContinuousDistrib](https://statmodels7.github.io/distributions7/reference/FixedContinuousDistrib.md)
and
[FixedDiscreteDistrib](https://statmodels7.github.io/distributions7/reference/FixedDiscreteDistrib.md)
for the two siblings,
[`mvgaussian_distrib()`](https://statmodels7.github.io/distributions7/reference/mvgaussian_distrib.md)
for a parent, and
[`mv_summary()`](https://statmodels7.github.io/distributions7/reference/mv_summary.md)
for the quantities a fit of one reports.

## Examples

``` r
# A centered two-dimensional gaussian: the matrix alone, which is what a
# random effect is distributed by.
d <- fixed(mvgaussian_distrib(2), mu1 = 0, mu2 = 0)
d@params
#> [1] "sigma_log_L1" "sigma_log_L2" "sigma_L2.1"  
d
#> Distribution: Fixed Multivariate Gaussian [2d, Sigma=log_cholesky] [mu1=0,mu2=0]
#> Type:         Continuous, 2-dimensional
#> Dimensions:   multivariate
#> 
#> Parameters:
#>   sigma_log_L1 (covariance)         | Link: identity   | Domain: (-Inf, Inf)
#>   sigma_log_L2 (covariance)         | Link: identity   | Domain: (-Inf, Inf)
#>   sigma_L2.1   (covariance)         | Link: identity   | Domain: (-Inf, Inf)
#> 
#> Fixed:
#>   mu1 = 0
#>   mu2 = 0

theta <- list(sigma_log_L1 = 0, sigma_log_L2 = 0, sigma_L2.1 = 0.5)

# The law is the parent's at the spliced theta.
y <- rbind(c(0, 0), c(1, -1))
all.equal(distrib_pdf(d, y, theta, log = TRUE),
          distrib_pdf(mvgaussian_distrib(2), y,
                      c(list(mu1 = 0, mu2 = 0), theta), log = TRUE))
#> [1] TRUE

# The multivariate contract survives the wrapper.
mv_sigma(d, theta)
#>     v1   v2
#> v1 1.0 0.50
#> v2 0.5 1.25
names(mv_derived(d, theta)$value)
#> [1] "sd_v1"     "sd_v2"     "cor_v1_v2"
mv_marginal(d, theta, 1)$distrib@n_dim
#> [1] 1

# And so do the refusals the base class registers.
try(distrib_cdf(d, y, theta))
#> Error : distrib_cdf() is not defined for 'fixed multivariate gaussian [2d, sigma=log_cholesky] [mu1=0,mu2=0]': the distribution function is an integral over an orthant, with no closed form and no one-dimensional fallback.
```
