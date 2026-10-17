# S7 Class for Distributions With Fixed Parameters (Multivariate)

This is the S7 class of a MULTIVARIATE distribution in which some
parameters of the wrapped distribution are held at known values. It
follows
[FixedContinuousDistrib](https://statmodels7.github.io/distributions7/reference/FixedContinuousDistrib.md)
in splicing the held values into `theta` and delegating to the parent.
It has no distribution function, quantile or atoms, and adds the
multivariate generics:
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
a multivariate family at zero leaves only the matrix parameter, and a
random effect is distributed according to a law with this matrix as its
only parameter.

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
  response unchanged.

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

The following methods are registered on this class, and each is linked
to the page that documents it:

- the distribution:
  [`distrib_pdf()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.md)
  and
  [`distrib_rng()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.md);

- the derivatives in the parameters:
  [`distrib_gradient()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.md),
  [`distrib_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.md),
  [`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.md),
  [`distrib_deriv3()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.md)
  and
  [`distrib_deriv4()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.md);

- the derivatives in the response:
  [`distrib_grad_y()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.md),
  [`distrib_hess_y()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.md),
  [`distrib_deriv3_y()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3_y.md)
  and
  [`distrib_deriv4_y()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3_y.md);

- the mixed derivatives in the response and the parameters:
  [`distrib_cross_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross_y.md),
  [`distrib_cross2_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross2_y.md),
  [`distrib_cross3_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross3_y.md),
  [`distrib_grad_y_hess()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y_hess.md)
  and
  [`distrib_hess_y_hess()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y_hess.md);

- the moments:
  [`variance()`](https://statmodels7.github.io/distributions7/reference/variance.md);

- the multivariate summaries:
  [`mv_location()`](https://statmodels7.github.io/distributions7/reference/mv_location.md),
  [`mv_sigma()`](https://statmodels7.github.io/distributions7/reference/mv_sigma.md),
  [`mv_marginal()`](https://statmodels7.github.io/distributions7/reference/mv_marginal.md),
  [`mv_support()`](https://statmodels7.github.io/distributions7/reference/mv_support.md),
  [`mv_derived()`](https://statmodels7.github.io/distributions7/reference/mv_derived.md)
  and
  [`mv_reference_draw()`](https://statmodels7.github.io/distributions7/reference/mv_reference_draw.md);

- the exactness flags:
  [`expected_hessian_exact()`](https://statmodels7.github.io/distributions7/reference/expected_hessian_exact.md)
  and
  [`has_exact_deriv4()`](https://statmodels7.github.io/distributions7/reference/has_exact_deriv4.md).

Every other generic is inherited from
[`multivariate_distrib()`](https://statmodels7.github.io/distributions7/reference/multivariate_distrib.md).

## How every method works

The free parameters are the parent's minus the fixed ones, IN THE
PARENT'S ORDER. Every method splices the fixed values back into `theta`
at their positions and delegates to the parent, so the parent's closed
forms are used wherever they exist, and a derivative method then keeps
only the components in which every index is a free parameter. Every
method of this class delegates its computation to the parent, so the
speed and the accuracy of each method are those of the parent's method.

## Why the components can be subset by name

The free set preserves the parent's ORDER, so a combination of free
parameters produces the same name string under the parent's enumeration
as under the wrapper's. Subsetting by name therefore cannot pair a
component with the wrong indices. Re-parsing a name by splitting on the
underscore makes this mistake for a parameter whose own name contains an
underscore.

## The methods carry no documentation pages

The continuous and discrete branches register their methods inside a
loop over the two classes, so roxygen has no top-level assignment to
which it can attach a block. The multivariate branch registers its
methods at the top level and leaves them undocumented to match, because
each of them is the same delegation. Read the parent's page for what a
method computes and this page for what the wrapper does to it.

## Generics that are rejected

The multivariate branch sits BESIDE the continuous and discrete ones
rather than under either, so the generics that a multivariate family
rejects by design, the distribution function and the quantile, are
inherited unregistered and go on rejecting through this wrapper.

## See also

[`fixed()`](https://statmodels7.github.io/distributions7/reference/fixed.md)
to build one,
[FixedContinuousDistrib](https://statmodels7.github.io/distributions7/reference/FixedContinuousDistrib.md)
and
[FixedDiscreteDistrib](https://statmodels7.github.io/distributions7/reference/FixedDiscreteDistrib.md)
for the two siblings,
[`mvgaussian1_distrib()`](https://statmodels7.github.io/distributions7/reference/mvgaussian1_distrib.md)
for a parent, and
[`mv_summary()`](https://statmodels7.github.io/distributions7/reference/mv_summary.md)
for the quantities that a fit of one reports.

## Examples

``` r
# A centered two-dimensional gaussian has only the matrix parameter, as the
# law of a random effect does.
d <- fixed(mvgaussian1_distrib(2), mu1 = 0, mu2 = 0)
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
          distrib_pdf(mvgaussian1_distrib(2), y,
                      c(list(mu1 = 0, mu2 = 0), theta), log = TRUE))
#> [1] TRUE

# The multivariate generics work through the wrapper.
mv_sigma(d, theta)
#>     v1   v2
#> v1 1.0 0.50
#> v2 0.5 1.25
names(mv_derived(d, theta)$value)
#> [1] "sd_v1"     "sd_v2"     "cor_v1_v2"
mv_marginal(d, theta, 1)$distrib@n_dim
#> [1] 1

# And so do the errors that the base class registers.
try(distrib_cdf(d, y, theta))
#> Error : distrib_cdf() is not defined for 'fixed multivariate gaussian [2d, sigma=log_cholesky] [mu1=0,mu2=0]': the distribution function is an integral over an orthant, with no closed form and no one-dimensional fallback.
```
