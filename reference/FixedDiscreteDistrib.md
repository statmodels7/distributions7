# S7 Class for Distributions With Fixed Parameters (Discrete)

This is the S7 class of a DISCRETE distribution in which some parameters
of the wrapped distribution are held at known values. It behaves exactly
as
[FixedContinuousDistrib](https://statmodels7.github.io/distributions7/reference/FixedContinuousDistrib.md)
does; the split into three classes exists so that the wrapper inherits
the right base class, and with it the right defaults for anything that
the parent does not register. Build one with
[`fixed()`](https://statmodels7.github.io/distributions7/reference/fixed.md),
which picks the class from the parent.

## Usage

``` r
FixedDiscreteDistrib(
  distrib_name = character(0),
  dimension = character(0),
  bounds = integer(0),
  params = character(0),
  params_interpretation = character(0),
  n_params = integer(0),
  params_bounds = list(),
  link_params = list(),
  params_smooth = logical(0),
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

- parent_distrib:

  The wrapped `discrete_distrib` object.

- fixed_params:

  A named list of the fixed values, one single finite number each,
  strictly inside the corresponding parameter's open domain.

## Value

An S7 object of class `FixedDiscreteDistrib`, inheriting from
`discrete_distrib` and from `distrib`. It carries `parent_distrib` and
`fixed_params` beside the parent's properties.

## Methods

The following methods are registered on this class, and each is linked
to the page that documents it:

- the distribution:
  [`distrib_pdf()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.md),
  [`distrib_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.md),
  [`distrib_quantile()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.md),
  [`distrib_rng()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.md)
  and
  [`distrib_atoms()`](https://statmodels7.github.io/distributions7/reference/distrib_atoms.md);

- the derivatives in the parameters:
  [`distrib_gradient()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.md),
  [`distrib_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.md),
  [`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.md),
  [`distrib_deriv3()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.md),
  [`distrib_deriv4()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.md)
  and
  [`distrib_deriv5()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv5.wrappers.md);

- the derivatives of the expected information:
  [`distrib_dexpected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_dexpected_hessian.md)
  and
  [`distrib_d2expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_d2expected_hessian.md);

- the derivatives in the response:
  [`distrib_grad_y()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.md)
  and
  [`distrib_hess_y()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.md);

- the mixed derivatives in the response and the parameters:
  [`distrib_cross_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross_y.md),
  [`distrib_cross2_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross2_y.md),
  [`distrib_grad_y_hess()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y_hess.md)
  and
  [`distrib_hess_y_hess()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y_hess.md);

- the derivatives of the distribution function:
  [`distrib_grad_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_cdf.md)
  and
  [`distrib_hess_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_cdf.md);

- the moments:
  [`variance()`](https://statmodels7.github.io/distributions7/reference/variance.md),
  [`std_dev()`](https://statmodels7.github.io/distributions7/reference/std_dev.md),
  [`skewness()`](https://statmodels7.github.io/distributions7/reference/skewness.md)
  and
  [`kurtosis()`](https://statmodels7.github.io/distributions7/reference/kurtosis.md);

- the exactness flags:
  [`expected_hessian_exact()`](https://statmodels7.github.io/distributions7/reference/expected_hessian_exact.md)
  and
  [`has_exact_deriv4()`](https://statmodels7.github.io/distributions7/reference/has_exact_deriv4.md);

- the fitting helpers:
  [`distrib_scalar_route()`](https://statmodels7.github.io/distributions7/reference/distrib_scalar_route.md).

Every other generic is inherited from
[`discrete_distrib()`](https://statmodels7.github.io/distributions7/reference/discrete_distrib.md).

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

## See also

[`fixed()`](https://statmodels7.github.io/distributions7/reference/fixed.md)
to build one,
[FixedContinuousDistrib](https://statmodels7.github.io/distributions7/reference/FixedContinuousDistrib.md)
and
[FixedMultivariateDistrib](https://statmodels7.github.io/distributions7/reference/FixedMultivariateDistrib.md)
for the two siblings, and
[`zero_inflated()`](https://statmodels7.github.io/distributions7/reference/zero_inflated.md),
one of whose own parameters can usefully be held.

## Examples

``` r
# A Poisson with its mean known has no free parameter at all.
d <- fixed(poisson_distrib(), mu = 3)
c(n_params = d@n_params, class = class(d)[1])
#>                               n_params                                  class 
#>                                    "0" "distributions7::FixedDiscreteDistrib" 
all.equal(distrib_pdf(d, 0:3, list()), dpois(0:3, 3))
#> [1] TRUE

# A wrapper's OWN parameter can be held, so a zero-inflated model can have
# a known inflation rate.
zi <- fixed(zero_inflated(poisson_distrib()), zi = 0.3)
zi@params
#> [1] "mu"
all.equal(distrib_pdf(zi, 0:3, list(mu = 3)),
          0.3 * (0:3 == 0) + 0.7 * dpois(0:3, 3))
#> [1] TRUE
```
