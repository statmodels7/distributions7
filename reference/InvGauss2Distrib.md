# Inverse Gaussian Distribution Class, Mean and Shape

The S7 class of the inverse Gaussian family on \\(0, \infty)\\ in its
classical parametrization, the mean \\\mu \> 0\\ and the shape \\\lambda
\> 0\\, so that \\\operatorname{Var}(Y) = \mu^3/\lambda\\. It inherits
from `continuous_distrib`, so it answers every generic of the `distrib`
contract; the nine methods listed below are registered on it in this
file and everything else comes from the parent.

Build one with
[`invgauss2_distrib()`](https://statmodels7.github.io/distributions7/reference/invgauss2_distrib.md),
which supplies the two link functions and fills the properties in. This
page documents the raw S7 constructor, which takes the parent's
properties and validates none of the relationships between them.

## Usage

``` r
InvGauss2Distrib(
  distrib_name = character(0),
  dimension = character(0),
  bounds = integer(0),
  params = character(0),
  params_interpretation = character(0),
  n_params = integer(0),
  params_bounds = list(),
  link_params = list(),
  params_smooth = logical(0)
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

## Value

An S7 object of class `InvGauss2Distrib`, inheriting from
`continuous_distrib` and from `distrib`. Its properties are the
parent's: `distrib_name`, `dimension`, `bounds`, `params`,
`params_interpretation`, `n_params`, `params_bounds`, `link_params` and
`params_smooth`. For an object built by
[`invgauss2_distrib()`](https://statmodels7.github.io/distributions7/reference/invgauss2_distrib.md)
they hold `"invgauss2"`, `"univariate"`, `c(0, Inf)`,
`c("mu", "lambda")`, the interpretations
`c(mu = "mean", lambda = "shape")`, `2`, the domain \\(0, \infty)\\ for
both parameters, and the two links.

## Methods

Registered on this class in this file:
[`distrib_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.InvGauss2Distrib.md),
[`distrib_deriv3()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.InvGauss2Distrib.md),
[`distrib_deriv4()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.InvGauss2Distrib.md),
[`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.InvGauss2Distrib.md),
[`distrib_gradient()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.InvGauss2Distrib.md),
[`distrib_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.InvGauss2Distrib.md),
[`distrib_pdf()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.InvGauss2Distrib.md),
[`distrib_quantile()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.InvGauss2Distrib.md),
[`distrib_rng()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.InvGauss2Distrib.md)

Three more are registered elsewhere in the package and are closed form
too:
[`distrib_grad_y()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.InvGauss2Distrib.md),
[`distrib_hess_y()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.InvGauss2Distrib.md)
and
[`distrib_cross_y()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.InvGauss2Distrib.md),
the derivatives in the response and the mixed block, and
[`distrib_grad_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_cdf.InvGauss2Distrib.md)
for the derivatives of the distribution function.

Everything else is inherited from
[`continuous_distrib()`](https://statmodels7.github.io/distributions7/reference/continuous_distrib.md).

## See also

[`invgauss2_distrib()`](https://statmodels7.github.io/distributions7/reference/invgauss2_distrib.md)
to build one;
[`invgauss1_distrib()`](https://statmodels7.github.io/distributions7/reference/invgauss1_distrib.md)
for the same law in mean and dispersion \\\phi = 1/\lambda\\;
[`distrib_pdf.InvGauss2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.InvGauss2Distrib.md)
and
[`distrib_gradient.InvGauss2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.InvGauss2Distrib.md)
for the closed forms this class supplies.

## Examples

``` r
d <- invgauss2_distrib()
S7::S7_inherits(d, continuous_distrib)
#> [1] TRUE

# The properties a consumer reads to drive the family without knowing it.
d@params
#> [1] "mu"     "lambda"
d@params_interpretation
#>      mu  lambda 
#>  "mean" "shape" 
d@bounds
#> [1]   0 Inf

# lambda is a shape, so a larger value is a tighter law at a fixed mean.
vapply(c(1, 3, 30), function(l) variance(d, list(mu = 2, lambda = l)),
       numeric(1))
#> [1] 8.0000000 2.6666667 0.2666667
```
