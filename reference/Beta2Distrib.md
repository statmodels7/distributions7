# Beta Distribution Class, the Two Shapes

The S7 class of the beta family on \\(0, 1)\\ in its canonical
parametrization, the two shapes \\\alpha \> 0\\ and \\\beta \> 0\\, with
density \\y^{\alpha-1}(1-y)^{\beta-1}/B(\alpha, \beta)\\. It inherits
from `continuous_distrib`, so it answers every generic of the `distrib`
contract; the eleven methods listed below are registered on it directly
and everything else comes from the parent.

Build one with
[`beta2_distrib()`](https://statmodels7.github.io/distributions7/reference/beta2_distrib.md),
which supplies the two link functions and fills the properties in. This
page documents the raw S7 constructor, which takes the parent's
properties and validates none of the relationships between them.

## Usage

``` r
Beta2Distrib(
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

An S7 object of class `Beta2Distrib`, inheriting from
`continuous_distrib` and from `distrib`. Its properties are the
parent's: `distrib_name`, `dimension`, `bounds`, `params`,
`params_interpretation`, `n_params`, `params_bounds`, `link_params` and
`params_smooth`. For an object built by
[`beta2_distrib()`](https://statmodels7.github.io/distributions7/reference/beta2_distrib.md)
they hold `"beta2"`, `"univariate"`, `c(0, 1)`, `c("alpha", "beta")`,
the interpretations `c(alpha = "shape", beta = "shape")`, `2`, the
domain \\(0, \infty)\\ for both parameters, and the two links.

## Methods

Registered on this class:
[`distrib_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.Beta2Distrib.md),
[`distrib_deriv3()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.Beta2Distrib.md),
[`distrib_deriv4()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.Beta2Distrib.md),
[`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.Beta2Distrib.md),
[`distrib_grad_y()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.Beta2Distrib.md),
[`distrib_hess_y()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.Beta2Distrib.md),
[`distrib_gradient()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.Beta2Distrib.md),
[`distrib_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.Beta2Distrib.md),
[`distrib_pdf()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.Beta2Distrib.md),
[`distrib_quantile()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.Beta2Distrib.md),
[`distrib_rng()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.Beta2Distrib.md)

Everything else is inherited from
[`continuous_distrib()`](https://statmodels7.github.io/distributions7/reference/continuous_distrib.md).

## See also

[`beta2_distrib()`](https://statmodels7.github.io/distributions7/reference/beta2_distrib.md)
to build one;
[`beta1_distrib()`](https://statmodels7.github.io/distributions7/reference/beta1_distrib.md)
for the same law in mean and precision;
[`distrib_pdf.Beta2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.Beta2Distrib.md)
and
[`distrib_gradient.Beta2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.Beta2Distrib.md)
for the closed forms this class supplies.

## Examples

``` r
d <- beta2_distrib()
S7::S7_inherits(d, continuous_distrib)
#> [1] TRUE

# The properties a consumer reads to drive the family without knowing it.
d@params
#> [1] "alpha" "beta" 
d@params_interpretation
#>   alpha    beta 
#> "shape" "shape" 
d@bounds
#> [1] 0 1

# Both parameters are shapes, so neither is a mean: the mean is their ratio.
mean(d, list(alpha = 2, beta = 5))
#> [1] 0.2857143
```
