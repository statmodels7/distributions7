# Weibull Distribution Class, Scale and Shape

The S7 class of the Weibull family parametrized by a scale \\\mu \> 0\\
and a shape \\\sigma \> 0\\, with density \\f(y) =
(\sigma/\mu)(y/\mu)^{\sigma-1}\exp\\-(y/\mu)^{\sigma}\\\\ on the
positive half line. It inherits from `continuous_distrib`, so it answers
every generic of the `distrib` contract; the methods listed below are
registered on it directly and everything else comes from the parent.

Build one with
[`weibull1_distrib()`](https://statmodels7.github.io/distributions7/reference/weibull1_distrib.md),
which supplies the two link functions and fills the properties in. This
page documents the raw S7 constructor, which takes the parent's
properties and validates none of the relationships between them.

## Usage

``` r
Weibull1Distrib(
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

An S7 object of class `Weibull1Distrib`, inheriting from
`continuous_distrib` and from `distrib`. Its properties are the
parent's: `distrib_name`, `dimension`, `bounds`, `params`,
`params_interpretation`, `n_params`, `params_bounds`, `link_params` and
`params_smooth`. For an object built by
[`weibull1_distrib()`](https://statmodels7.github.io/distributions7/reference/weibull1_distrib.md)
they hold `"weibull1"`, `"univariate"`, `c(0, Inf)`, `c("mu", "sigma")`,
the interpretations `c(mu = "scale", sigma = "shape")`, `2`, the domain
\\(0, \infty)\\ for both parameters, and the two links.

## Methods

Registered on this class in this file:
[`distrib_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.Weibull1Distrib.md),
[`distrib_deriv3()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.Weibull1Distrib.md),
[`distrib_deriv4()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.Weibull1Distrib.md),
[`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.Weibull1Distrib.md),
[`distrib_grad_y()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.Weibull1Distrib.md),
[`distrib_gradient()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.Weibull1Distrib.md),
[`distrib_hess_y()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.Weibull1Distrib.md),
[`distrib_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.Weibull1Distrib.md),
[`distrib_pdf()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.Weibull1Distrib.md),
[`distrib_quantile()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.Weibull1Distrib.md),
[`distrib_rng()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.Weibull1Distrib.md).

Four more are registered on the class from other files:
[`distrib_cross_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross_y.Weibull1Distrib.md)
in `cross_derivatives_families.R`,
[`distrib_grad_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_cdf.Weibull1Distrib.md)
in `cdf_survival_higher.R`, and the four moments
[`mean()`](https://statmodels7.github.io/distributions7/reference/mean.Weibull1Distrib.md),
[`variance()`](https://statmodels7.github.io/distributions7/reference/variance.Weibull1Distrib.md),
[`skewness()`](https://statmodels7.github.io/distributions7/reference/skewness.Weibull1Distrib.md)
and
[`kurtosis()`](https://statmodels7.github.io/distributions7/reference/kurtosis.Weibull1Distrib.md)
in `moments.R`.

Everything else is inherited from
[`continuous_distrib()`](https://statmodels7.github.io/distributions7/reference/continuous_distrib.md).

## See also

[`weibull1_distrib()`](https://statmodels7.github.io/distributions7/reference/weibull1_distrib.md)
to build one;
[`weibull3_distrib()`](https://statmodels7.github.io/distributions7/reference/weibull3_distrib.md)
for the same law parametrized by a quantile;
[`gumbel_distrib()`](https://statmodels7.github.io/distributions7/reference/gumbel_distrib.md),
which is this family on a negative log scale;
[`distrib_pdf.Weibull1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.Weibull1Distrib.md)
and
[`distrib_gradient.Weibull1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.Weibull1Distrib.md)
for the closed forms it supplies.

## Examples

``` r
d <- weibull1_distrib()
S7::S7_inherits(d, continuous_distrib)
#> [1] TRUE

# The properties a consumer reads to drive the family without knowing it.
d@params
#> [1] "mu"    "sigma"
d@params_interpretation
#>      mu   sigma 
#> "scale" "shape" 
d@bounds
#> [1]   0 Inf

# Both parameters are positive, so both ride a log by default.
vapply(d@link_params, function(l) l@link_name, character(1))
#>    mu sigma 
#> "log" "log" 
```
