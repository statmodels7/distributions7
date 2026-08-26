# Gumbel Distribution Class

The S7 class of the Gumbel (type I extreme value) family in the form for
**maxima**, a location-scale family on the whole real line with location
\\\mu\\ and scale \\\sigma \> 0\\ and density \\\sigma^{-1}\exp\\-z -
e^{-z}\\\\ at \\z = (y-\mu)/\sigma\\. Its shape is fixed: the skewness
and the excess kurtosis are constants and cannot be fitted. It inherits
from `continuous_distrib`, so it answers every generic of the `distrib`
contract; the eleven methods listed below are registered on it directly
and everything else comes from the parent.

Build one with
[`gumbel_distrib()`](https://statmodels7.github.io/distributions7/reference/gumbel_distrib.md),
which supplies the two link functions and fills the properties in. This
page documents the raw S7 constructor, which takes the parent's
properties and validates none of the relationships between them.

## Usage

``` r
GumbelDistrib(
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

An S7 object of class `GumbelDistrib`, inheriting from
`continuous_distrib` and from `distrib`. Its properties are the
parent's: `distrib_name`, `dimension`, `bounds`, `params`,
`params_interpretation`, `n_params`, `params_bounds`, `link_params` and
`params_smooth`. For an object built by
[`gumbel_distrib()`](https://statmodels7.github.io/distributions7/reference/gumbel_distrib.md)
they hold `"gumbel"`, `"univariate"`, `c(-Inf, Inf)`,
`c("mu", "sigma")`, `2`, the domains \\(-\infty, \infty)\\ and \\(0,
\infty)\\, and the two links.

## Methods

Registered on this class:
[`distrib_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.GumbelDistrib.md),
[`distrib_deriv3()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.GumbelDistrib.md),
[`distrib_deriv4()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.GumbelDistrib.md),
[`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.GumbelDistrib.md),
[`distrib_grad_y()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.GumbelDistrib.md),
[`distrib_gradient()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.GumbelDistrib.md),
[`distrib_hess_y()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.GumbelDistrib.md),
[`distrib_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.GumbelDistrib.md),
[`distrib_pdf()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.GumbelDistrib.md),
[`distrib_quantile()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.GumbelDistrib.md),
[`distrib_rng()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.GumbelDistrib.md)

Everything else is inherited from
[`continuous_distrib()`](https://statmodels7.github.io/distributions7/reference/continuous_distrib.md).

## See also

[`gumbel_distrib()`](https://statmodels7.github.io/distributions7/reference/gumbel_distrib.md)
to build one;
[`weibull1_distrib()`](https://statmodels7.github.io/distributions7/reference/weibull1_distrib.md),
which is this family under \\e^{-Y}\\;
[`gpd_distrib()`](https://statmodels7.github.io/distributions7/reference/gpd_distrib.md)
for the other limit law of extremes;
[`distrib_pdf.GumbelDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.GumbelDistrib.md)
and
[`distrib_gradient.GumbelDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.GumbelDistrib.md)
for the closed forms this class supplies.

## Examples

``` r
d <- gumbel_distrib()
S7::S7_inherits(d, continuous_distrib)
#> [1] TRUE

# The properties a consumer reads to drive the family without knowing it.
d@params
#> [1] "mu"    "sigma"
d@params_bounds
#> $mu
#> [1] -Inf  Inf
#> 
#> $sigma
#> [1]   0 Inf
#> 
d@bounds
#> [1] -Inf  Inf

# The shape is fixed: skewness and excess kurtosis do not move with the
# parameters.
rbind(a = c(skew = skewness(d, list(mu = 0, sigma = 1)),
            kurt = kurtosis(d, list(mu = 0, sigma = 1))),
      b = c(skewness(d, list(mu = 5, sigma = 3)),
            kurtosis(d, list(mu = 5, sigma = 3))))
#>       skew kurt
#> a 1.139547  2.4
#> b 1.139547  2.4
```
