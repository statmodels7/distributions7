# Negative Binomial Distribution Class, NB2

The S7 class of the negative binomial family on the non-negative
integers in the NB2 parametrization: the mean \\\mu \> 0\\ and a
dispersion \\\theta \> 0\\, so that \\\operatorname{Var}(Y) = \mu +
\mu^2/\theta\\ and the variance is quadratic in the mean. It inherits
from `discrete_distrib`, so it answers every generic of the `distrib`
contract; the nine methods listed below are registered on it directly
and everything else comes from the parent.

Build one with
[`negbin2_distrib()`](https://statmodels7.github.io/distributions7/reference/negbin2_distrib.md),
which supplies the two link functions and fills the properties in. This
page documents the raw S7 constructor, which takes the parent's
properties and validates none of the relationships between them.

## Usage

``` r
NegBin2Distrib(
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

An S7 object of class `NegBin2Distrib`, inheriting from
`discrete_distrib` and from `distrib`. Its properties are the parent's:
`distrib_name`, `dimension`, `bounds`, `params`,
`params_interpretation`, `n_params`, `params_bounds`, `link_params` and
`params_smooth`. For an object built by
[`negbin2_distrib()`](https://statmodels7.github.io/distributions7/reference/negbin2_distrib.md)
they hold `"negbin2"`, `"univariate"`, `c(0, Inf)`, `c("mu", "theta")`,
the interpretations `c(mu = "mean", theta = "dispersion")`, `2`, the
domain \\(0, \infty)\\ for both parameters, and the two links.

## Methods

Registered on this class:
[`distrib_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.NegBin2Distrib.md),
[`distrib_deriv3()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.NegBin2Distrib.md),
[`distrib_deriv4()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.NegBin2Distrib.md),
[`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.NegBin2Distrib.md),
[`distrib_gradient()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.NegBin2Distrib.md),
[`distrib_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.NegBin2Distrib.md),
[`distrib_pdf()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.NegBin2Distrib.md),
[`distrib_quantile()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.NegBin2Distrib.md),
[`distrib_rng()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.NegBin2Distrib.md)

A discrete family has no derivatives in the response, so there are nine
here where a continuous family has eleven. Everything else is inherited
from
[`discrete_distrib()`](https://statmodels7.github.io/distributions7/reference/discrete_distrib.md).

## See also

[`negbin2_distrib()`](https://statmodels7.github.io/distributions7/reference/negbin2_distrib.md)
to build one;
[`negbin1_distrib()`](https://statmodels7.github.io/distributions7/reference/negbin1_distrib.md)
for the NB1 parametrization, whose variance is linear in the mean;
[`poisson_distrib()`](https://statmodels7.github.io/distributions7/reference/poisson_distrib.md)
for the limit as \\\theta\\ grows;
[`distrib_pdf.NegBin2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.NegBin2Distrib.md)
and
[`distrib_gradient.NegBin2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.NegBin2Distrib.md)
for the closed forms this class supplies.

## Examples

``` r
d <- negbin2_distrib()
S7::S7_inherits(d, discrete_distrib)
#> [1] TRUE

# The properties a consumer reads to drive the family without knowing it.
d@params
#> [1] "mu"    "theta"
d@params_interpretation
#>           mu        theta 
#>       "mean" "dispersion" 
d@bounds
#> [1]   0 Inf

# theta is a dispersion read the other way round from a variance: the
# smaller it is, the more overdispersed the counts.
vapply(c(0.5, 2, 1e6), function(t) variance(d, list(mu = 4, theta = t)),
       numeric(1))
#> [1] 36.000000 12.000000  4.000016
```
