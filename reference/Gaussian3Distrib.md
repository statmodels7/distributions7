# Gaussian Distribution Class, Mean and Precision

The S7 class of the Gaussian (normal) family parametrized by its mean
\\\mu\\ and its precision \\\tau = 1/\sigma^2 \> 0\\, with density
\\f(y) = (\tau/2\pi)^{1/2}\exp\\-\tau(y-\mu)^2/2\\\\ on the whole real
line. It inherits from `continuous_distrib`, so it answers every generic
of the `distrib` contract; the eleven methods listed below are
registered on it directly and everything else comes from the parent.

Build one with
[`gaussian3_distrib()`](https://statmodels7.github.io/distributions7/reference/gaussian3_distrib.md),
which supplies the two link functions and fills the properties in. This
page documents the raw S7 constructor, which takes the parent's
properties and validates none of the relationships between them.

## Usage

``` r
Gaussian3Distrib(
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

An S7 object of class `Gaussian3Distrib`, inheriting from
`continuous_distrib` and from `distrib`. Its properties are the
parent's: `distrib_name`, `dimension`, `bounds`, `params`,
`params_interpretation`, `n_params`, `params_bounds`, `link_params` and
`params_smooth`. For an object built by
[`gaussian3_distrib()`](https://statmodels7.github.io/distributions7/reference/gaussian3_distrib.md)
they hold `"gaussian3"`, `"univariate"`, `c(-Inf, Inf)`,
`c("mu", "tau")`, the interpretations
`c(mu = "mean", tau = "precision")`, `2`, the domains \\(-\infty,
\infty)\\ and \\(0, \infty)\\, and the two links.

## Methods

Registered on this class:
[`distrib_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.Gaussian3Distrib.md),
[`distrib_deriv3()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.Gaussian3Distrib.md),
[`distrib_deriv4()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.Gaussian3Distrib.md),
[`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.Gaussian3Distrib.md),
[`distrib_grad_y()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.Gaussian3Distrib.md),
[`distrib_gradient()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.Gaussian3Distrib.md),
[`distrib_hess_y()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.Gaussian3Distrib.md),
[`distrib_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.Gaussian3Distrib.md),
[`distrib_pdf()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.Gaussian3Distrib.md),
[`distrib_quantile()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.Gaussian3Distrib.md),
[`distrib_rng()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.Gaussian3Distrib.md)

Everything else is inherited from
[`continuous_distrib()`](https://statmodels7.github.io/distributions7/reference/continuous_distrib.md).

## See also

[`gaussian3_distrib()`](https://statmodels7.github.io/distributions7/reference/gaussian3_distrib.md)
to build one;
[`gaussian1_distrib()`](https://statmodels7.github.io/distributions7/reference/gaussian1_distrib.md)
for the same law in mean and standard deviation and
[`gaussian2_distrib()`](https://statmodels7.github.io/distributions7/reference/gaussian2_distrib.md)
for mean and variance;
[`distrib_pdf.Gaussian3Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.Gaussian3Distrib.md)
and
[`distrib_gradient.Gaussian3Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.Gaussian3Distrib.md)
for the closed forms this class supplies.

## Examples

``` r
d <- gaussian3_distrib()
S7::S7_inherits(d, continuous_distrib)
#> [1] TRUE

# The properties a consumer reads to drive the family without knowing it.
d@params
#> [1] "mu"  "tau"
d@params_interpretation
#>          mu         tau 
#>      "mean" "precision" 
d@params_bounds
#> $mu
#> [1] -Inf  Inf
#> 
#> $tau
#> [1]   0 Inf
#> 

# The precision is the reciprocal of the variance, so a large tau is a
# tight distribution.
variance(d, list(mu = 0, tau = 100))
#> [1] 0.01
```
