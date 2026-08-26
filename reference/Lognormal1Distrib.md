# Lognormal Distribution Class

The S7 class of the lognormal family on \\(0, \infty)\\: the law of a
variable whose logarithm is Gaussian with mean \\\mu\\ and variance
\\\sigma^2 \> 0\\. Both parameters live **on the log scale**, so neither
is the mean or the variance of \\Y\\ itself. It inherits from
`continuous_distrib`, so it answers every generic of the `distrib`
contract; the eleven methods listed below are registered on it directly
and everything else comes from the parent.

Build one with
[`lognormal1_distrib()`](https://statmodels7.github.io/distributions7/reference/lognormal1_distrib.md),
which supplies the two link functions and fills the properties in. This
page documents the raw S7 constructor, which takes the parent's
properties and validates none of the relationships between them.

## Usage

``` r
Lognormal1Distrib(
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

An S7 object of class `Lognormal1Distrib`, inheriting from
`continuous_distrib` and from `distrib`. Its properties are the
parent's: `distrib_name`, `dimension`, `bounds`, `params`,
`params_interpretation`, `n_params`, `params_bounds`, `link_params` and
`params_smooth`. For an object built by
[`lognormal1_distrib()`](https://statmodels7.github.io/distributions7/reference/lognormal1_distrib.md)
they hold `"lognormal1"`, `"univariate"`, `c(0, Inf)`,
`c("mu", "sigma2")`, `2`, the domains \\(-\infty, \infty)\\ and \\(0,
\infty)\\, and the two links.

## Methods

Registered on this class:
[`distrib_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.Lognormal1Distrib.md),
[`distrib_deriv3()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.Lognormal1Distrib.md),
[`distrib_deriv4()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.Lognormal1Distrib.md),
[`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.Lognormal1Distrib.md),
[`distrib_grad_y()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.Lognormal1Distrib.md),
[`distrib_gradient()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.Lognormal1Distrib.md),
[`distrib_hess_y()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.Lognormal1Distrib.md),
[`distrib_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.Lognormal1Distrib.md),
[`distrib_pdf()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.Lognormal1Distrib.md),
[`distrib_quantile()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.Lognormal1Distrib.md),
[`distrib_rng()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.Lognormal1Distrib.md)

Everything else is inherited from
[`continuous_distrib()`](https://statmodels7.github.io/distributions7/reference/continuous_distrib.md).

## See also

[`lognormal1_distrib()`](https://statmodels7.github.io/distributions7/reference/lognormal1_distrib.md)
to build one;
[`lognormal2_distrib()`](https://statmodels7.github.io/distributions7/reference/lognormal2_distrib.md)
for the same law parametrized by the mean and the variance of \\Y\\;
[`gaussian2_distrib()`](https://statmodels7.github.io/distributions7/reference/gaussian2_distrib.md),
which this becomes at \\\log y\\;
[`distrib_pdf.Lognormal1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.Lognormal1Distrib.md)
and
[`distrib_gradient.Lognormal1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.Lognormal1Distrib.md)
for the closed forms this class supplies.

## Examples

``` r
d <- lognormal1_distrib()
S7::S7_inherits(d, continuous_distrib)
#> [1] TRUE

# The properties a consumer reads to drive the family without knowing it.
d@params
#> [1] "mu"     "sigma2"
d@params_bounds
#> $mu
#> [1] -Inf  Inf
#> 
#> $sigma2
#> [1]   0 Inf
#> 
d@bounds
#> [1]   0 Inf

# mu is the mean of log(Y), so exp(mu) is the MEDIAN of Y and the mean sits
# above it.
th <- list(mu = 0.5, sigma2 = 0.36)
c(exp_mu = exp(0.5), median = distrib_quantile(d, 0.5, th),
  mean = mean(d, th))
#>   exp_mu   median     mean 
#> 1.648721 1.648721 1.973878 
```
