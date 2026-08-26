# Cauchy Distribution Class

The S7 class of the Cauchy family with location \\\mu\\ and scale
\\\sigma \> 0\\, the Student's t on one degree of freedom, with density
\\f(y) = \[\pi\sigma(1 + ((y-\mu)/\sigma)^2)\]^{-1}\\ on the whole real
line. It inherits from `continuous_distrib`, so it answers every generic
of the `distrib` contract; the eleven methods listed below are
registered on it directly and everything else comes from the parent.

No moment of this family exists, so
[`mean()`](https://rdrr.io/r/base/mean.html),
[`variance()`](https://statmodels7.github.io/distributions7/reference/variance.md),
[`skewness()`](https://statmodels7.github.io/distributions7/reference/skewness.md)
and
[`kurtosis()`](https://statmodels7.github.io/distributions7/reference/kurtosis.md)
return `NaN`. Its location and scale are the median and the
half-interquartile range, both available from
[`distrib_quantile()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.md).

Build one with
[`cauchy_distrib()`](https://statmodels7.github.io/distributions7/reference/cauchy_distrib.md).
This page documents the raw S7 constructor, which takes the parent's
properties and validates none of the relationships between them.

## Usage

``` r
CauchyDistrib(
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

An S7 object of class `CauchyDistrib`, inheriting from
`continuous_distrib` and from `distrib`. Its properties are the
parent's: `distrib_name`, `dimension`, `bounds`, `params`,
`params_interpretation`, `n_params`, `params_bounds`, `link_params` and
`params_smooth`. For an object built by
[`cauchy_distrib()`](https://statmodels7.github.io/distributions7/reference/cauchy_distrib.md)
they hold `"cauchy"`, `"univariate"`, `c(-Inf, Inf)`,
`c("mu", "sigma")`, the interpretations
`c(mu = "location", sigma = "scale")`, `2`, the domains \\(-\infty,
\infty)\\ and \\(0, \infty)\\, and the two links.

## Methods

Registered on this class:
[`distrib_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.CauchyDistrib.md),
[`distrib_deriv3()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.CauchyDistrib.md),
[`distrib_deriv4()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.CauchyDistrib.md),
[`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.CauchyDistrib.md),
[`distrib_grad_y()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.CauchyDistrib.md),
[`distrib_gradient()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.CauchyDistrib.md),
[`distrib_hess_y()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.CauchyDistrib.md),
[`distrib_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.CauchyDistrib.md),
[`distrib_pdf()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.CauchyDistrib.md),
[`distrib_quantile()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.CauchyDistrib.md),
[`distrib_rng()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.CauchyDistrib.md)

Everything else is inherited from
[`continuous_distrib()`](https://statmodels7.github.io/distributions7/reference/continuous_distrib.md).

## See also

[`cauchy_distrib()`](https://statmodels7.github.io/distributions7/reference/cauchy_distrib.md)
to build one;
[`student_t1_distrib()`](https://statmodels7.github.io/distributions7/reference/student_t1_distrib.md),
of which this is the case \\\nu = 1\\;
[`skewness.CauchyDistrib()`](https://statmodels7.github.io/distributions7/reference/skewness.CauchyDistrib.md)
for why the moments are `NaN`;
[`distrib_gradient.CauchyDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.CauchyDistrib.md)
for the redescending score.

## Examples

``` r
d <- cauchy_distrib()
S7::S7_inherits(d, continuous_distrib)
#> [1] TRUE
d@params
#> [1] "mu"    "sigma"
d@params_interpretation
#>         mu      sigma 
#> "location"    "scale" 

# The interpretations say "location" and "scale", not "mean" and "standard
# deviation": no moment of this family exists.
th <- list(mu = 0.4, sigma = 1.5)
c(mean(d, th), variance(d, th))
#> [1] NaN NaN

# What does exist is the median and the half-interquartile range.
q <- distrib_quantile(d, c(0.25, 0.5, 0.75), th)
c(median = q[2], half_iqr = (q[3] - q[1]) / 2)
#>   median half_iqr 
#>      0.4      1.5 
```
