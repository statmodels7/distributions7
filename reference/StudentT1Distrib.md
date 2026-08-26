# Student t Distribution Class, Location, Scale and Degrees of Freedom

The S7 class of the location-scale Student t family, parametrized by a
location \\\mu\\, a scale \\\sigma \> 0\\ and degrees of freedom \\\nu
\> 0\\, with density proportional to \\\\1 +
(y-\mu)^2/(\nu\sigma^2)\\^{-(\nu+1)/2}\\ on the whole real line. It
inherits from `continuous_distrib`, so it answers every generic of the
`distrib` contract; the eleven methods listed below are registered on it
directly and everything else comes from the parent.

Build one with
[`student_t1_distrib()`](https://statmodels7.github.io/distributions7/reference/student_t1_distrib.md),
which supplies the three link functions and fills the properties in.
This page documents the raw S7 constructor, which takes the parent's
properties and validates none of the relationships between them.

## Usage

``` r
StudentT1Distrib(
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

An S7 object of class `StudentT1Distrib`, inheriting from
`continuous_distrib` and from `distrib`. Its properties are the
parent's: `distrib_name`, `dimension`, `bounds`, `params`,
`params_interpretation`, `n_params`, `params_bounds`, `link_params` and
`params_smooth`. For an object built by
[`student_t1_distrib()`](https://statmodels7.github.io/distributions7/reference/student_t1_distrib.md)
they hold `"student t1"`, `"univariate"`, `c(-Inf, Inf)`,
`c("mu", "sigma", "nu")`, the interpretations
`c(mu = "location", sigma = "scale", nu = "shape")`, `3`, the domains
\\(-\infty, \infty)\\ and \\(0, \infty)\\ twice, and the three links.

## Methods

Registered on this class in this file:
[`distrib_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.StudentT1Distrib.md),
[`distrib_deriv3()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.StudentT1Distrib.md),
[`distrib_deriv4()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.StudentT1Distrib.md),
[`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.StudentT1Distrib.md),
[`distrib_grad_y()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.StudentT1Distrib.md),
[`distrib_gradient()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.StudentT1Distrib.md),
[`distrib_hess_y()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.StudentT1Distrib.md),
[`distrib_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.StudentT1Distrib.md),
[`distrib_pdf()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.StudentT1Distrib.md),
[`distrib_quantile()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.StudentT1Distrib.md),
[`distrib_rng()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.StudentT1Distrib.md).

Nine more are registered on the class from other files: the mixed
derivatives
[`distrib_cross_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross_y.StudentT1Distrib.md)
and
[`distrib_cross2_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross2_y.StudentT1Distrib.md),
the distribution-function derivatives
[`distrib_grad_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_cdf.StudentT1Distrib.md)
and
[`distrib_hess_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_cdf.StudentT1Distrib.md),
[`distrib_grad_y_hess()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y_hess.md)
and
[`distrib_hess_y_hess()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y_hess.md)
in `theta2_families.R`, and the four moments
[`mean()`](https://statmodels7.github.io/distributions7/reference/mean.StudentT1Distrib.md),
[`variance()`](https://statmodels7.github.io/distributions7/reference/variance.StudentT1Distrib.md),
[`skewness()`](https://statmodels7.github.io/distributions7/reference/skewness.StudentT1Distrib.md)
and
[`kurtosis()`](https://statmodels7.github.io/distributions7/reference/kurtosis.StudentT1Distrib.md)
in `moments.R`.

Everything else is inherited from
[`continuous_distrib()`](https://statmodels7.github.io/distributions7/reference/continuous_distrib.md).

## See also

[`student_t1_distrib()`](https://statmodels7.github.io/distributions7/reference/student_t1_distrib.md)
to build one;
[`student_t2_distrib()`](https://statmodels7.github.io/distributions7/reference/student_t2_distrib.md)
for the same law parametrized by the standard deviation;
[`gaussian1_distrib()`](https://statmodels7.github.io/distributions7/reference/gaussian1_distrib.md)
for the limit at \\\nu \to \infty\\;
[`cauchy_distrib()`](https://statmodels7.github.io/distributions7/reference/cauchy_distrib.md)
for the case \\\nu = 1\\;
[`distrib_gradient.StudentT1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.StudentT1Distrib.md)
for the redescending score.

## Examples

``` r
d <- student_t1_distrib()
S7::S7_inherits(d, continuous_distrib)
#> [1] TRUE

# The properties a consumer reads to drive the family without knowing it.
d@params
#> [1] "mu"    "sigma" "nu"   
d@params_interpretation
#>         mu      sigma         nu 
#> "location"    "scale"    "shape" 
d@params_bounds
#> $mu
#> [1] -Inf  Inf
#> 
#> $sigma
#> [1]   0 Inf
#> 
#> $nu
#> [1]   0 Inf
#> 

# The location is already free; the scale and the degrees of freedom ride
# a log, both being positive.
vapply(d@link_params, function(l) l@link_name, character(1))
#>         mu      sigma         nu 
#> "identity"      "log"      "log" 
```
