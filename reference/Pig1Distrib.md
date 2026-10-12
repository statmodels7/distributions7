# Poisson-Inverse Gaussian Distribution Class

The S7 class of the Poisson-inverse Gaussian on \\\\0, 1, 2, \dots\\\\
in its mean-dispersion parametrization, gamlss's `PIG`: the mean is
\\\mu\\ and the variance \\\mu + \sigma\mu^2\\. The family is a Poisson
mixed over an inverse Gaussian rate, an overdispersed count model with a
heavier tail than the negative binomial at the same variance.

Build one with
[`pig1_distrib()`](https://statmodels7.github.io/distributions7/reference/pig1_distrib.md),
which supplies the two link functions. This page documents the raw S7
constructor, which takes the parent's properties; the validator of
`distrib` checks them, and the constructor adds no check of its own.

## Usage

``` r
Pig1Distrib(
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

An S7 object of class `Pig1Distrib`, inheriting from `discrete_distrib`
and from `distrib`. For an object built by
[`pig1_distrib()`](https://statmodels7.github.io/distributions7/reference/pig1_distrib.md)
the properties hold `"poisson-inverse gaussian"`, `"univariate"`,
`c(0, Inf)`, `c("mu", "sigma")`, the interpretations
`c(mu = "mean", sigma = "dispersion")`, `2`, and two domains \\(0,
\infty)\\.

## Methods

Registered on this class, each linked to the page that documents it:

- the distribution:
  [`distrib_pdf()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.Pig1Distrib.md)
  and
  [`distrib_rng()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.Pig1Distrib.md);

- the derivatives in the parameters:
  [`distrib_gradient()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.Pig1Distrib.md),
  [`distrib_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.Pig1Distrib.md),
  [`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.Pig2Distrib.md),
  [`distrib_deriv3()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.Pig1Distrib.md)
  and
  [`distrib_deriv4()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.Pig1Distrib.md);

- the derivatives of the expected information:
  [`distrib_dexpected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.Pig2Distrib.md)
  and
  [`distrib_d2expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.Pig2Distrib.md);

- the moments:
  [`mean()`](https://statmodels7.github.io/distributions7/reference/mean.Pig1Distrib.md),
  [`variance()`](https://statmodels7.github.io/distributions7/reference/variance.Pig1Distrib.md),
  [`skewness()`](https://statmodels7.github.io/distributions7/reference/skewness.Pig1Distrib.md)
  and
  [`kurtosis()`](https://statmodels7.github.io/distributions7/reference/kurtosis.Pig1Distrib.md);

- the fitting helpers:
  [`distrib_start()`](https://statmodels7.github.io/distributions7/reference/distrib_start.Pig1Distrib.md)
  and
  [`expected_hessian_costly()`](https://statmodels7.github.io/distributions7/reference/expected_hessian_costly.md).

Every other generic is inherited from
[`discrete_distrib()`](https://statmodels7.github.io/distributions7/reference/discrete_distrib.md).
The distribution function and the quantile are the inherited ones,
documented at
[`distrib_cdf.discrete_distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.discrete_distrib.md)
and
[`distrib_quantile.discrete_distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.discrete_distrib.md).

## See also

[`pig1_distrib()`](https://statmodels7.github.io/distributions7/reference/pig1_distrib.md)
to build one;
[`pig2_distrib()`](https://statmodels7.github.io/distributions7/reference/pig2_distrib.md)
for the parametrization whose two parameters are orthogonal;
[`negbin2_distrib()`](https://statmodels7.github.io/distributions7/reference/negbin2_distrib.md)
for the other overdispersed count family.

## Examples

``` r
d <- pig1_distrib()
S7::S7_inherits(d, discrete_distrib)
#> [1] TRUE

d@params
#> [1] "mu"    "sigma"
d@params_interpretation
#>           mu        sigma 
#>       "mean" "dispersion" 

# The mean is a parameter and the variance follows from both.
th <- list(mu = 3, sigma = 0.8)
c(mean = mean(d, th), variance = variance(d, th),
  formula = 3 + 0.8 * 3^2)
#>     mean variance  formula 
#>      3.0     10.2     10.2 

# The mass sums to one over the support.
sum(distrib_pdf(d, 0:300, th))
#> [1] 1
```
