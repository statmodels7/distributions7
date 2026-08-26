# Exponential Distribution Class

The S7 class of the exponential family parametrized by its **mean**
\\\mu \> 0\\, with density \\f(y) = \mu^{-1}e^{-y/\mu}\\ on \\\[0,
\infty)\\. It inherits from `continuous_distrib`, so it answers every
generic of the `distrib` contract; the eleven methods listed below are
registered on it directly.

This is the only single-parameter continuous family in the package, so
every derivative array it returns has one component per order and its
information is a number.

The parametrization is by the mean, not by the rate: R's own `dexp`
takes `rate`, and the methods pass `rate = 1/mu`. The mean is the scale
of the distribution, so \\\mu\\ is also its standard deviation.

Build one with
[`exponential_distrib()`](https://statmodels7.github.io/distributions7/reference/exponential_distrib.md).
This page documents the raw S7 constructor, which takes the parent's
properties and validates none of the relationships between them.

## Usage

``` r
ExponentialDistrib(
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

An S7 object of class `ExponentialDistrib`, inheriting from
`continuous_distrib` and from `distrib`. Its properties are the
parent's: `distrib_name`, `dimension`, `bounds`, `params`,
`params_interpretation`, `n_params`, `params_bounds`, `link_params` and
`params_smooth`. For an object built by
[`exponential_distrib()`](https://statmodels7.github.io/distributions7/reference/exponential_distrib.md)
they hold `"exponential"`, `"univariate"`, `c(0, Inf)`, `"mu"`,
`c(mu = "mean")`, `1`, the domain \\(0, \infty)\\, and the one link.

## Methods

Registered on this class:
[`distrib_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.ExponentialDistrib.md),
[`distrib_deriv3()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.ExponentialDistrib.md),
[`distrib_deriv4()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.ExponentialDistrib.md),
[`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.ExponentialDistrib.md),
[`distrib_grad_y()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.ExponentialDistrib.md),
[`distrib_gradient()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.ExponentialDistrib.md),
[`distrib_hess_y()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.ExponentialDistrib.md),
[`distrib_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.ExponentialDistrib.md),
[`distrib_pdf()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.ExponentialDistrib.md),
[`distrib_quantile()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.ExponentialDistrib.md),
[`distrib_rng()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.ExponentialDistrib.md)

Everything else is inherited from
[`continuous_distrib()`](https://statmodels7.github.io/distributions7/reference/continuous_distrib.md).

## See also

[`exponential_distrib()`](https://statmodels7.github.io/distributions7/reference/exponential_distrib.md)
to build one;
[`gamma1_distrib()`](https://statmodels7.github.io/distributions7/reference/gamma1_distrib.md)
and
[`weibull1_distrib()`](https://statmodels7.github.io/distributions7/reference/weibull1_distrib.md),
both of which contain this family at a unit shape;
[`geometric_distrib()`](https://statmodels7.github.io/distributions7/reference/geometric_distrib.md)
for its discrete analogue.

## Examples

``` r
d <- exponential_distrib()
S7::S7_inherits(d, continuous_distrib)
#> [1] TRUE

# One parameter, the mean, on the positive half line.
d@params
#> [1] "mu"
d@n_params
#> [1] 1
d@bounds
#> [1]   0 Inf

# The mean is also the standard deviation, and the shape is fixed: the
# skewness is 2 and the excess kurtosis 6 at every mu.
th <- list(mu = 2)
c(mean = mean(d, th), sd = std_dev(d, th),
  skew = skewness(d, th), kurt = kurtosis(d, th))
#> mean   sd skew kurt 
#>    2    2    2    6 
```
