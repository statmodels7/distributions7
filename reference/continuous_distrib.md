# S7 Class for Continuous Distributions

A subclass of `distrib` specifically for continuous probability
distributions.

## Usage

``` r
continuous_distrib(
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
  [`distrib_expected_hessian`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.md)).

## Value

An object of class `continuous_distrib`.

## Methods

Defaults for continuous distributions, built from the density alone: the
cdf by quadrature, the quantile function by root finding, and the
generator by Generalized Ratio-of-Uniforms
([`rng_grou`](https://statmodels7.github.io/distributions7/reference/rng_grou.md))
or inverse transform when an analytical quantile is available.

[`distrib_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.continuous_distrib.md),
[`distrib_grad_y()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.continuous_distrib.md),
[`distrib_hess_y()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.continuous_distrib.md),
[`distrib_quantile()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.continuous_distrib.md),
[`distrib_rng()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.continuous_distrib.md),
[`expectation()`](https://statmodels7.github.io/distributions7/reference/expectation.md),
[`plot()`](https://statmodels7.github.io/distributions7/reference/plot.continuous_distrib.md)

Everything else is inherited from
[`distrib`](https://statmodels7.github.io/distributions7/reference/distrib.md).

## See also

[`distrib`](https://statmodels7.github.io/distributions7/reference/distrib.md),
[`discrete_distrib`](https://statmodels7.github.io/distributions7/reference/discrete_distrib.md),
[`multivariate_distrib`](https://statmodels7.github.io/distributions7/reference/multivariate_distrib.md)

## Examples

``` r
S7::S7_inherits(gaussian1_distrib(), continuous_distrib)
#> [1] TRUE
S7::S7_inherits(poisson_distrib(), continuous_distrib)
#> [1] FALSE
```
