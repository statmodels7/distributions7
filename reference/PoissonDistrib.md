# S7 Class for Poisson Distribution

A subclass of `discrete_distrib` representing the Poisson distribution.

## Usage

``` r
PoissonDistrib(
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

An object of class `PoissonDistrib`.

## Methods

Methods implemented for this class:
[`distrib_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.PoissonDistrib.md),
[`distrib_deriv3()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.PoissonDistrib.md),
[`distrib_deriv4()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.PoissonDistrib.md),
[`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.PoissonDistrib.md),
[`distrib_gradient()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.PoissonDistrib.md),
[`distrib_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.PoissonDistrib.md),
[`distrib_pdf()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.PoissonDistrib.md),
[`distrib_quantile()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.PoissonDistrib.md),
[`distrib_rng()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.PoissonDistrib.md)

Everything else is inherited from
[`discrete_distrib`](https://statmodels7.github.io/distributions7/reference/discrete_distrib.md).

## See also

[`poisson_distrib`](https://statmodels7.github.io/distributions7/reference/poisson_distrib.md)
