# S7 Class for the Elastic-Net Distribution

A subclass of `continuous_distrib` for the density whose negative
logarithm is the elastic-net penalty, \\f(y) \propto
\exp\\-\lambda\alpha\|y-\mu\| - \lambda(1-\alpha)(y-\mu)^2/2\\\\: the
product of a Laplace and a Gaussian, normalized. At \\\alpha \to 1\\ it
is the Laplace of
[`Laplace2Distrib`](https://statmodels7.github.io/distributions7/reference/Laplace2Distrib.md)
and at \\\alpha \to 0\\ the Gaussian, both of which remain families of
their own; \\\alpha\\ is confined to the open interval, as every bounded
parameter in this package is. Like the Laplace, its log-likelihood is
not differentiable in \\\mu\\.

## Usage

``` r
EnetDistrib(
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

An object of class `EnetDistrib`.

## Methods

Methods implemented for this class:
[`distrib_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.EnetDistrib.md),
[`distrib_cross_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross_y.EnetDistrib.md),
[`distrib_grad_y()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.EnetDistrib.md),
[`distrib_gradient()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.EnetDistrib.md),
[`distrib_hess_y()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.EnetDistrib.md),
[`distrib_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.EnetDistrib.md),
[`distrib_pdf()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.EnetDistrib.md),
[`distrib_quantile()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.EnetDistrib.md),
[`distrib_rng()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.EnetDistrib.md),
[`mean()`](https://statmodels7.github.io/distributions7/reference/mean.distrib.md),
[`skewness()`](https://statmodels7.github.io/distributions7/reference/skewness.md),
[`variance()`](https://statmodels7.github.io/distributions7/reference/variance.md)

Everything else is inherited from
[`continuous_distrib`](https://statmodels7.github.io/distributions7/reference/continuous_distrib.md).

## See also

[`enet_distrib`](https://statmodels7.github.io/distributions7/reference/enet_distrib.md),
[`laplace2_distrib`](https://statmodels7.github.io/distributions7/reference/laplace2_distrib.md)
