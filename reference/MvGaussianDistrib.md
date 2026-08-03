# Multivariate Gaussian Distribution

The S7 class of multivariate gaussian distributions, parametrised by a
mean vector and by a covstructs7 structure for the covariance or the
precision. Constructed by
[`mvgaussian_distrib`](https://statmodels7.github.io/distributions7/reference/mvgaussian_distrib.md).

## Usage

``` r
MvGaussianDistrib(
  distrib_name = character(0),
  dimension = character(0),
  bounds = integer(0),
  params = character(0),
  params_interpretation = character(0),
  n_params = integer(0),
  params_bounds = list(),
  link_params = list(),
  params_smooth = logical(0),
  n_dim = integer(0),
  struct = covstructs7::covstruct(),
  inverted = logical(0)
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

- n_dim:

  The dimension \\p\\ of an observation.

- struct:

  The covstructs7 structure carrying the matrix.

- inverted:

  Whether the structure parametrises the precision rather than the
  covariance.

## Value

An object of class `MvGaussianDistrib`.

## See also

[`mvgaussian_distrib`](https://statmodels7.github.io/distributions7/reference/mvgaussian_distrib.md)

## Examples

``` r
S7::S7_inherits(mvgaussian_distrib(2), MvGaussianDistrib)
#> [1] TRUE
```
