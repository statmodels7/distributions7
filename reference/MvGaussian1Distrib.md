# S7 Classes for the Two Multivariate Gaussian Parametrizations

`MvGaussian1Distrib` is the family whose matrix parametrization is the
**covariance** and `MvGaussian2Distrib` the family whose matrix
parametrization is the **precision**. Both inherit from
[MvGaussianDistrib](https://statmodels7.github.io/distributions7/reference/MvGaussianDistrib.md),
which carries every method: the two differ in what their free values
describe, not in how anything is computed.

## Usage

``` r
MvGaussian1Distrib(
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
  param = parameters7::parameter(),
  inverted = logical(0)
)

MvGaussian2Distrib(
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
  param = parameters7::parameter(),
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
  [`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.md)).

- n_dim:

  The dimension \\p\\ of one observation. A single positive integer.

- param:

  A parameters7 parametrization of the matrix, inheriting from
  [`parameters7::parameter`](https://statmodels7.github.io/parameters7/reference/parameter.html).
  Its `n_free` free values are flattened into scalar parameters of the
  distribution.

- inverted:

  Logical of length 1. `TRUE` when `param` carries the precision
  \\\Omega = \Sigma^{-1}\\ and `FALSE` when it carries the covariance
  \\\Sigma\\. Nothing but the sign of the log-determinant and one matrix
  inversion depends on it; the law is the same either way.

## Value

An object of class `MvGaussian1Distrib` or `MvGaussian2Distrib`, a
subclass of
[MvGaussianDistrib](https://statmodels7.github.io/distributions7/reference/MvGaussianDistrib.md)
adding no properties of its own.

## Details

The two are separate families for the reason
[`gaussian1_distrib()`](https://statmodels7.github.io/distributions7/reference/gaussian1_distrib.md),
[`gaussian2_distrib()`](https://statmodels7.github.io/distributions7/reference/gaussian2_distrib.md)
and
[`gaussian3_distrib()`](https://statmodels7.github.io/distributions7/reference/gaussian3_distrib.md)
are separate families: a parametrization is a modeling decision, and a
predictor imposed on the coordinates of \\\Sigma\\ is a different model
from the same predictor imposed on the coordinates of \\\Omega\\. Where
the matrix parametrization is closed under inversion the two describe
the same laws in different coordinates; where it is not, they describe
different laws.

## See also

[`mvgaussian1_distrib()`](https://statmodels7.github.io/distributions7/reference/mvgaussian1_distrib.md)
and
[`mvgaussian2_distrib()`](https://statmodels7.github.io/distributions7/reference/mvgaussian1_distrib.md),
the constructors, and
[MvGaussianDistrib](https://statmodels7.github.io/distributions7/reference/MvGaussianDistrib.md)
for the properties and the methods.

## Examples

``` r
c(covariance = S7::S7_inherits(mvgaussian1_distrib(2), MvGaussian1Distrib),
  precision = S7::S7_inherits(mvgaussian2_distrib(2), MvGaussian2Distrib))
#> covariance  precision 
#>       TRUE       TRUE 

# Both are multivariate gaussians, which is what carries the methods.
S7::S7_inherits(mvgaussian2_distrib(2), MvGaussianDistrib)
#> [1] TRUE
```
