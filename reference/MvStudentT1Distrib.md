# S7 Classes for the Two Multivariate Student t Parametrizations

`MvStudentT1Distrib` is the family whose matrix parametrization is the
**scale matrix** \\\Sigma\\ and `MvStudentT2Distrib` the family whose
parametrization is its **inverse** \\\Sigma^{-1}\\. Both inherit from
[MvStudentTDistrib](https://statmodels7.github.io/distributions7/reference/MvStudentTDistrib.md),
which carries every method: the two differ in what their free values
describe, not in how anything is computed.

## Usage

``` r
MvStudentT1Distrib(
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

MvStudentT2Distrib(
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

  A parameters7 parametrization of the scale matrix or of its inverse,
  according to `inverted`, inheriting from
  [`parameters7::parameter`](https://statmodels7.github.io/parameters7/reference/parameter.html).
  Its `n_free` free values are flattened into scalar parameters of the
  distribution.

- inverted:

  Logical of length 1. `TRUE` when `param` carries the inverse scale
  matrix \\\Sigma^{-1}\\ and `FALSE` when it carries the scale matrix
  \\\Sigma\\. Nothing but the sign of the log-determinant and one matrix
  inversion depends on it; the law is the same either way, and neither
  matrix is a moment of the response.

## Value

An object of class `MvStudentT1Distrib` or `MvStudentT2Distrib`, a
subclass of
[MvStudentTDistrib](https://statmodels7.github.io/distributions7/reference/MvStudentTDistrib.md)
adding no properties of its own.

## Details

Neither matrix is a covariance or a precision of the response, and the
difference is a factor in \\\nu\\ that
[`mvstudent_t1_distrib()`](https://statmodels7.github.io/distributions7/reference/mvstudent_t1_distrib.md)
writes out. The numbering is the package's convention, the one that
separates
[`mvgaussian1_distrib()`](https://statmodels7.github.io/distributions7/reference/mvgaussian1_distrib.md)
from
[`mvgaussian2_distrib()`](https://statmodels7.github.io/distributions7/reference/mvgaussian1_distrib.md)
and
[`gaussian1_distrib()`](https://statmodels7.github.io/distributions7/reference/gaussian1_distrib.md)
from
[`gaussian3_distrib()`](https://statmodels7.github.io/distributions7/reference/gaussian3_distrib.md):
a parametrization decides what a linear predictor acts on, so it gets a
name.

## See also

[`mvstudent_t1_distrib()`](https://statmodels7.github.io/distributions7/reference/mvstudent_t1_distrib.md)
and
[`mvstudent_t2_distrib()`](https://statmodels7.github.io/distributions7/reference/mvstudent_t1_distrib.md),
the constructors, and
[MvStudentTDistrib](https://statmodels7.github.io/distributions7/reference/MvStudentTDistrib.md)
for the properties and the methods.

## Examples

``` r
c(scale = S7::S7_inherits(mvstudent_t1_distrib(2), MvStudentT1Distrib),
  inverse = S7::S7_inherits(mvstudent_t2_distrib(2), MvStudentT2Distrib))
#>   scale inverse 
#>    TRUE    TRUE 
```
