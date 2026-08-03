# S7 Class for Multivariate Distributions

A subclass of `distrib` for distributions whose observations are
vectors. The response is an \\n \times p\\ matrix, one row per
observation, and `n_dim` records \\p\\.

## Usage

``` r
multivariate_distrib(
  distrib_name = character(0),
  dimension = character(0),
  bounds = integer(0),
  params = character(0),
  params_interpretation = character(0),
  n_params = integer(0),
  params_bounds = list(),
  link_params = list(),
  params_smooth = logical(0),
  n_dim = integer(0)
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

## Value

An object of class `multivariate_distrib`.

## Details

The parameters of a multivariate distribution are still **scalars** as
far as the rest of the package is concerned. A mean vector contributes
\\p\\ of them and a covariance matrix contributes the free values of the
covstructs7 structure that parametrises it, so `theta` remains the named
list of numbers every generic already understands, and the derivative
bookkeeping –
[`deriv_names`](https://statmodels7.github.io/distributions7/reference/deriv_names.md),
the Hessian keys, the link scale,
[`fit_distrib`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
– needs no special case. The constraint on the matrix lives inside the
structure rather than in a link, which is why the links of a
multivariate distribution are all the identity: the free values are
already unconstrained.

This class is a sibling of
[`continuous_distrib`](https://statmodels7.github.io/distributions7/reference/continuous_distrib.md)
rather than a subclass of it, and deliberately. The defaults registered
there are one-dimensional – a cdf by quadrature over an interval, a
quantile by root finding, a generator by ratio-of-uniforms on a scalar
density – and none of them means anything in \\p\\ dimensions.
Inheriting them would offer answers that do not exist.

Parameters that vary from observation to observation are not supported
here. A distribution whose parameters depend on covariates is a model,
which is the layer above this one.

## See also

[`mvgaussian_distrib`](https://statmodels7.github.io/distributions7/reference/mvgaussian_distrib.md),
[`n_obs`](https://statmodels7.github.io/distributions7/reference/n_obs.md)

## Examples

``` r
d <- mvgaussian_distrib(2)
S7::S7_inherits(d, multivariate_distrib)
#> [1] TRUE
c(n_dim = d@n_dim, n_params = d@n_params)
#>    n_dim n_params 
#>        2        5 
```
