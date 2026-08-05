# S7 Class for Probability Distributions

The base S7 class for probability distributions. It carries the name,
the parameters with their domains and links, and the support.

## Usage

``` r
distrib(
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

An object of class `distrib`.

## Methods

Registered on the base class, so every distribution inherits them unless
it registers something more specific. Those that compute derivatives do
so by finite differences, which is why a subclass that implements
nothing but
[`distrib_pdf`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.md)
is already fully functional. Note that `distrib_pdf` itself has no
default: the density is the one thing a distribution must supply.

[`distrib_deriv3()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.distrib.md),
[`distrib_deriv4()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.distrib.md),
[`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.distrib.md),
[`distrib_gradient()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.distrib.md),
[`distrib_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.distrib.md),
[`generate_random_theta()`](https://statmodels7.github.io/distributions7/reference/generate_random_theta.distrib.md),
[`kurtosis()`](https://statmodels7.github.io/distributions7/reference/kurtosis.md),
[`mean()`](https://statmodels7.github.io/distributions7/reference/mean.distrib.md),
[`print()`](https://statmodels7.github.io/distributions7/reference/print.distrib.md),
[`skewness()`](https://statmodels7.github.io/distributions7/reference/skewness.md),
[`std_dev()`](https://statmodels7.github.io/distributions7/reference/std_dev.md),
[`variance()`](https://statmodels7.github.io/distributions7/reference/variance.md)

## Examples

``` r
d <- gaussian1_distrib()
S7::S7_inherits(d, distrib)
#> [1] TRUE
d@params
#> [1] "mu"    "sigma"
d@params_bounds
#> $mu
#> [1] -Inf  Inf
#> 
#> $sigma
#> [1]   0 Inf
#> 
```
