# Geometric Distribution Class

The S7 class of the geometric family parametrized by its **mean** \\\mu
\> 0\\: the number of failures before the first success in independent
trials, each succeeding with probability \\p = 1/(1+\mu)\\. It inherits
from `discrete_distrib`, so its support is a set of isolated points, so
expectations are sums and no derivative with respect to the response is
defined.

The parametrization is by the mean, the quantity the modeling layer
above models: a link carries \\\mu\\ to the unconstrained scale and the
probability follows. R's own `dgeom` takes `prob`, and the methods
convert through
[`geom_prob()`](https://statmodels7.github.io/distributions7/reference/geom_prob.md).

The variance is \\\mu(1+\mu)\\, always above the mean, so this family is
overdispersed relative to a Poisson of the same mean. It is the negative
binomial at a dispersion of one, and it is the only discrete law that is
memoryless.

Build one with
[`geometric_distrib()`](https://statmodels7.github.io/distributions7/reference/geometric_distrib.md).
This page documents the raw S7 constructor, which takes the parent's
properties and validates none of the relationships between them.

## Usage

``` r
GeometricDistrib(
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

An S7 object of class `GeometricDistrib`, inheriting from
`discrete_distrib` and from `distrib`. Its properties are the parent's:
`distrib_name`, `dimension`, `bounds`, `params`,
`params_interpretation`, `n_params`, `params_bounds`, `link_params` and
`params_smooth`. For an object built by
[`geometric_distrib()`](https://statmodels7.github.io/distributions7/reference/geometric_distrib.md)
they hold `"geometric"`, `"univariate"`, `c(0, Inf)`, `"mu"`,
`c(mu = "mean")`, `1`, the domain \\(0, \infty)\\, and the one link.

## Methods

Registered on this class:
[`distrib_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.GeometricDistrib.md),
[`distrib_deriv3()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.GeometricDistrib.md),
[`distrib_deriv4()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.GeometricDistrib.md),
[`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.GeometricDistrib.md),
[`distrib_gradient()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.GeometricDistrib.md),
[`distrib_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.GeometricDistrib.md),
[`distrib_pdf()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.GeometricDistrib.md),
[`distrib_quantile()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.GeometricDistrib.md),
[`distrib_rng()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.GeometricDistrib.md)

Everything else is inherited from
[`discrete_distrib()`](https://statmodels7.github.io/distributions7/reference/discrete_distrib.md).

## See also

[`geometric_distrib()`](https://statmodels7.github.io/distributions7/reference/geometric_distrib.md)
to build one;
[`negbin2_distrib()`](https://statmodels7.github.io/distributions7/reference/negbin2_distrib.md),
of which this is the case \\\theta = 1\\;
[`exponential_distrib()`](https://statmodels7.github.io/distributions7/reference/exponential_distrib.md),
its memoryless continuous counterpart;
[`poisson_distrib()`](https://statmodels7.github.io/distributions7/reference/poisson_distrib.md)
for the equidispersed alternative.

## Examples

``` r
d <- geometric_distrib()
S7::S7_inherits(d, discrete_distrib)
#> [1] TRUE
d@params
#> [1] "mu"
d@params_interpretation
#>     mu 
#> "mean" 

# Overdispersed by construction: the variance is mu(1+mu), above the mean.
vapply(c(0.5, 3, 20), function(m) {
  th <- list(mu = m)
  c(mean = mean(d, th), var = variance(d, th))
}, numeric(2))
#>      [,1] [,2] [,3]
#> mean 0.50    3   20
#> var  0.75   12  420

# The mass is R's own at prob = 1/(1+mu).
all.equal(distrib_pdf(d, c(0, 2, 7), list(mu = 3)),
          dgeom(c(0, 2, 7), prob = 1 / 4))
#> [1] TRUE
```
