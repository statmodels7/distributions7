# Poisson Distribution Class

The S7 class of the Poisson family parametrized by its mean \\\mu \>
0\\, with mass \\P(Y = y) = e^{-\mu}\mu^y/y!\\ on the non-negative
integers. It inherits from `discrete_distrib`, so expectations over its
support are exact sums and no derivative with respect to the response is
defined.

The mean and the variance are both \\\mu\\, which is the constraint a
count model most often has to relax;
[`negbin2_distrib()`](https://statmodels7.github.io/distributions7/reference/negbin2_distrib.md)
and
[`negbin1_distrib()`](https://statmodels7.github.io/distributions7/reference/negbin1_distrib.md)
do that, and both contain this family in the limit.

The default link is the logarithm, which is the **canonical** link here.
On its scale the observed and the expected information coincide, so
Fisher scoring and Newton's method take the same step.

Build one with
[`poisson_distrib()`](https://statmodels7.github.io/distributions7/reference/poisson_distrib.md).
This page documents the raw S7 constructor, which takes the parent's
properties and validates none of the relationships between them.

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
  [`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.md)).

## Value

An S7 object of class `PoissonDistrib`, inheriting from
`discrete_distrib` and from `distrib`. Its properties are the parent's:
`distrib_name`, `dimension`, `bounds`, `params`,
`params_interpretation`, `n_params`, `params_bounds`, `link_params` and
`params_smooth`. For an object built by
[`poisson_distrib()`](https://statmodels7.github.io/distributions7/reference/poisson_distrib.md)
they hold `"poisson"`, `"univariate"`, `c(0, Inf)`, `"mu"`,
`c(mu = "mean")`, `1`, the domain \\(0, \infty)\\, and the one link.

## Methods

Registered on this class:
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
[`discrete_distrib()`](https://statmodels7.github.io/distributions7/reference/discrete_distrib.md).
Derivatives with respect to the response are refused by that parent: the
support is a lattice, so there is nothing to differentiate along.

## See also

[`poisson_distrib()`](https://statmodels7.github.io/distributions7/reference/poisson_distrib.md)
to build one;
[`negbin2_distrib()`](https://statmodels7.github.io/distributions7/reference/negbin2_distrib.md)
and
[`negbin1_distrib()`](https://statmodels7.github.io/distributions7/reference/negbin1_distrib.md)
for overdispersed counts;
[`zero_inflated()`](https://statmodels7.github.io/distributions7/reference/zero_inflated.md)
and
[`zero_adjusted()`](https://statmodels7.github.io/distributions7/reference/zero_adjusted.md)
for excess zeros;
[`distrib_expected_hessian.PoissonDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.PoissonDistrib.md)
for the information.

## Examples

``` r
d <- poisson_distrib()
S7::S7_inherits(d, discrete_distrib)
#> [1] TRUE
d@params
#> [1] "mu"
d@bounds
#> [1]   0 Inf

# Equidispersion: the mean is the variance at every parameter value.
vapply(c(0.5, 3, 20), function(m) {
  th <- list(mu = m)
  c(mean = mean(d, th), var = variance(d, th))
}, numeric(2))
#>      [,1] [,2] [,3]
#> mean  0.5    3   20
#> var   0.5    3   20

# The log link is canonical, so on its scale the observed and the expected
# information are the same number.
th <- list(mu = 3)
rbind(observed = distrib_hessian(d, c(0, 2, 7), th, scale = "link")$mu_mu,
      expected = distrib_expected_hessian(d, c(0, 2, 7), th,
                                          scale = "link")$mu_mu)
#>          [,1] [,2] [,3]
#> observed   -3   -3   -3
#> expected   -3   -3   -3
```
