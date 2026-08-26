# Negative Binomial Distribution Class, NB1

The S7 class of the negative binomial family on the non-negative
integers whose variance is **linear** in the mean: with a mean \\\mu \>
0\\ and a dispersion \\\theta \> 0\\, \\\operatorname{Var}(Y) =
\mu(1+\theta)\\, so the variance-to-mean ratio is \\1+\theta\\ at every
mean. It inherits from `discrete_distrib`, so it answers every generic
of the `distrib` contract; the seven methods listed below are registered
on it in this file and everything else comes from the parent.

Build one with
[`negbin1_distrib()`](https://statmodels7.github.io/distributions7/reference/negbin1_distrib.md),
which supplies the two link functions and fills the properties in. This
page documents the raw S7 constructor, which takes the parent's
properties and validates none of the relationships between them.

## Usage

``` r
NegBin1Distrib(
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

An S7 object of class `NegBin1Distrib`, inheriting from
`discrete_distrib` and from `distrib`. Its properties are the parent's:
`distrib_name`, `dimension`, `bounds`, `params`,
`params_interpretation`, `n_params`, `params_bounds`, `link_params` and
`params_smooth`. For an object built by
[`negbin1_distrib()`](https://statmodels7.github.io/distributions7/reference/negbin1_distrib.md)
they hold `"negbin1"`, `"univariate"`, `c(0, Inf)`, `c("mu", "theta")`,
the interpretations `c(mu = "mean", theta = "dispersion")`, `2`, the
domain \\(0, \infty)\\ for both parameters, and the two links.

## Methods

Registered on this class in this file:
[`distrib_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.NegBin1Distrib.md),
[`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.NegBin1Distrib.md),
[`distrib_gradient()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.NegBin1Distrib.md),
[`distrib_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.NegBin1Distrib.md),
[`distrib_pdf()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.NegBin1Distrib.md),
[`distrib_quantile()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.NegBin1Distrib.md),
[`distrib_rng()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.NegBin1Distrib.md)

The third and fourth orders,
[`distrib_deriv3()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.NegBin1Distrib.md)
and
[`distrib_deriv4()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.NegBin1Distrib.md),
are registered elsewhere in the package, on top of
[`negbin1_components()`](https://statmodels7.github.io/distributions7/reference/negbin1_components.md).
A discrete family has no derivatives in the response. Everything else is
inherited from
[`discrete_distrib()`](https://statmodels7.github.io/distributions7/reference/discrete_distrib.md).

## See also

[`negbin1_distrib()`](https://statmodels7.github.io/distributions7/reference/negbin1_distrib.md)
to build one;
[`negbin2_distrib()`](https://statmodels7.github.io/distributions7/reference/negbin2_distrib.md)
for the quadratic-variance family, which is a different family and not a
reparametrization of this one;
[`poisson_distrib()`](https://statmodels7.github.io/distributions7/reference/poisson_distrib.md)
for the limit as \\\theta\\ goes to zero;
[`distrib_pdf.NegBin1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.NegBin1Distrib.md)
and
[`distrib_gradient.NegBin1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.NegBin1Distrib.md)
for the closed forms this class supplies.

## Examples

``` r
d <- negbin1_distrib()
S7::S7_inherits(d, discrete_distrib)
#> [1] TRUE

# The properties a consumer reads to drive the family without knowing it.
d@params
#> [1] "mu"    "theta"
d@params_interpretation
#>           mu        theta 
#>       "mean" "dispersion" 
d@bounds
#> [1]   0 Inf

# The variance-to-mean ratio is 1 + theta at every mean, which is what
# separates this family from the quadratic one.
vapply(c(1, 10, 100),
       function(m) variance(d, list(mu = m, theta = 4)) / m, numeric(1))
#> [1] 5 5 5
```
