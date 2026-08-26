# Beta-Binomial Distribution Class, Two Shapes

The S7 class of the beta-binomial family in its canonical
parametrization, the two beta shapes \\\alpha \> 0\\ and \\\beta \> 0\\,
on the finite support \\\\0, 1, \dots, n\\\\. It inherits from
`discrete_distrib`, so it answers every generic of the `distrib`
contract; the seven methods listed below are registered on it directly
and everything else comes from the parent.

The class carries an extra property beyond the parent's, `size`: the
number of trials \\n\\, fixed at construction as it is for
[`BinomialDistrib()`](https://statmodels7.github.io/distributions7/reference/BinomialDistrib.md).
Build one with
[`betabinom2_distrib()`](https://statmodels7.github.io/distributions7/reference/betabinom2_distrib.md),
which validates `size`, supplies the two link functions and fills the
properties in. This page documents the raw S7 constructor, which
validates none of the relationships between them.

## Usage

``` r
BetaBinom2Distrib(
  distrib_name = character(0),
  dimension = character(0),
  bounds = integer(0),
  params = character(0),
  params_interpretation = character(0),
  n_params = integer(0),
  params_bounds = list(),
  link_params = list(),
  params_smooth = logical(0),
  size = integer(0)
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

- size:

  The number of trials \\n\\, a single positive integer. It belongs to
  the object, so an object cannot be reused across data sets whose group
  sizes differ.

## Value

An S7 object of class `BetaBinom2Distrib`, inheriting from
`discrete_distrib` and from `distrib`. Beyond `size` its properties are
the parent's: `distrib_name`, `dimension`, `bounds`, `params`,
`params_interpretation`, `n_params`, `params_bounds`, `link_params` and
`params_smooth`. For an object built by
[`betabinom2_distrib()`](https://statmodels7.github.io/distributions7/reference/betabinom2_distrib.md)
they hold `"betabinom2 [size=n]"`, `"univariate"`, `c(0, size)`,
`c("alpha", "beta")`, the interpretations
`c(alpha = "shape", beta = "shape")`, `2`, the domain \\(0, \infty)\\
for both, and the two links.

## Methods

Registered on this class:
[`distrib_deriv3()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.BetaBinom2Distrib.md),
[`distrib_deriv4()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.BetaBinom2Distrib.md),
[`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.BetaBinom2Distrib.md),
[`distrib_gradient()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.BetaBinom2Distrib.md),
[`distrib_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.BetaBinom2Distrib.md),
[`distrib_pdf()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.BetaBinom2Distrib.md),
[`distrib_rng()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.BetaBinom2Distrib.md).

The four moments
[`mean()`](https://statmodels7.github.io/distributions7/reference/mean.BetaBinom2Distrib.md),
[`variance()`](https://statmodels7.github.io/distributions7/reference/variance.BetaBinom2Distrib.md),
[`skewness()`](https://statmodels7.github.io/distributions7/reference/skewness.BetaBinom2Distrib.md)
and
[`kurtosis()`](https://statmodels7.github.io/distributions7/reference/kurtosis.BetaBinom2Distrib.md)
are registered in `moments.R`.

The distribution and quantile functions come from
[`discrete_distrib()`](https://statmodels7.github.io/distributions7/reference/discrete_distrib.md)
rather than from this class, and on a finite support the parent's
cumulative sum of the mass is exact.

## See also

[`betabinom2_distrib()`](https://statmodels7.github.io/distributions7/reference/betabinom2_distrib.md)
to build one;
[`betabinom1_distrib()`](https://statmodels7.github.io/distributions7/reference/betabinom1_distrib.md)
for the same law in a mean proportion and a dispersion;
[`beta1_distrib()`](https://statmodels7.github.io/distributions7/reference/beta1_distrib.md)
for the mixing law;
[`distrib_pdf.BetaBinom2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.BetaBinom2Distrib.md)
for the mass function.

## Examples

``` r
d <- betabinom2_distrib(size = 10)
S7::S7_inherits(d, discrete_distrib)
#> [1] TRUE

# The trial count is a property of the object, not an entry of theta.
d@size
#> [1] 10
d@bounds
#> [1]  0 10

# Both shapes are positive, so both ride a log by default.
d@params
#> [1] "alpha" "beta" 
vapply(d@link_params, function(l) l@link_name, character(1))
#> alpha  beta 
#> "log" "log" 
```
