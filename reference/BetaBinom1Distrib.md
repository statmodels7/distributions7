# Beta-Binomial Distribution Class, Mean Proportion and Dispersion

The S7 class of the beta-binomial family parametrized by a mean
proportion \\\mu \in (0, 1)\\ and a dispersion \\\sigma \> 0\\, on the
finite support \\\\0, 1, \dots, n\\\\. It inherits from
`discrete_distrib`, so it answers every generic of the `distrib`
contract; the seven methods listed below are registered on it in this
file, two more in `betabinom1_higher.R`, and everything else comes from
the parent.

The class carries an extra property beyond the parent's, `size`: the
number of trials \\n\\, fixed at construction as it is for
[`BinomialDistrib()`](https://statmodels7.github.io/distributions7/reference/BinomialDistrib.md).
Build one with
[`betabinom1_distrib()`](https://statmodels7.github.io/distributions7/reference/betabinom1_distrib.md),
which validates `size`, supplies the two link functions and fills the
properties in. This page documents the raw S7 constructor, which
validates none of the relationships between them.

## Usage

``` r
BetaBinom1Distrib(
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

  The number of trials \\n\\, a single positive integer stored as a
  numeric. It belongs to the object, so an object cannot be reused
  across data sets whose group sizes differ.

## Value

An S7 object of class `BetaBinom1Distrib`, inheriting from
`discrete_distrib` and from `distrib`. Beyond `size` its properties are
the parent's: `distrib_name`, `dimension`, `bounds`, `params`,
`params_interpretation`, `n_params`, `params_bounds`, `link_params` and
`params_smooth`. For an object built by
[`betabinom1_distrib()`](https://statmodels7.github.io/distributions7/reference/betabinom1_distrib.md)
they hold `"beta-binomial [size=n]"`, `"univariate"`, `c(0, size)`,
`c("mu", "sigma")`, the interpretations
`c(mu = "mean proportion", sigma = "dispersion")`, `2`, the domains
\\(0, 1)\\ and \\(0, \infty)\\, and the two links.

## Methods

Registered on this class in this file:
[`distrib_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.BetaBinom1Distrib.md),
[`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.BetaBinom1Distrib.md),
[`distrib_gradient()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.BetaBinom1Distrib.md),
[`distrib_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.BetaBinom1Distrib.md),
[`distrib_pdf()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.BetaBinom1Distrib.md),
[`distrib_quantile()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.BetaBinom1Distrib.md),
[`distrib_rng()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.BetaBinom1Distrib.md).

Six more are registered on the class from other files: the closed-form
[`distrib_deriv3()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.BetaBinom1Distrib.md)
and
[`distrib_deriv4()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.BetaBinom1Distrib.md)
in `betabinom1_higher.R`, and the four moments
[`mean()`](https://statmodels7.github.io/distributions7/reference/mean.BetaBinom1Distrib.md),
[`variance()`](https://statmodels7.github.io/distributions7/reference/variance.BetaBinom1Distrib.md),
[`skewness()`](https://statmodels7.github.io/distributions7/reference/skewness.BetaBinom1Distrib.md)
and
[`kurtosis()`](https://statmodels7.github.io/distributions7/reference/kurtosis.BetaBinom1Distrib.md)
in `moments.R`.

The response derivatives are refused, as for every discrete family: a
mass function has no derivative in its argument.

## See also

[`betabinom1_distrib()`](https://statmodels7.github.io/distributions7/reference/betabinom1_distrib.md)
to build one;
[`betabinom2_distrib()`](https://statmodels7.github.io/distributions7/reference/betabinom2_distrib.md)
for the same law in its two beta shapes;
[`binomial_distrib()`](https://statmodels7.github.io/distributions7/reference/binomial_distrib.md)
for the limit at \\\sigma \to 0\\;
[`distrib_pdf.BetaBinom1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.BetaBinom1Distrib.md)
for the mass function.

## Examples

``` r
d <- betabinom1_distrib(size = 10)
S7::S7_inherits(d, discrete_distrib)
#> [1] TRUE

# The trial count is a property of the object, not an entry of theta.
d@size
#> [1] 10
d@bounds
#> [1]  0 10

# The properties a consumer reads to drive the family without knowing it.
d@params
#> [1] "mu"    "sigma"
d@params_interpretation
#>                mu             sigma 
#> "mean proportion"      "dispersion" 

# The mean proportion rides a logit and the dispersion a log.
vapply(d@link_params, function(l) l@link_name, character(1))
#>      mu   sigma 
#> "logit"   "log" 
```
