# Poisson-Inverse Gaussian Distribution Class

The S7 class of the Poisson-inverse Gaussian on \\\\0, 1, 2, \dots\\\\
in its mean-dispersion parametrization, gamlss's `PIG`: the mean is
\\\mu\\ and the variance \\\mu + \sigma\mu^2\\. The family is a Poisson
mixed over an inverse Gaussian rate, an overdispersed count model with a
heavier tail than the negative binomial at the same variance.

Build one with
[`pig1_distrib()`](https://statmodels7.github.io/distributions7/reference/pig1_distrib.md),
which supplies the two link functions. This page documents the raw S7
constructor, which validates none of the relationships between its
properties.

## Usage

``` r
Pig1Distrib(
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

An S7 object of class `Pig1Distrib`, inheriting from `discrete_distrib`
and from `distrib`. For an object built by
[`pig1_distrib()`](https://statmodels7.github.io/distributions7/reference/pig1_distrib.md)
the properties hold `"poisson-inverse gaussian"`, `"univariate"`,
`c(0, Inf)`, `c("mu", "sigma")`, the interpretations
`c(mu = "mean", sigma = "dispersion")`, `2`, and two domains \\(0,
\infty)\\.

## Methods

Registered in this file, the middle five reading one compiled kernel:
[`distrib_pdf()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.Pig1Distrib.md),
[`distrib_gradient()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.Pig1Distrib.md),
[`distrib_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.Pig1Distrib.md),
[`distrib_deriv3()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.Pig1Distrib.md),
[`distrib_deriv4()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.Pig1Distrib.md),
[`distrib_rng()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.Pig1Distrib.md).

Registered elsewhere: the four moments in `moments.R` and the data-based
starting value in `starting_values.R`.

The **distribution function** and the **quantile** come from
[`discrete_distrib()`](https://statmodels7.github.io/distributions7/reference/discrete_distrib.md),
where both are exact sums over the support. The **expected information**
has no closed form and goes through
[`expected_derivative_methods()`](https://statmodels7.github.io/distributions7/reference/expected_derivative_methods.md).

## See also

[`pig1_distrib()`](https://statmodels7.github.io/distributions7/reference/pig1_distrib.md)
to build one;
[`pig2_distrib()`](https://statmodels7.github.io/distributions7/reference/pig2_distrib.md)
for the parametrization whose two parameters are orthogonal;
[`negbin2_distrib()`](https://statmodels7.github.io/distributions7/reference/negbin2_distrib.md)
for the other overdispersed count family;
[`pig_hd_block()`](https://statmodels7.github.io/distributions7/reference/pig_hd_block.md)
for the kernel all five derivative methods read.

## Examples

``` r
d <- pig1_distrib()
S7::S7_inherits(d, discrete_distrib)
#> [1] TRUE

d@params
#> [1] "mu"    "sigma"
d@params_interpretation
#>           mu        sigma 
#>       "mean" "dispersion" 

# The mean is a parameter and the variance follows from both.
th <- list(mu = 3, sigma = 0.8)
c(mean = mean(d, th), variance = variance(d, th),
  formula = 3 + 0.8 * 3^2)
#>     mean variance  formula 
#>      3.0     10.2     10.2 

# The mass sums to one over the support.
sum(distrib_pdf(d, 0:300, th))
#> [1] 1
```
