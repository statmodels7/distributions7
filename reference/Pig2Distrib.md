# Poisson-Inverse Gaussian Distribution Class, Orthogonal Parametrization

The S7 class of the Poisson-inverse Gaussian in the parametrization
whose two parameters are orthogonal, gamlss's `PIG2`. It is the same law
as
[Pig1Distrib](https://statmodels7.github.io/distributions7/reference/Pig1Distrib.md);
the mean stays \\\mu\\ and the dispersion is replaced by \\\alpha\\, the
argument the mass function's Bessel function is evaluated at.

Orthogonal means the expected information is diagonal. Measured at \\\mu
= 3\\, its mixed entry summed over the support is \\-8.8\times10^{-15}\\
here against 7.39 in the mean-dispersion parametrization.

Build one with
[`pig2_distrib()`](https://statmodels7.github.io/distributions7/reference/pig2_distrib.md),
which supplies the two link functions. This page documents the raw S7
constructor, which validates none of the relationships between its
properties.

## Usage

``` r
Pig2Distrib(
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

An S7 object of class `Pig2Distrib`, inheriting from `discrete_distrib`
and from `distrib`. For an object built by
[`pig2_distrib()`](https://statmodels7.github.io/distributions7/reference/pig2_distrib.md)
the properties hold `"poisson-inverse gaussian (orthogonal)"`,
`"univariate"`, `c(0, Inf)`, `c("mu", "alpha")`, the interpretations
`c(mu = "mean", alpha = "bessel argument")`, `2`, and two domains \\(0,
\infty)\\.

## Methods

Registered in this file, the middle five reading one compiled kernel:
[`distrib_pdf()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.Pig2Distrib.md),
[`distrib_gradient()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.Pig2Distrib.md),
[`distrib_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.Pig2Distrib.md),
[`distrib_deriv3()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.Pig2Distrib.md),
[`distrib_deriv4()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.Pig2Distrib.md),
[`distrib_rng()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.Pig2Distrib.md).

Registered elsewhere: the four moments in `moments.R` and the data-based
starting value in `starting_values.R`. The distribution function and the
quantile come from
[`discrete_distrib()`](https://statmodels7.github.io/distributions7/reference/discrete_distrib.md),
both exact sums over the support.

## See also

[`pig2_distrib()`](https://statmodels7.github.io/distributions7/reference/pig2_distrib.md)
to build one;
[`pig1_distrib()`](https://statmodels7.github.io/distributions7/reference/pig1_distrib.md)
for the mean-dispersion parametrization;
[`pig2_sigma()`](https://statmodels7.github.io/distributions7/reference/pig2_sigma.md)
for the map between them;
[`pig_hd_block()`](https://statmodels7.github.io/distributions7/reference/pig_hd_block.md)
for the kernel all five derivative methods read.

## Examples

``` r
d <- pig2_distrib()
S7::S7_inherits(d, discrete_distrib)
#> [1] TRUE

d@params
#> [1] "mu"    "alpha"
d@params_interpretation
#>                mu             alpha 
#>            "mean" "bessel argument" 

# The same law as pig1, at the dispersion alpha implies.
al <- 3.010399
d1 <- pig1_distrib()
rbind(pig2 = distrib_pdf(d, 0:5, list(mu = 3, alpha = al)),
      pig1 = distrib_pdf(d1, 0:5,
                         list(mu = 3,
                              sigma = distributions7:::pig2_sigma(3, al))))
#>           [,1]      [,2]      [,3]      [,4]       [,5]       [,6]
#> pig2 0.1719763 0.2142278 0.1777529 0.1289567 0.08968701 0.06196187
#> pig1 0.1719763 0.2142278 0.1777529 0.1289567 0.08968701 0.06196187
```
