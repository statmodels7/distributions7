# Poisson-Inverse Gaussian Distribution Class, Orthogonal Parametrization

This is the S7 class of the Poisson-inverse Gaussian in the
parametrization whose two parameters are orthogonal, gamlss's `PIG2`. It
is the same law as
[Pig1Distrib](https://statmodels7.github.io/distributions7/reference/Pig1Distrib.md);
the mean stays \\\mu\\ and the dispersion is replaced by \\\alpha\\, the
argument the mass function's Bessel function is evaluated at.

Orthogonal means the expected information is diagonal: at \\\mu = 3\\
and \\\sigma = 0.8\\ its mixed entry is zero to rounding error here,
against 0.0368 in the mean-dispersion parametrization.

Build one with
[`pig2_distrib()`](https://statmodels7.github.io/distributions7/reference/pig2_distrib.md),
which supplies the two link functions. This page documents the raw S7
constructor, which takes the parent's properties; the validator of
`distrib` checks them, and the constructor adds no check of its own.

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

The following methods are registered on this class, and each is linked
to the page that documents it:

- the distribution:
  [`distrib_pdf()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.Pig2Distrib.md)
  and
  [`distrib_rng()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.Pig2Distrib.md);

- the derivatives in the parameters:
  [`distrib_gradient()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.Pig2Distrib.md),
  [`distrib_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.Pig2Distrib.md),
  [`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.Pig2Distrib.md),
  [`distrib_deriv3()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.Pig2Distrib.md)
  and
  [`distrib_deriv4()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.Pig2Distrib.md);

- the derivatives of the expected information:
  [`distrib_dexpected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.Pig2Distrib.md)
  and
  [`distrib_d2expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.Pig2Distrib.md);

- the moments:
  [`mean()`](https://statmodels7.github.io/distributions7/reference/mean.Pig2Distrib.md),
  [`variance()`](https://statmodels7.github.io/distributions7/reference/variance.Pig2Distrib.md),
  [`skewness()`](https://statmodels7.github.io/distributions7/reference/skewness.Pig2Distrib.md)
  and
  [`kurtosis()`](https://statmodels7.github.io/distributions7/reference/kurtosis.Pig2Distrib.md);

- the fitting helpers:
  [`distrib_start()`](https://statmodels7.github.io/distributions7/reference/distrib_start.Pig2Distrib.md)
  and
  [`expected_hessian_costly()`](https://statmodels7.github.io/distributions7/reference/expected_hessian_costly.md).

Every other generic is inherited from
[`discrete_distrib()`](https://statmodels7.github.io/distributions7/reference/discrete_distrib.md).
The distribution function and the quantile are the inherited ones,
documented at
[`distrib_cdf.discrete_distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.discrete_distrib.md)
and
[`distrib_quantile.discrete_distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.discrete_distrib.md).

## See also

[`pig2_distrib()`](https://statmodels7.github.io/distributions7/reference/pig2_distrib.md)
to build one;
[`pig1_distrib()`](https://statmodels7.github.io/distributions7/reference/pig1_distrib.md)
for the mean-dispersion parametrization;
[`pig2_sigma()`](https://statmodels7.github.io/distributions7/reference/pig2_sigma.md)
for the map between them.

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
