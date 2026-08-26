# Log-CDF Gradient for Discrete Distributions

Exact, and nothing is differenced. The partial expectation of the score
is a finite sum over the support up to \\q\\, \\\partial^I F(q) =
\sum\_{y \le q} f(y)\\\partial^I f/f\\(y)\\, evaluated by
[`discrete_cdf_deriv()`](https://statmodels7.github.io/distributions7/reference/discrete_cdf_deriv.md)
and put on the requested tail by
[`cdf_tail_scale()`](https://statmodels7.github.io/distributions7/reference/cdf_tail_scale.md).

## Arguments

- distrib:

  A `discrete_distrib` object.

- q:

  A numeric vector of quantiles. Values below the support give a
  derivative of zero, the sum being empty.

- theta:

  A named list of parameters on the parameter scale.

- lower.tail:

  Is the lower tail wanted? A single logical, `TRUE` by default.

- log:

  Are derivatives of the log probability wanted? A single logical,
  `TRUE` by default.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A named list of numeric vectors, one per parameter, each the length of
`q` recycled against `theta`.

## Details

The sum is finite because the discrete class requires a finite lower
bound on the support. The cost therefore grows with the largest quantile
asked for and not with the number of parameters: on a Poisson of mean 30
at 200 quantiles the whole gradient takes 0.14 ms.

Six of the shipped discrete families use this method; the rest register
a formula, several of which are exact simplifications of this same sum.
A Poisson's, for instance, telescopes to \\\partial F(k)/\partial\mu =
-f(k)\\.

## See also

[`discrete_cdf_deriv()`](https://statmodels7.github.io/distributions7/reference/discrete_cdf_deriv.md)
for the sum;
[`distrib_grad_cdf.PoissonDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_cdf.PoissonDistrib.md)
for a family that simplifies it;
[`distrib_hess_cdf.discrete_distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_cdf.discrete_distrib.md)
for the second order.

## Examples

``` r
# A beta-binomial reaches this method, and the sum is exact.
d <- betabinom1_distrib(size = 10)
distrib_grad_cdf(d, c(2, 5, 8), list(mu = 0.3, sigma = 0.5))
#> $mu
#> [1] -2.8015558 -1.3321275 -0.4226089
#> 
#> $sigma
#> [1]  0.2456977 -0.1557802 -0.1469920
#> 
```
