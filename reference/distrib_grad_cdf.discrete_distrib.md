# Log-CDF Gradient for Discrete Distributions

Computes the gradient exactly, without differencing. The partial
expectation of the score is a finite sum over the support up to \\q\\,
\\\partial^I F(q) = \sum\_{y \le q} f(y)\\\partial^I f/f\\(y)\\,
evaluated by
[`discrete_cdf_deriv()`](https://statmodels7.github.io/distributions7/reference/discrete_cdf_deriv.md)
and put on the requested tail by
[`cdf_tail_scale()`](https://statmodels7.github.io/distributions7/reference/cdf_tail_scale.md).

## Arguments

- distrib:

  A `discrete_distrib` object.

- q:

  A numeric vector of quantiles. Values below the support give a
  derivative of zero, because the sum is empty.

- theta:

  A named list of parameters on the parameter scale.

- lower.tail:

  A single logical, `TRUE` by default; `TRUE` selects the lower tail and
  `FALSE` the upper tail.

- log:

  A single logical, `TRUE` by default; `TRUE` selects the derivatives of
  the log probability and `FALSE` those of the probability.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A named list of numeric vectors, one per parameter, in which each vector
has the length of `q` recycled against `theta`.

## Details

The sum is finite because the discrete class requires a finite lower
bound on the support. The cost therefore grows with the largest
requested quantile and not with the number of parameters, and each
quantile is summed separately.

Some of the shipped discrete families use this method; the others
register a formula, several of which are exact simplifications of this
same sum. A Poisson's, for instance, telescopes to \\\partial
F(k)/\partial\mu = -f(k)\\.

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
