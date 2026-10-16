# Default Log-CDF Gradient for Continuous Distributions

The fallback for a continuous family that registers no closed form: one
central difference of
[`distrib_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.md)
per parameter, through
[`numerical_cdf_deriv()`](https://statmodels7.github.io/distributions7/reference/numerical_cdf_deriv.md),
put on the requested tail and scale by
[`cdf_tail_scale()`](https://statmodels7.github.io/distributions7/reference/cdf_tail_scale.md).
The step is \\\varepsilon^{1/3}\\ relative, about \\6.1\times10^{-6}\\.

## Arguments

- distrib:

  A `continuous_distrib` object.

- q:

  A numeric vector of quantiles.

- theta:

  A named list of parameters on the parameter scale. Components may be
  vectors, in which case the step is chosen elementwise.

- lower.tail:

  A single logical, `TRUE` by default; `TRUE` selects the lower tail and
  `FALSE` the upper tail.

- log:

  A single logical, `TRUE` by default; `TRUE` selects the derivatives of
  the log probability and `FALSE` those of the probability. Far into a
  tail the probability underflows and the result is `-Inf` or `NaN`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A named list of numeric vectors, one per parameter, each the length of
`q` recycled against `theta`.

## Details

No shipped univariate family reaches this method. Every one registers
its own gradient: a closed form or, for the families whose derivative in
a shape parameter has no elementary form, the compiled quadrature of
[`compiled_cdf()`](https://statmodels7.github.io/distributions7/reference/compiled_cdf.md).
The method serves continuous families defined by the user and the
wrappers of continuous families.

## See also

[`numerical_cdf_deriv()`](https://statmodels7.github.io/distributions7/reference/numerical_cdf_deriv.md)
for the stencil and its step;
[`distrib_grad_cdf.discrete_distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_cdf.discrete_distrib.md),
which is exact;
[`distrib_hess_cdf.continuous_distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_cdf.continuous_distrib.md)
for the second order.

## Examples

``` r
# A truncated gamma registers no gradient of its own, so it reaches this
# method and the cdf is differenced.
d <- truncated(gamma2_distrib(), upper = 6)
distrib_grad_cdf(d, c(1, 2, 4), list(mu = 2, sigma2 = 1))
#> $mu
#> [1] -2.38235834 -0.74685542 -0.05165589
#> 
#> $sigma2
#> [1]  1.13004209  0.06780059 -0.05727458
#> 
```
