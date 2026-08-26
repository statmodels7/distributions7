# Default Log-CDF Gradient for Continuous Distributions

The fallback for a continuous family that registers no closed form: one
central difference of
[`distrib_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.md)
per parameter, through
[`numerical_cdf_deriv()`](https://statmodels7.github.io/distributions7/reference/numerical_cdf_deriv.md),
put on the requested tail and scale by
[`cdf_tail_scale()`](https://statmodels7.github.io/distributions7/reference/cdf_tail_scale.md).
The step is \\\varepsilon^{1/3}\\ relative, about \\6.1\times10^{-6}\\,
and the accuracy measured against a family's own closed form is
\\6.1\times10^{-11}\\ relative.

## Arguments

- distrib:

  A `continuous_distrib` object.

- q:

  A numeric vector of quantiles.

- theta:

  A named list of parameters on the parameter scale. Components may be
  vectors, in which case the step is chosen elementwise.

- lower.tail:

  Is the lower tail wanted? A single logical, `TRUE` by default.

- log:

  Are derivatives of the log probability wanted? A single logical,
  `TRUE` by default. Far into a tail the probability underflows and the
  result is `-Inf` or `NaN`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A named list of numeric vectors, one per parameter, each the length of
`q` recycled against `theta`.

## Details

Eight of the 42 univariate families reach this method: beta1, beta2,
chisq, gamma1, gamma2, gengamma1 and the two von Mises. In the first six
the derivative wanted is that of an incomplete gamma or beta function
with respect to its shape, which is hypergeometric and has no elementary
form; in the last two the distribution function is itself a quadrature.
Every other continuous family registers a formula.

## See also

[`numerical_cdf_deriv()`](https://statmodels7.github.io/distributions7/reference/numerical_cdf_deriv.md)
for the stencil and its step;
[`distrib_grad_cdf.discrete_distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_cdf.discrete_distrib.md),
which is exact;
[`distrib_hess_cdf.continuous_distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_cdf.continuous_distrib.md)
for the second order.

## Examples

``` r
# A gamma reaches this method: its cdf derivative in the shape has no
# elementary form, so the cdf is differenced instead.
d <- gamma2_distrib()
distrib_grad_cdf(d, c(1, 2, 4), list(mu = 2, sigma2 = 1))
#> $mu
#> [1] -2.38389173 -0.74838881 -0.05318929
#> 
#> $sigma2
#> [1]  1.12093390  0.05869240 -0.06638277
#> 
```
