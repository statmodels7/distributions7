# Exponential Log-CDF Derivatives

Closed form at every order from one to four, from the survival function
\\S = e^{-q/\mu}\\. Its logarithm is \\L = -q/\mu\\, whose partial
derivatives are \\\partial^j L/\partial\mu^j = -q(-1)^j j!/\mu^{j+1}\\,
and
[`surv_cdf_deriv_k()`](https://statmodels7.github.io/distributions7/reference/surv_cdf_deriv_k.md)
turns those into the derivatives of \\F\\.

## Arguments

- distrib:

  An `ExponentialDistrib` object, from
  [`exponential_distrib()`](https://statmodels7.github.io/distributions7/reference/exponential_distrib.md).

- q:

  A numeric vector of quantiles. Values at or below zero give
  derivatives of exactly zero.

- theta:

  A named list with one component, `mu` (positive), a numeric vector of
  length 1 or `n`.

- lower.tail:

  Is the lower tail wanted? A single logical, `TRUE` by default.

- log:

  Are derivatives of the log probability wanted? A single logical,
  `TRUE` by default.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A named list with one numeric vector per component of the order the
generic asked for, which for a one-parameter family is one at every
order.

## Details

Against a product stencil on the same cdf: \\2.7\times10^{-11}\\ at
order 1 and \\2.7\times10^{-5}\\ at order 4. Below the support every
derivative is exactly zero, the mask in
[`surv_cdf_deriv_k()`](https://statmodels7.github.io/distributions7/reference/surv_cdf_deriv_k.md)
suppressing the finite value \\L\\ would otherwise give there.

## Notation

\\\mu \> 0\\ is the mean, \\F\\ the distribution function and \\S = 1 -
F\\ the survival function.

## See also

[`surv_cdf_deriv_k()`](https://statmodels7.github.io/distributions7/reference/surv_cdf_deriv_k.md)
for the identity;
[`distrib_grad_cdf.Weibull1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_cdf.Weibull1Distrib.md)
and
[`distrib_grad_cdf.GPDDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_cdf.GPDDistrib.md),
the two families that contain this one;
[`exponential_distrib()`](https://statmodels7.github.io/distributions7/reference/exponential_distrib.md).

## Examples

``` r
d <- exponential_distrib()
q <- c(0.5, 2, 5)

# Against a central difference of the cdf, which shares no arithmetic.
fd <- numerical_cdf_deriv(d, q, list(mu = 3), order = 1)
max(abs(distrib_grad_cdf(d, q, list(mu = 3), log = FALSE)$mu / fd$mu - 1))
#> [1] 2.722822e-11

# Exactly zero below the support.
distrib_grad_cdf(d, c(-1, 0.5), list(mu = 3), log = FALSE)$mu
#> [1]  0.00000000 -0.04702676
```
