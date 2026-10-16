# Default Log-CDF Hessian for Continuous Distributions

The fallback for a continuous family that registers no closed second
derivative:
[`numerical_cdf_deriv()`](https://statmodels7.github.io/distributions7/reference/numerical_cdf_deriv.md)
at order 2, with the first-order part taken from
[`distrib_grad_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_cdf.md),
so a closed gradient the family does register is used here too. The step
is \\\varepsilon^{1/4}\\ relative, about \\1.2\times10^{-4}\\.

## Arguments

- distrib:

  A `continuous_distrib` object.

- q:

  A numeric vector of quantiles.

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

A named list of numeric vectors keyed as
[`hess_names()`](https://statmodels7.github.io/distributions7/reference/hess_names.md),
each the length of `q` recycled against `theta`. The gradient is not
returned alongside.

## Details

A diagonal component is the three-point second difference and an
off-diagonal one the four-point mixed stencil, which differences two
different variables and is therefore a single stencil. The log-scale
correction \\\partial^2 P/P - (\partial P/P)^2\\ is applied afterwards
by
[`cdf_tail_scale()`](https://statmodels7.github.io/distributions7/reference/cdf_tail_scale.md),
which is why the gradient is fetched even when only the Hessian was
requested.

## See also

[`distrib_grad_cdf.continuous_distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_cdf.continuous_distrib.md)
for the first order;
[`numerical_cdf_deriv()`](https://statmodels7.github.io/distributions7/reference/numerical_cdf_deriv.md)
for the stencils;
[`cdf_tail_scale()`](https://statmodels7.github.io/distributions7/reference/cdf_tail_scale.md).

## Examples

``` r
# A truncated gamma reaches this method at both orders.
d <- truncated(gamma2_distrib(), upper = 6)
distrib_hess_cdf(d, c(1, 2), list(mu = 2, sigma2 = 1))
#> $mu_mu
#> [1] -2.1595748 -0.8712839
#> 
#> $sigma2_sigma2
#> [1] -1.93036989 -0.01444246
#> 
#> $mu_sigma2
#> [1] 1.930828 0.379176
#> 
```
