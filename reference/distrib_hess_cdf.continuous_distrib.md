# Default Log-CDF Hessian for Continuous Distributions

The fallback for a continuous family that registers no closed second
derivative:
[`numerical_cdf_deriv()`](https://statmodels7.github.io/distributions7/reference/numerical_cdf_deriv.md)
at order 2, with the first-order part taken from
[`distrib_grad_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_cdf.md),
so a closed gradient the family does register is used here too. The step
is \\\varepsilon^{1/4}\\ relative, about \\1.2\times10^{-4}\\, and the
accuracy measured against a family's own closed form is
\\1.7\times10^{-7}\\ relative.

## Arguments

- distrib:

  A `continuous_distrib` object.

- q:

  A numeric vector of quantiles.

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
asked for.

## See also

[`distrib_grad_cdf.continuous_distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_cdf.continuous_distrib.md)
for the first order;
[`numerical_cdf_deriv()`](https://statmodels7.github.io/distributions7/reference/numerical_cdf_deriv.md)
for the stencils;
[`cdf_tail_scale()`](https://statmodels7.github.io/distributions7/reference/cdf_tail_scale.md).

## Examples

``` r
# A gamma reaches this method at both orders.
d <- gamma2_distrib()
distrib_hess_cdf(d, c(1, 2), list(mu = 2, sigma2 = 1))
#> $mu_mu
#> [1] -2.1633205 -0.8750296
#> 
#> $sigma2_sigma2
#> [1] -1.94894361 -0.03301618
#> 
#> $mu_sigma2
#> [1] 1.9255857 0.3739333
#> 
```
